"""
The base experiment file which implements the aggregation MD simulation using OpenMM.
"""
import os
import sys
import time
import datetime
import logging
from typing import Dict, List, Tuple, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pycomex.functional.experiment import Experiment
from pycomex.utils import folder_path, file_namespace

from rdkit import Chem
import openmm
import parmed
from openmm import Platform
from parmed.openmm import MdcrdReporter
from openff.toolkit.typing.engines.smirnoff import ForceField as OFF_ForceField
from openff.toolkit.topology import Molecule, Topology
from openff.interchange import Interchange
from openff.interchange.components._packmol import RHOMBIC_DODECAHEDRON, pack_box
from openff.units import unit as off_unit


# == SIMULATION PARAMETERS ==

BOX_VECTORS: np.ndarray = np.array([
    [1.0, 0.0, 0.0],
    [0.5, np.sqrt(3)/2, 0.0],
    [0.5, 1/(2*np.sqrt(3)), np.sqrt(2/3)]
])
BOX_SIZE: float = 180 # in angstroms

MOLECULE_CONFIGURATION: Dict[str, dict] = {
    'ligand': {
        'num': 11,
        'residue_name': 'LIG',
        'smiles': 'C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O',
        'use_rdkit': True,
    },
    'dmso': {
        'num': 320,
        'residue_name': 'DMS',
        'smiles': 'CS(=O)C',
        'use_rdkit': True,
    },
    'water': {
        'num': 35000,
        'residue_name': 'HOH',
        'smiles': '[O]([H])[H]',
        'use_rdkit': False,
    },
    'na': {
        'num': 32,
        'residue_name': 'NA',
        'smiles': '[Na+]',
        'use_rdkit': False,
    },
    'cl': {
        'num': 32,
        'residue_name': 'CL',
        'smiles': '[Cl-]',
        'use_rdkit': False,
    }
}

NON_BONDED_CUTOFF: float = 10.0  # in angstroms

TEMPERATURE: float = 300.0  # in Kelvin
PRESSURE: float = 1.0  # in bar

MINIMIZATION_ITERATIONS: int = 5000

RAMP_STEPS: int = 50_000
EQUILIBRATION_STEPS: int = 100_000

PRODUCTION_STEPS: int = 50_000_000

FRAME_STRIDE: int = 10_000

DEVICE_INDEX: int = 0

__PREFIX__ = 'org'

# == EXPERIMENT PARAMETERS ==

__DEBUG__ = False
__TESTING__ = False

experiment = Experiment(
    base_path=folder_path(__file__),
    namespace=file_namespace(__file__),
    glob=globals(),
)


@experiment.hook('generate_molecule', default=True, replace=False)
def generate_molecule(e: Experiment,
                      smiles: str,
                      residue_name: str,
                      use_rdkit: bool = False
                      ) -> Molecule:
    
    # For complex molecules we first use RDKIT to process the SMILES
    if use_rdkit:
        
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        mol = Chem.AddHs(mol)
        
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS)
        molecule = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
        molecule.generate_conformers(n_conformers=1)
    
    # For simpler molecules we can just load them from SMILES directly
    else:
        
        molecule = Molecule.from_smiles(smiles)
        
    # Set the residue name - this will later be important for the visualization
    for atom in molecule.atoms:
        atom.metadata["residue_name"] = residue_name

    return molecule


class LoggerWriter:
    
    def __init__(self, logger, level=logging.INFO):
        self.logger = logger
        self.level = level
        self._buffer = ""

    def write(self, message):
        # Buffer partial lines until newline
        self._buffer += message
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            if line.strip():  # avoid logging empty lines
                self.logger.log(self.level, line)

    def flush(self):
        if self._buffer.strip():
            self.logger.log(self.level, self._buffer.strip())
        self._buffer = ""


@experiment
def experiment(e: Experiment):
    
    e.log('starting experiment...')
    
    if e.__TESTING__:
        e.log('running in testing mode...')
        e.MINIMIZATION_ITERATIONS = 100
        e.RAMP_STEPS = 10
        e.EQUILIBRATION_STEPS = 100
        e.PRODUCTION_STEPS = 1001
        e.FRAME_STRIDE = 10
    
    # ~ Platform information
    # Here we want to make sure that the OpenMM simulation will actually use the CUDA acceleration
    platform = Platform.getPlatform(0)
    e.log(f"Default OpenMM platform: {platform.getName()}")

    try:
        platform = Platform.getPlatformByName('CUDA')
        #platform.setPropertyValue()
        e.log(f"Using platform: {platform.getName()} - Device Index: {e.DEVICE_INDEX}")
    except Exception as exc:
        e.log(str(exc))
        e.log("CUDA platform not found!")

    # ~ Setting up the molecules for the simulation
    e.log('Setting up the MOLECULES...')
    
    for key, data in e.MOLECULE_CONFIGURATION.items():
        
        molecule: Molecule = e.apply_hook(
            'generate_molecule',
            smiles=data['smiles'],
            residue_name=data['residue_name'],
            use_rdkit=data.get('use_rdkit', False)
        )
        data['molecule'] = molecule
        
        num_atoms = molecule.n_atoms
        num_bonds = molecule.n_bonds
        e.log(f' * {key} - {num_atoms} atoms - {num_bonds} bonds - {data["num"]}x')

    # ~ Setting up the topology
    # Using the "Pack_box" method which uses PACKMOL in the background to place the required number 
    # of molecule in the given simulation box.
    e.log('Setting up the TOPOLOGY...')
    
    molecules = [data['molecule'] for data in e.MOLECULE_CONFIGURATION.values()]
    nums = [data['num'] for data in e.MOLECULE_CONFIGURATION.values()]
    
    topology = pack_box(
        molecules=molecules,
        number_of_copies=nums,
        box_vectors=np.eye(3) * 180 * off_unit.angstrom,
    )

    e.log('Saving packed topology...')
    topology.to_file(
        os.path.join(e.path, "packed.pdb"),
    )
    
    # ~ Setting up the Interchange
    # The Interchange object is the main object that will be used for the simulation.
    e.log('Setting up the INTERCHANGE...')
    sage = OFF_ForceField("openff-2.1.0.offxml")
    interchange: Interchange = Interchange.from_smirnoff(force_field=sage, topology=topology)

    interchange["vdW"].cutoff = e.NON_BONDED_CUTOFF * off_unit.angstrom
    interchange["Electrostatics"].cutoff = e.NON_BONDED_CUTOFF * off_unit.angstrom
    
    e.log('saving the topology file...')
    topology_path = os.path.join(e.path, 'topology.pdb')
    #interchange.to_pdb(topology_path)

    # ~ Equilibrium Simulation
    # The first part of the simulation consists of a minimization of the energy of the packed system and then 
    # a few ps of equilibration using the NVT ensemble (constant volume simulation).
    
    e.log('Setting up the simulation...')
    integrator = openmm.LangevinIntegrator(
        0 * openmm.unit.kelvin,
        1 / openmm.unit.picosecond,
        2 * openmm.unit.femtoseconds,
    )

    simulation = interchange.to_openmm_simulation(
        combine_nonbonded_forces=False,
        integrator=integrator,
    )
    
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(topology_path, 'w') as file:
        openmm.app.PDBFile.writeFile(simulation.topology, positions, file)
    
    e.log('minimizing energy...')
    # https://github.com/openmm/openmm/issues/3736#issuecomment-1217250635
    time_start = time.time()
    simulation.minimizeEnergy(maxIterations=e.MINIMIZATION_ITERATIONS)
    duration = time.time() - time_start
    e['minimization_duration'] = duration
    e.log(f' * minimization done after {duration:.2f} seconds')
    
    e.log('starting temperature ramp...')
    time_start = time.time()
    
    for i in range(e.RAMP_STEPS):
        
        simulation.step(1)
        
        temperature = (i / e.RAMP_STEPS) * 500
        integrator.setTemperature(temperature * openmm.unit.kelvin)
        if i % 100 == 0:
            e.log(f' * ramping temperature to {temperature:.2f} K')

    for i in range(e.RAMP_STEPS):
        
        simulation.step(1)
        
        temperature = 500 - (i / e.RAMP_STEPS) * (500 - e.TEMPERATURE)
        integrator.setTemperature(temperature * openmm.unit.kelvin)
        if i % 100 == 0:
            e.log(f' * ramping temperature to {temperature:.2f} K')
    
    duration = time.time() - time_start
    e['ramp_duration'] = duration
    e.log(f' * ramp done after {duration:.2f} seconds')
    
    e.log('Starting equilibration simulation...')
    barostat = openmm.MonteCarloBarostat(
        e.PRESSURE * openmm.unit.bar,
        e.TEMPERATURE * openmm.unit.kelvin,
        50,
    )
    simulation.system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)
    
    time_start = time.time()
    simulation.step(e.EQUILIBRATION_STEPS)
    duration = time.time() - time_start
    e['equilibration_duration'] = duration
    e.log(f' * equilibration done after {duration:.2f} seconds')
    
    # Save the equilibrated system
    e.log('Saving equilibrated system...')
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(os.path.join(e.path, "equilibrated.pdb"), 'w') as file:
        openmm.app.PDBFile.writeFile(simulation.topology, state.getPositions(), file)
    
    # ~ Production Simulation
    # After the equilibration, we will run the production simulation using the NPT ensemble (constant pressure simulation).
    e.log('Setting up production simulation...')
    
    production_simulation = simulation
    
    trajectory_path = os.path.join(e.path, "trajectory.pdb")
    pdb_reporter = openmm.app.PDBReporter(trajectory_path, e.FRAME_STRIDE)
    production_simulation.reporters.append(pdb_reporter)
    
    trajectory_dcd_path = os.path.join(e.path, "trajectory.dcd")
    dcd_reporter = openmm.app.DCDReporter(trajectory_dcd_path, e.FRAME_STRIDE)
    production_simulation.reporters.append(dcd_reporter)
    
    state_data_path = os.path.join(e.path, "state_data.csv")
    state_data_reporter = openmm.app.StateDataReporter(
        open(state_data_path, 'w'),
        e.FRAME_STRIDE if not e.__TESTING__ else 1,
        step=True,
        time=True,
        potentialEnergy=True,
        temperature=True,
        density=True,
        volume=True,
    )
    production_simulation.reporters.append(state_data_reporter)

    # Run the production simulation
    e.log('Starting production simulation...')
    step_batch = 1000
    steps_total = 0
    time_start = time.time()
    for _ in range(e.PRODUCTION_STEPS // step_batch):
        production_simulation.step(step_batch)
        
        steps_total += step_batch
        duration = time.time() - time_start
        time_per_step = duration / steps_total
        time_per_ns = time_per_step * 500_000
        time_remaining = (e.PRODUCTION_STEPS - steps_total) * time_per_step
        eta = datetime.datetime.now() + datetime.timedelta(seconds=time_remaining)
        e.log(f' * {steps_total} steps done'
              f' - {duration:.2f} seconds'
              f' - {time_per_step:.5f} seconds/step'
              f' - {time_per_ns/3600:.2f} hrs/ns'
              f' - {time_remaining/3600:.2f} hours remaining'
              f' - ETA: {eta.strftime("%Y-%m-%d %H:%M:%S")}')
        
    duration = time.time() - time_start
    e['production_duration'] = duration
    e.log(f' * production done after {duration:.2f} seconds')
    
    e.log('Saving production simulation...')
    state = production_simulation.context.getState(getPositions=True, getVelocities=True)
    with open(os.path.join(e.path, "production.pdb"), 'w') as file:
        openmm.app.PDBFile.writeFile(production_simulation.topology, state.getPositions(), file)

experiment.run_if_main()
