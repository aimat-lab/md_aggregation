from typing import Dict

from pycomex.functional.experiment import Experiment
from pycomex.utils import folder_path, file_namespace


MOLECULE_CONFIGURATION: Dict[str, dict] = {
    'ligand': {
        'num': 11,
        'residue_name': 'LIG',
        'smiles': 'C1(C(=O)NC2SC=CC=2C(=O)N)C=C2CCCCC2=CC=1',
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

experiment = Experiment.extend(
    'aggregation_simulation.py',
    base_path=folder_path(__file__),
    namespace=file_namespace(__file__),
    glob=globals()
)

experiment.run_if_main()