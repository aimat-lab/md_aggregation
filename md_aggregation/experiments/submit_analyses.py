import os
import re
import json
import pathlib
import rich_click as click
import tempfile
import jinja2 as j2
from rich.pretty import pprint
from typing import List, Dict, Any
from pycomex.functional.experiment import Experiment

# This is the path to the folder in which the current file resides and which 
# will subsequently be used as the anchor to find the other relevant locations in 
# the filesystem.
PATH = pathlib.Path(__file__).parent.absolute()


def sanitize_molecule_fields(input_str):
    """
    Removes all occurrences of the pattern starting with ,'molecule' and ending
    at the next top-level comma (not within brackets/braces).
    """
    pattern = r",\s*'molecule'\s*:\s*[^,}]*"

    def replacer(match):
        s = match.group(0)
        # Find where to cut: keep last character if it's } or ,
        for i in range(len(s)):
            if s[i] in ',}':
                return ''  # don't include any of the matched portion
        return ''  # fallback: remove entire matched portion

    # Apply all matches
    return re.sub(pattern, replacer, input_str)


@click.command()
@click.option('-c', '--config', type=click.STRING, default='haicore_1gpu', show_default=True,
              help='AutoSlurm configuration under wich to submit the jobs.')
@click.option('-p', '--bash-path', type=click.Path(), show_default=True,
              default=os.path.join(tempfile.gettempdir(), 'submit_analyses.sh'),
              help='Path to the folder in which the results are stored.')
@click.option('-r', '--cutoff', type=click.FLOAT, default=3.0, show_default=True,
              help='The cutoff distance in angstroms for the clustering algorithm.')
@click.option('-s', '--submit', is_flag=True, default=False,
              help='Submit the job to the SLURM queue after contructing the bash script.')
@click.option('-e', '--environment', type=click.STRING, default='md_aggregation', show_default=True,
              help='The conda environment to be activated to execute the analysis script.')
def main(config: str, bash_path: str, cutoff: float, submit: bool, environment: str):
    """
    This script will iterate over all the individual experiment archives in the "results" folder 
    and if there are any experiments which are derived from the "aggregation_simulation.py" base 
    experiment, which have terminated but not been analyzed yet, this script will submit the 
    analysis job to the SLURM queue using AutoSlurm.
    """
    
    results_path: str = PATH / 'results'

    template_env = j2.Environment(
        loader=j2.FileSystemLoader(PATH),
    )
    template: j2.Template = template_env.get_template('submit_analyses.sh.j2')

    click.echo('Scanning results folder: {}'.format(results_path))
    
    # This list will be used to store the paths to all the arichive folders which need to be processed
    # with the analysis script later on. This list will be populated as the individual experiment folders 
    # are scanned.
    relevant_archive_folders: List[str] = []
    
    # The results folder itself contains a series of subfolders which are named the same as the experiment 
    # modules. Each of these folders then contains a series of subfolders which represent the individual 
    # runs of the corresponding experiment.
    experiment_folders: List[str] = [
        path 
        for file in os.listdir(results_path) 
        if os.path.isdir(path := results_path / file) and 
            path.name.startswith('aggregation_simulation_')
    ]
    for folder_path in experiment_folders:
        
        click.echo(f' * experiment: "{folder_path}"')
        archive_folders = [
            path.absolute()
            for path in pathlib.Path(folder_path).iterdir()
            if path.is_dir() 
                and os.path.exists(path / 'trajectory.dcd')
                and os.path.exists(path / 'topology.pdb')
        ]
        
        for archive_path in archive_folders:
            
            # Check if the analysis has already been performed
            if os.path.exists(archive_path / 'per_frame_clusters.csv'):
                click.echo(f'   ! analysis already done for {archive_path.name}.')
                continue
            
            # Check if the job is still running
            if os.path.exists(archive_path / 'job_running.txt'):
                click.echo(f'   ! not done for {archive_path.name}.')
                continue
            
            # If the analysis has not been performed and the job is not running, add it to the list
            click.echo(f'   - adding {archive_path.name}')
            relevant_archive_folders.append(archive_path)

    # This will be a list of dictionaries where each dictionary will contain all the necessary information
    # about one of the relevant archive folders to be processed. This information will include 
    relevant_infos: List[dict] = []
    
    for path in relevant_archive_folders:
        
        metadata_path = path / 'experiment_meta.json'
        with open(metadata_path, 'r') as file:
            content = file.read()
            metadata = json.loads(content)
        
        molecule_config_raw = sanitize_molecule_fields(metadata['parameters']['MOLECULE_CONFIGURATION']['value'])
        molecule_config = eval(molecule_config_raw)
        
        info = {
            'path': path,
            'topology_path': path / 'topology.pdb',
            'trajectory_path': path / 'trajectory.dcd',
            'solvent': molecule_config['water']['residue_name'],
            'exclude': " ".join([
                conf['residue_name']
                for name, conf in molecule_config.items()
                if name != 'water' and name != 'ligand'
            ]),
            'end': 5000,
        }
        relevant_infos.append(info)
        
    # Finally, using the list of the relevant archive folders we use the given jinja template to 
    # construct the bash script that defines the correct commands to schedule all the analysis scripts 
    # to the SLURM queue.
    content = template.render({
        'environment': environment,
        'project_path': str(PATH.parent.parent.absolute()),
        'config': config,
        'archive_folders': relevant_archive_folders,
        'infos': relevant_infos,
        'cutoff': cutoff,
    })
    
    with open(bash_path, mode='w') as file:
        file.write(content)

    click.echo(f'DONE - bash script written @ {bash_path}!')
    
    if submit:
        click.echo('Submitting the job to the SLURM queue...')
        os.system(f'sbatch {bash_path}')
    else:
        click.echo('skipping submission, youll have to do it manually!')


if __name__ == "__main__":

    # Call the main function
    main()