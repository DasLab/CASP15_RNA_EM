from casp_rna_em.process_files import prepare_pdbs, concat_all_result_files, center_map_and_native, reduce_df, write_and_sbatch_scoring
from os import system, path
import numpy as np
import pandas as pd
from glob import glob

# TODO
#   for all targets
#   combine output after run -- new script


usalign_exec = '/home/anthony/USalign/USalign'
chimerax_exec = '/home/anthony/chimerax-1.6-rc2023.03.31//bin/ChimeraX'
phenix_location = '/home/anthony/phenix-1.20.1-4487/build/bin/'

cwd = '/home/anthony/Desktop/CASP/script_test'


def run_all_on_target(tar_file, natives, models_folder, name,
                      prepare, run, EM=False, resolution=None, threshold=None, center=False, N=None,
                      base_sbatch=None, sbatch=False):
    if prepare:
        print(f'Preparing models in {models_folder}.')
        if path.isdir(models_folder):
            system(f'rm -r {models_folder}')
        system(f'tar -xf {tar_file}')
        if EM:
            for i, emmap in enumerate(natives):
                for j, native in enumerate(natives[emmap]):
                    system(f'cp {native} {models_folder}/{name}TSex{i}_{j}')
        else:
            for i, native in enumerate(natives):
                system(f'cp {native} {models_folder}/{name}TSexp_{i}')
        pdbs = prepare_pdbs(f'{models_folder}/{name}TS???_?')
    else:
        pdbs = glob(f'{models_folder}/{name}TS???_?.pdb')
    if run:
        if EM:
            for i, emmap in enumerate(natives):
                for j, native in enumerate(natives[emmap]):
                    native = f'{models_folder}/{name}TSex{i}_{j}.pdb'
                    # only run EM comparison once
                    if j == 0:
                        if center:
                            center_map_and_native(emmap, native)
                        print(f'Scoring models (EM) in {models_folder}, {native}, {emmap}.')
                        write_and_sbatch_scoring(pdbs, N, f'run_files/run_{name}', name, base_sbatch, f"result_files/{name}_{i}-{j}", native, usalign_exec=usalign_exec, chimerax_exec=chimerax_exec, EM=EM, emmap=emmap, resolution=resolution[emmap], threshold=threshold[emmap], phenix_location=phenix_location, sbatch=sbatch)
                    else:
                        print(f'Scoring models in {models_folder}, {native}, {emmap}.')
                        write_and_sbatch_scoring(pdbs, N, f'run_files/run_{name}', name, base_sbatch, f"result_files/{name}_{i}-{j}", native, usalign_exec=usalign_exec, EM=False, phenix_location=phenix_location, sbatch=sbatch)
        else:
            for i, native in enumerate(natives):
                native = f'{models_folder}/{name}TSexp_{i}.pdb'
                print(f'Scoring models in {models_folder}, {native}.')
                write_and_sbatch_scoring(pdbs, N, f'run_files/run_{name}', name, base_sbatch, f"result_files/{name}_{i}", native, usalign_exec=usalign_exec, EM=EM, phenix_location=phenix_location, sbatch=sbatch)


###############################################################################
# R1107 ~9.2s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1107.tar.gz',
                  natives=[f'{cwd}/D_1292119758_model-annotate_P1human.pdb'],
                  models_folder=f'{cwd}/R1107',
                  name="R1107",
                  prepare=False,
                  run=True,
                  base_sbatch=f'{cwd}/example.sbatch',
                  sbatch=False,
                  N=2)
