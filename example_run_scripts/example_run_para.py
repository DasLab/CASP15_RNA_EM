from casp_rna_em.process_files import prepare_pdbs, concat_all_result_files, center_map_and_native, reduce_df, write_and_sbatch_scoring
from os import system, path
import numpy as np
import pandas as pd
from glob import glob

# TODO
#   for all targets
#   combine output after run -- new script


usalign_location = '/home/groups/rhiju/philpham/casp15/us-align/'
chimerax_location = '/home/groups/rhiju/rkretsch/chimerax-1.3/bin/'
phenix_location = '/home/groups/rhiju/rkretsch/phenix/phenix-1.18.2-3874/build/bin/'

cwd = '/scratch/users/rkretsch/230403_casp/data'


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
                        write_and_sbatch_scoring(pdbs, N, f'run_files/run_{name}', name, base_sbatch, f"result_files/{name}_{i}-{j}", native, usalign_location=usalign_location, chimerax_location=chimerax_location, EM=EM, emmap=emmap, resolution=resolution[emmap], threshold=threshold[emmap], phenix_location=phenix_location, sbatch=sbatch)
                    else:
                        print(f'Scoring models in {models_folder}, {native}, {emmap}.')
                        write_and_sbatch_scoring(pdbs, N, f'run_files/run_{name}', name, base_sbatch, f"result_files/{name}_{i}-{j}", native, usalign_location=usalign_location, EM=False, phenix_location=phenix_location, sbatch=sbatch)
        else:
            for i, native in enumerate(natives):
                native = f'{models_folder}/{name}TSexp_{i}.pdb'
                print(f'Scoring models in {models_folder}, {native}.')
                write_and_sbatch_scoring(pdbs, N, f'run_files/run_{name}', name, base_sbatch, f"result_files/{name}_{i}", native, usalign_location=usalign_location, EM=EM, phenix_location=phenix_location, sbatch=sbatch)


###############################################################################
# R1107 ~9.2s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1107.tar.gz',
                  natives=[f'{cwd}/D_1292119758_model-annotate_P1human.pdb'],
                  models_folder=f'{cwd}/R1107',
                  name="R1107",
                  prepare=False,
                  run=False,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=4)

###############################################################################
# R1108 ~9.2s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1108.tar.gz',
                  natives=[f'{cwd}/D_1292119797_model-annotate_P1chimp_A.pdb',
                           f'{cwd}/D_1292119797_model-annotate_P1chimp_A.pdb'],
                  models_folder=f'{cwd}/R1108',
                  name="R1108",
                  prepare=False,
                  run=False,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=8)

###############################################################################
# R1116 ~11.5s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1116.tar.gz',
                  natives=[f'{cwd}/PV_CL_092022_refine_26.pdb'],
                  models_folder=f'{cwd}/R1116',
                  name="R1116",
                  prepare=False,
                  run=False,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=8)


###############################################################################
# R1117 ~8.2s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1117.tar.gz',
                  natives=[f'{cwd}/R1117_structure_31.pdb'],
                  models_folder=f'{cwd}/R1117',
                  name="R1117",
                  prepare=False,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  run=False,
                  N=4)

run_all_on_target(tar_file=f'{cwd}/R1117v2.tar.gz',
                  natives=[f'{cwd}/R1117_structure_31.pdb'],
                  models_folder=f'{cwd}/R1117v2',
                  name="R1117v2",
                  prepare=False,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  run=False,
                  N=4)

###############################################################################
# R1126 ~212s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1126.tar.gz',
                  natives={f'{cwd}/R1126_structure_21.mrc': [f'{cwd}/R1126_structure_21_Traptamer_rsr.pdb']},
                  models_folder=f'{cwd}/R1126',
                  name="R1126",
                  EM=True,
                  prepare=False,
                  run=False,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  resolution={f'{cwd}/R1126_structure_21.mrc': 5.6},
                  threshold={f'{cwd}/R1126_structure_21.mrc': 0.11},
                  N=50)

###############################################################################
# R1128 ~180s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1128.tar.gz',
                  natives={f'{cwd}/R1128_structure_24.mrc': [f'{cwd}/R1128_structure_24_6wj_refined_046.pdb']},
                  models_folder=f'{cwd}/R1128',
                  name="R1128",
                  EM=True,
                  prepare=False,
                  run=False,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  resolution={f'{cwd}/R1128_structure_24.mrc': 5.3},
                  threshold={f'{cwd}/R1128_structure_24.mrc': 0.124},
                  N=50)

###############################################################################
# R1136 ~155s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1136.tar.gz',
                  natives={f'{cwd}/R1136v1_structure_22.mrc': [f'{cwd}/R1136v1_Apta_FRET-Bound_rsr15-withAt2OP1.pdb'],
                           f'{cwd}/R1136v2_Apta-FRET-unbound.mrc': [f'{cwd}/R1136v2_structure_22v2_Apta_FRET-unbound_rsr_deposit.pdb']},
                  models_folder=f'{cwd}/R1136',
                  name="R1136",
                  EM=True,
                  prepare=True,
                  run=True,
                  resolution={f'{cwd}/R1136v1_structure_22.mrc': 4.4,
                              f'{cwd}/R1136v2_Apta-FRET-unbound.mrc': 4.5},
                  threshold={f'{cwd}/R1136v1_structure_22.mrc': 0.11,
                              f'{cwd}/R1136v2_Apta-FRET-unbound.mrc': 0.11},
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=20)

###############################################################################
# R1138 ~460s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1138.tar.gz',
                  natives={f'{cwd}/R1138v1_structure_23.mrc': [f'{cwd}/R1138v1_structure_23v1_6hbc_young.pdb'],
                           f'{cwd}/R1138v2_6HBC-Mature.mrc': [f'{cwd}/R1138v2_structure_23v2_6hbc_mature.pdb']},
                  models_folder=f'{cwd}/R1138',
                  name="R1138",
                  EM=True,
                  prepare=True,
                  run=True,
                  resolution={f'{cwd}/R1138v1_structure_23.mrc': 5.2,
                              f'{cwd}/R1138v2_6HBC-Mature.mrc': 4.9},
                  threshold={f'{cwd}/R1138v1_structure_23.mrc': 0.1,
                              f'{cwd}/R1138v2_6HBC-Mature.mrc': 0.44},
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=50)

###############################################################################
# R1149 ~490s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1149.tar.gz',
                  natives={f'{cwd}/R11149.mrc': [f'{cwd}/R11149_{i}.pdb' for i in range(10)]},
                  models_folder=f'{cwd}/R1149',
                  name="R1149",
                  EM=True,
                  prepare=True,
                  run=True,
                  center=True,
                  resolution={f'{cwd}/R11149.mrc': 4.74},
                  threshold={f'{cwd}/R11149.mrc': 0.345},
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=40)

###############################################################################
# R1156 ~270s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1156.tar.gz',
                  natives={f'{cwd}/R1156_conformation1_resolution5.83.mrc': [f'{cwd}/R1156_conformation1_model{i+1}.pdb' for i in range(10)],
                           f'{cwd}/R1156_conformation2_resolution6.59.mrc': [f'{cwd}/R1156_conformation2_model{i+1}.pdb' for i in range(10)],
                           f'{cwd}/R1156_conformation3_resolution7.48.mrc': [f'{cwd}/R1156_conformation3_model{i+1}.pdb' for i in range(10)],
                           f'{cwd}/R1156_conformation4_resolution7.61.mrc': [f'{cwd}/R1156_conformation4_model{i+1}.pdb' for i in range(10)]},
                  models_folder=f'{cwd}/R1156',
                  name="R1156",
                  EM=True,
                  prepare=True,
                  run=True,
                  center=True,
                  resolution={f'{cwd}/R1156_conformation1_resolution5.83.mrc': 5.83,
                              f'{cwd}/R1156_conformation2_resolution6.59.mrc': 6.59,
                              f'{cwd}/R1156_conformation3_resolution7.48.mrc': 7.48,
                              f'{cwd}/R1156_conformation4_resolution7.61.mrc': 7.61},
                  threshold={f'{cwd}/R1156_conformation1_resolution5.83.mrc': 0.217,
                             f'{cwd}/R1156_conformation2_resolution6.59.mrc': 0.189,
                             f'{cwd}/R1156_conformation3_resolution7.48.mrc': 0.171,
                             f'{cwd}/R1156_conformation4_resolution7.61.mrc': 0.220},
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=20)

###############################################################################
# R1189 ~10.6s/pdb
# EM not run due to poor fits
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1189.tar.gz',
                  natives=[f'{cwd}/RT1189_A-6B.pdb'],
                  models_folder=f'{cwd}/R1189',
                  name="R1189",
                  prepare=True,
                  run=True,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=3)

run_all_on_target(tar_file=f'{cwd}/RT1189.tar.gz',
                  natives=[f'{cwd}/RT1189_A-6B.pdb'],
                  models_folder=f'{cwd}/RT1189',
                  name="RT1189",
                  prepare=True,
                  run=True,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=3)


###############################################################################
# R1190 ~10.6s/pdb
# EM not run due to poor fits
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1190.tar.gz',
                  natives=[f'{cwd}/RT1190_A-4B.pdb'],
                  models_folder=f'{cwd}/R1190',
                  name="R1190",
                  prepare=True,
                  run=True,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=3)

run_all_on_target(tar_file=f'{cwd}/RT1190.tar.gz',
                  natives=[f'{cwd}/RT1190_A-4B.pdb'],
                  models_folder=f'{cwd}/RT1190',
                  name="RT1190",
                  prepare=True,
                  run=True,
                  base_sbatch=f'{cwd}/../CASP15_RNA_EM/example_run_scripts/example.sbatch',
                  sbatch=True,
                  N=2)

