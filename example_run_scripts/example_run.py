from casp_rna_em.process_files import prepare_pdbs, center_map_and_native
from casp_rna_em.process_files import concat_all_result_files, reduce_df
from casp_rna_em.process_files import get_zscore, get_weighted_sum_z
from casp_rna_em.run_metric_programs import score_all
from os import system, path
import pandas as pd
from glob import glob

###############################################################################
# Preparation
###############################################################################
# Create folder with all model tars, native, and native maps in it

###############################################################################
# Results
###############################################################################
# you will run example_run_combine.py to find combined results in all*csv

###############################################################################
# TODO
###############################################################################
# code
#   move in visualization code for all_per_thershold and all_per_residue

# Documentation
#   docstring

# Later:
#   options to run: ddsr, qscore, fsc
#   qscore: fork mapq with parallelization and sigfig fixes
#   add in segmentation code
#   Full?
#       lddt and inf not calculated, not sure which rmsd as well?
#       Z_rna = {"lddt":1/8,"tm_score":1/3,"gdt_ts":1/3,'clashscore':1/12,"inf_all":1/8}
#       Z_local = {"inf_all":1/2,"lddt":1/2}
#       temp_df["Z_rna"] = get_weighted_sum_z(temp_df,Z_rna,"Z_")
#       temp_df["Z_local"] = get_weighted_sum_z(temp_df,Z_local,"Z_")

###############################################################################


# Locations
usalign_location = '/home/anthony/USalign/'
chimerax_location = '/home/anthony/chimerax-1.6-rc2023.03.31//bin/'
phenix_location = '/home/anthony/phenix-1.20.1-4487/build/bin/'

cwd = '/home/anthony/Desktop/CASP/data'


def run_all_on_target(tar_file, natives, models_folder, name,
                      prepare, run, EM=False, resolution=None,
                      threshold=None, center=False, N=None):
    '''
    Prepare and score models for a given target

    Args:
        tar_files (str): tar file of all models that unzips to models_folder
        natives (list of str of dict str:list of str): if not EM, list of native pdbs, if EM dictionary of map:list of natives
        models_folder (str): location of models, after untarred, normally cwd/name
        name (str): name of target (eg R1107 or RT1189)
        prepare (bool): whether to prepare the pdbs, if False assume they are
            already prepared in models_folder
        EM (bool): if to run EM metrics (default False)
        resolution (dict str:float): only if EM, resolution of maps in format map:resolution (default None)
        threshold (dict str:float): only if EM, thresholds to report atomic inclusion in format map:threshold (default None)
        center (bool): only if EM, whether to center map and native (default False)
        N (int): just for testing, run only on N pdbs, if None run on all (default None)
    '''
    if prepare:
        print(f'Preparing models in {models_folder}.')

        # get files ready
        if not path.isdir(f'{cwd}/result_files'):
            system(f'mkdir {cwd}/result_files')
        if path.isdir(models_folder):
            system(f'rm -r {models_folder}')
        system(f'tar -xf {tar_file}')

        # take care of pdb names to match {name}TS???_? pattern
        for x in glob(f'{models_folder}/*TS*'):
            xname = x.rsplit("/", 1)[-1]
            if xname[:len(name)] != name:
                xname = f'{name}TS{xname.rsplit("TS",1)[-1]}'
            if len(xname) == len(name)+5:
                # missing model
                xname = xname+"_1"
            y = f'{x.rsplit("/",1)[0]}/{xname[:len(name)+7]}'
            if x != y:
                system(f'mv {x} {y}')
        for x in glob(f'{models_folder}/*'):
            y = f'{x.rsplit("/",1)[0]}/{x.rsplit("/",1)[-1][:12]}'
            if x != y:
                system(f'mv {x} {y}')

        # copy native to the run folder
        if EM:
            for i, emmap in enumerate(natives):
                for j, native in enumerate(natives[emmap]):
                    system(f'cp {native} {models_folder}/{name}TSex{i}_{j}')
        else:
            for i, native in enumerate(natives):
                system(f'cp {native} {models_folder}/{name}TSexp_{i}')

        # prepare all pdbs
        pdbs = prepare_pdbs(f'{models_folder}/{name}TS???_?')

    else:
        pdbs = glob(f'{models_folder}/{name}TS???_?.pdb')

    if run:
        if N is not None:
            pdbs = pdbs[:N]
        if EM:
            for i, emmap in enumerate(natives):
                for j, native in enumerate(natives[emmap]):
                    native = f'{models_folder}/{name}TSex{i}_{j}.pdb'
                    # only run EM comparison once
                    if j == 0:
                        if center:
                            center_map_and_native(emmap, native)
                        print(f'Scoring models (EM) in {models_folder}, {native}, {emmap}.')
                        score_all(pdbs, f"{cwd}/result_files/{name}_{i}-{j}",
                                  native, usalign_location=usalign_location,
                                  EM=EM, emmap=emmap,
                                  resolution=resolution[emmap],
                                  threshold=threshold[emmap],
                                  chimerax_location=chimerax_location,
                                  phenix_location=phenix_location)
                    else:
                        print(f'Scoring models in {models_folder}, {native}, {emmap}.')
                        score_all(pdbs, f"{cwd}/result_files/{name}_{i}-{j}",
                                  native, usalign_location=usalign_location,
                                  EM=False, phenix_location=phenix_location)
        else:
            for i, native in enumerate(natives):
                native = f'{models_folder}/{name}TSexp_{i}.pdb'
                print(f'Scoring models in {models_folder}, {native}.')
                score_all(pdbs, f"{cwd}/result_files/{name}_{i}", native,
                          usalign_location=usalign_location, EM=EM,
                          phenix_location=phenix_location)


###############################################################################
# R1107 ~9.2s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1107.tar.gz',
                  natives=[f'{cwd}/D_1292119758_model-annotate_P1human.pdb'],
                  models_folder=f'{cwd}/R1107',
                  name="R1107",
                  prepare=True,
                  run=True,
                  N=None)

###############################################################################
# R1108 ~9.2s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1108.tar.gz',
                  natives=[f'{cwd}/D_1292119797_model-annotate_P1chimp_A.pdb',
                           f'{cwd}/D_1292119797_model-annotate_P1chimp_B.pdb'],
                  models_folder=f'{cwd}/R1108',
                  name="R1108",
                  prepare=True,
                  run=True,
                  N=None)

###############################################################################
# R1116 ~11.5s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1116.tar.gz',
                  natives=[f'{cwd}/PV_CL_092022_refine_26.pdb'],
                  models_folder=f'{cwd}/R1116',
                  name="R1116",
                  prepare=True,
                  run=True,
                  N=None)

###############################################################################
# R1117 ~8.2s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1117.tar.gz',
                  natives=[f'{cwd}/R1117_structure_31.pdb'],
                  models_folder=f'{cwd}/R1117',
                  name="R1117",
                  prepare=True,
                  run=True,
                  N=2)

run_all_on_target(tar_file=f'{cwd}/R1117v2.tar.gz',
                  natives=[f'{cwd}/R1117_structure_31.pdb'],
                  models_folder=f'{cwd}/R1117v2',
                  name="R1117v2",
                  prepare=True,
                  run=True,
                  N=2)

###############################################################################
# R1126 ~212s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1126.tar.gz',
                  natives={f'{cwd}/R1126_structure_21.mrc':
                           [f'{cwd}/R1126_structure_21_Traptamer_rsr.pdb']},
                  models_folder=f'{cwd}/R1126',
                  name="R1126",
                  EM=True,
                  prepare=True,
                  run=True,
                  resolution={f'{cwd}/R1126_structure_21.mrc': 5.6},
                  threshold={f'{cwd}/R1126_structure_21.mrc': 0.11},
                  N=2)

###############################################################################
# R1128 ~180s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1128.tar.gz',
                  natives={f'{cwd}/R1128_structure_24.mrc':
                           [f'{cwd}/R1128_structure_24_6wj_refined_046.pdb']},
                  models_folder=f'{cwd}/R1128',
                  name="R1128",
                  EM=True,
                  prepare=True,
                  run=True,
                  resolution={f'{cwd}/R1128_structure_24.mrc': 5.3},
                  threshold={f'{cwd}/R1128_structure_24.mrc': 0.124},
                  N=2)

###############################################################################
# R1136 ~155s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1136.tar.gz',
                  natives={f'{cwd}/R1136v1_structure_22.mrc':
                           [f'{cwd}/R1136v1_Apta_FRET-Bound_rsr15-withAt2OP1.pdb'],
                           f'{cwd}/R1136v2_Apta-FRET-unbound.mrc':
                           [f'{cwd}/R1136v2_structure_22v2_Apta_FRET-unbound_rsr_deposit.pdb']},
                  models_folder=f'{cwd}/R1136',
                  name="R1136",
                  EM=True,
                  prepare=True,
                  run=True,
                  resolution={f'{cwd}/R1136v1_structure_22.mrc': 4.4,
                              f'{cwd}/R1136v2_Apta-FRET-unbound.mrc': 4.5},
                  threshold={f'{cwd}/R1136v1_structure_22.mrc': 0.11,
                              f'{cwd}/R1136v2_Apta-FRET-unbound.mrc': 0.11},
                  N=2)

###############################################################################
# R1138 ~460s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1138.tar.gz',
                  natives={f'{cwd}/R1138v1_structure_23.mrc':
                           [f'{cwd}/R1138v1_structure_23v1_6hbc_young.pdb'],
                           f'{cwd}/R1138v2_6HBC-Mature.mrc':
                           [f'{cwd}/R1138v2_structure_23v2_6hbc_mature.pdb']},
                  models_folder=f'{cwd}/R1138',
                  name="R1138",
                  EM=True,
                  prepare=True,
                  run=True,
                  resolution={f'{cwd}/R1138v1_structure_23.mrc': 5.2,
                              f'{cwd}/R1138v2_6HBC-Mature.mrc': 4.9},
                  threshold={f'{cwd}/R1138v1_structure_23.mrc': 0.1,
                              f'{cwd}/R1138v2_6HBC-Mature.mrc': 0.44},
                  N=2)

###############################################################################
# R1149 ~490s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1149.tar.gz',
                  natives={f'{cwd}/R11149.mrc':
                           [f'{cwd}/R11149_{i}.pdb' for i in range(10)]},
                  models_folder=f'{cwd}/R1149',
                  name="R1149",
                  EM=True,
                  prepare=True,
                  run=True,
                  center=True,
                  resolution={f'{cwd}/R11149.mrc': 4.74},
                  threshold={f'{cwd}/R11149.mrc': 0.345},
                  N=2)

###############################################################################
# R1156 ~270s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1156.tar.gz',
                  natives={f'{cwd}/R1156_conformation1_resolution5.83.mrc':
                           [f'{cwd}/R1156_conformation1_model{i+1}.pdb' for i in range(10)],
                           f'{cwd}/R1156_conformation2_resolution6.59.mrc':
                           [f'{cwd}/R1156_conformation2_model{i+1}.pdb' for i in range(10)],
                           f'{cwd}/R1156_conformation3_resolution7.48.mrc':
                           [f'{cwd}/R1156_conformation3_model{i+1}.pdb' for i in range(10)],
                           f'{cwd}/R1156_conformation4_resolution7.61.mrc':
                           [f'{cwd}/R1156_conformation4_model{i+1}.pdb' for i in range(10)]},
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
                  N=2)

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
                  N=2)

run_all_on_target(tar_file=f'{cwd}/RT1189.tar.gz',
                  natives=[f'{cwd}/RT1189_A-6B.pdb'],
                  models_folder=f'{cwd}/RT1189',
                  name="RT1189",
                  prepare=True,
                  run=True,
                  N=2)

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
                  N=2)

run_all_on_target(tar_file=f'{cwd}/RT1190.tar.gz',
                  natives=[f'{cwd}/RT1190_A-4B.pdb'],
                  models_folder=f'{cwd}/RT1190',
                  name="RT1190",
                  prepare=True,
                  run=True,
                  N=2)
