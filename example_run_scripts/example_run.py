from casp_rna_em.process_files import prepare_pdbs, concat_all_result_files, center_map_and_native, reduce_df, get_zscore, get_weighted_sum_z
from casp_rna_em.run_metric_programs import score_all
from os import system, path
import pandas as pd
from glob import glob


###############################################################################
# TODO
###############################################################################
# code
#   push analysis edits

# Testing
#   test sbatch on sherlock
#   production run

# Documentation
#   pipeline

# Later:
#   options to run: ddsr, qscore, fsc
#   qscore: fork mapq with parralelizaiton and sigfig fixes
#   add in segmenetation code
#   Full?
#       lddt and inf not calculated, not sure which rmsd as well?
#       Z_rna = {"lddt":1/8,"tm_score":1/3,"gdt_ts":1/3,'clashscore':1/12,"inf_all":1/8}
#       Z_local = {"inf_all":1/2,"lddt":1/2}
#       temp_df["Z_rna"] = get_weighted_sum_z(temp_df,Z_rna,"Z_")
#       temp_df["Z_local"] = get_weighted_sum_z(temp_df,Z_local,"Z_")

###############################################################################


usalign_location = '/home/anthony/USalign/'
chimerax_location = '/home/anthony/chimerax-1.6-rc2023.03.31//bin/'
phenix_location = '/home/anthony/phenix-1.20.1-4487/build/bin/'

cwd = '/home/anthony/Desktop/CASP/data'


def run_all_on_target(tar_file, natives, models_folder, name,
                      prepare, run, EM=False, resolution=None, threshold=None, center=False, N=None):
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
                        score_all(pdbs, f"result_files/{name}_{i}-{j}", native, usalign_location=usalign_location, EM=EM, emmap=emmap, resolution=resolution[emmap], threshold=threshold[emmap], chimerax_location=chimerax_location, phenix_location=phenix_location)
                    else:
                        print(f'Scoring models in {models_folder}, {native}, {emmap}.')
                        score_all(pdbs, f"result_files/{name}_{i}-{j}", native, usalign_location=usalign_location, EM=False, phenix_location=phenix_location)
        else:
            for i, native in enumerate(natives):
                native = f'{models_folder}/{name}TSexp_{i}.pdb'
                print(f'Scoring models in {models_folder}, {native}.')
                score_all(pdbs, f"result_files/{name}_{i}", native, usalign_location=usalign_location, EM=EM, phenix_location=phenix_location)


###############################################################################
# R1107 ~9.2s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1107.tar.gz',
                  natives=[f'{cwd}/D_1292119758_model-annotate_P1human.pdb'],
                  models_folder=f'{cwd}/R1107',
                  name="R1107",
                  prepare=False,
                  run=False,
                  N=None)

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
                  N=None)

###############################################################################
# R1116 ~11.5s/pdb
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1116.tar.gz',
                  natives=[f'{cwd}/PV_CL_092022_refine_26.pdb'],
                  models_folder=f'{cwd}/R1116',
                  name="R1116",
                  prepare=False,
                  run=False,
                  N=None)

###############################################################################
# R1117
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1117.tar.gz',
                  natives=[f'{cwd}/R1117_structure_31.pdb'],
                  models_folder=f'{cwd}/R1117',
                  name="R1117",
                  prepare=False,
                  run=False,
                  N=2)

run_all_on_target(tar_file=f'{cwd}/R1117v2.tar.gz',
                  natives=[f'{cwd}/R1117_structure_31.pdb'],
                  models_folder=f'{cwd}/R1117v2',
                  name="R1117v2",
                  prepare=False,
                  run=False,
                  N=2)

###############################################################################
# R1126
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1126.tar.gz',
                  natives={f'{cwd}/R1126_structure_21.mrc': [f'{cwd}/R1126_structure_21_Traptamer_rsr.pdb']},
                  models_folder=f'{cwd}/R1126',
                  name="R1126",
                  EM=True,
                  prepare=False,
                  run=False,
                  resolution={f'{cwd}/R1126_structure_21.mrc': 5.6},
                  threshold={f'{cwd}/R1126_structure_21.mrc': 0.11},
                  N=3)

###############################################################################
# R1128
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1128.tar.gz',
                  natives={f'{cwd}/R1128_structure_24.mrc': [f'{cwd}/R1128_structure_24_6wj_refined_046.pdb']},
                  models_folder=f'{cwd}/R1128',
                  name="R1128",
                  EM=True,
                  prepare=False,
                  run=False,
                  resolution={f'{cwd}/R1128_structure_24.mrc': 5.3},
                  threshold={f'{cwd}/R1128_structure_24.mrc': 0.124},
                  N=2)

###############################################################################
# R1136
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1136.tar.gz',
                  natives={f'{cwd}/R1136v1_structure_22.mrc': [f'{cwd}/R1136v1_Apta_FRET-Bound_rsr15-withAt2OP1.pdb'],
                           f'{cwd}/R1136v2_Apta-FRET-unbound.mrc': [f'{cwd}/R1136v2_structure_22v2_Apta_FRET-unbound_rsr_deposit.pdb']},
                  models_folder=f'{cwd}/R1136',
                  name="R1136",
                  EM=True,
                  prepare=False,
                  run=False,
                  resolution={f'{cwd}/R1136v1_structure_22.mrc': 4.4,
                              f'{cwd}/R1136v2_Apta-FRET-unbound.mrc': 4.5},
                  threshold={f'{cwd}/R1136v1_structure_22.mrc': 0.11,
                              f'{cwd}/R1136v2_Apta-FRET-unbound.mrc': 0.11},
                  N=2)

###############################################################################
# R1138
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1138.tar.gz',
                  natives={f'{cwd}/R1138v1_structure_23.mrc': [f'{cwd}/R1138v1_structure_23v1_6hbc_young.pdb'],
                           f'{cwd}/R1138v2_6HBC-Mature.mrc': [f'{cwd}/R1138v2_structure_23v2_6hbc_mature.pdb']},
                  models_folder=f'{cwd}/R1138',
                  name="R1138",
                  EM=True,
                  prepare=False,
                  run=False,
                  resolution={f'{cwd}/R1138v1_structure_23.mrc': 5.2,
                              f'{cwd}/R1138v2_6HBC-Mature.mrc': 4.9},
                  threshold={f'{cwd}/R1138v1_structure_23.mrc': 0.1,
                              f'{cwd}/R1138v2_6HBC-Mature.mrc': 0.44},
                  N=2)

###############################################################################
# R1149
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1149.tar.gz',
                  natives={f'{cwd}/R11149.mrc': [f'{cwd}/R11149_{i}.pdb' for i in range(10)]},
                  models_folder=f'{cwd}/R1149',
                  name="R1149",
                  EM=True,
                  prepare=False,
                  run=False,
                  center=True,
                  resolution={f'{cwd}/R11149.mrc': 4.74},
                  threshold={f'{cwd}/R11149.mrc': 0.345},
                  N=2)

###############################################################################
# R1156
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1156.tar.gz',
                  natives={f'{cwd}/R1156_conformation1_resolution5.83.mrc': [f'{cwd}/R1156_conformation1_model{i+1}.pdb' for i in range(10)],
                           f'{cwd}/R1156_conformation2_resolution6.59.mrc': [f'{cwd}/R1156_conformation2_model{i+1}.pdb' for i in range(10)],
                           f'{cwd}/R1156_conformation3_resolution7.48.mrc': [f'{cwd}/R1156_conformation3_model{i+1}.pdb' for i in range(10)],
                           f'{cwd}/R1156_conformation4_resolution7.61.mrc': [f'{cwd}/R1156_conformation4_model{i+1}.pdb' for i in range(10)]},
                  models_folder=f'{cwd}/R1156',
                  name="R1156",
                  EM=True,
                  prepare=False,
                  run=False,
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
# R1189
# EM not run due to poor fits
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1189.tar.gz',
                  natives=[f'{cwd}/RT1189_A-6B.pdb'],
                  models_folder=f'{cwd}/R1189',
                  name="R1189",
                  prepare=False,
                  run=False,
                  N=2)

run_all_on_target(tar_file=f'{cwd}/RT1189.tar.gz',
                  natives=[f'{cwd}/RT1189_A-6B.pdb'],
                  models_folder=f'{cwd}/RT1189',
                  name="RT1189",
                  prepare=False,
                  run=False,
                  N=2)

###############################################################################
# R1190
# EM not run due to poor fits
###############################################################################

run_all_on_target(tar_file=f'{cwd}/R1190.tar.gz',
                  natives=[f'{cwd}/RT1190_A-4B.pdb'],
                  models_folder=f'{cwd}/R1190',
                  name="R1190",
                  prepare=False,
                  run=False,
                  N=2)

run_all_on_target(tar_file=f'{cwd}/RT1190.tar.gz',
                  natives=[f'{cwd}/RT1190_A-4B.pdb'],
                  models_folder=f'{cwd}/RT1190',
                  name="RT1190",
                  prepare=False,
                  run=False,
                  N=2)

###############################################################################
# Combine all
###############################################################################

metrics = {'angle_outlier': 'min', 'avg_suitness': 'max', 'suite_outlier': 'min',
           'bond_outlier': 'min', 'clashscore': 'min', 'pucker_outlier': 'min',
           'rmsd_lga': 'min', 'gdt_ts': 'max', 'rmsd_usalign': 'min', 'tm': 'max',
           'cc_volume': 'max', 'cc_mask': 'max', 'cc_peaks': 'max',
           'ai': 'max', 'ai_bb': 'max', 'ai_base': 'max', 'do': 'max',
           'do_bb': 'max', 'do_base': 'max',
           'tempy_ccc': 'max', 'tempy_mi': 'max',
           'tempy_lsf': 'max', 'tempy_env': 'max',
           'tempy_sccc': 'max', 'tempy_smoc': 'max',
           'tempy_smoc_1': 'max'}

concat_all_result_files(glob('result_files/*scores.csv'), "all_scores.csv")
df = pd.read_csv("all_scores.csv")
columns = ['target', 'gr_code', 'model', 'emmap', ]
df = reduce_df(df, score_to_choice_best=None, static_columns=columns, metric_dict=metrics)
df.to_csv('all_scores_best_per_conformation_per_model.pdb', index=False)
columns = ['target', 'gr_code', 'model']
df = reduce_df(df, score_to_choice_best=None, static_columns=columns, metric_dict=metrics)
df.to_csv('all_scores_best_per_model.pdb', index=False)
clean_columns = ['target', 'gr_code', 'model', 'cc_mask', 'cc_peaks', 'tempy_mi', 'tempy_smoc', 'ai', "tm", "gdt_ts"]
df[clean_columns].to_csv('all_scores_best_per_model_clean.pdb', index=False)

columns = ['target', 'gr_code']
df = reduce_df(df, score_to_choice_best=None, static_columns=columns, metric_dict=metrics)
df.to_csv('all_scores_best_per_group.pdb', index=False)

# for each metric and puzzle, calcualte Z-scores
for metric in metrics.keys():
    df["Z_" + metric] = get_zscore(df[metric].to_numpy(copy=True), negative=(metrics[metric] == "min"), threshold=-2)

# get weighted sum Zs of interest
Z_em = {'cc_mask': 1 / 5, 'cc_peaks': 1 / 5, 'tempy_mi': 1 / 5, 'tempy_smoc': 1 / 5, 'ai': 1 / 5}
Z_topo = {"tm": 1 / 2, "gdt_ts": 1 / 2}
df["Z_em"] = get_weighted_sum_z(df, Z_em, "Z_")
df["Z_topo"] = get_weighted_sum_z(df, Z_topo, "Z_")
df.to_csv('all_scores_best_per_group_withZ.pdb', index=False)
clean_columns = ['target', 'gr_code', 'cc_mask', 'cc_peaks', 'tempy_mi', 'tempy_smoc', 'ai', "tm", "gdt_ts",
                 'Z_cc_mask', 'Z_cc_peaks', 'Z_tempy_mi', 'Z_tempy_smoc', 'Z_ai', "Z_tm", "Z_gdt_ts",
                 'Z_em', 'Z_topo']
df[clean_columns].to_csv('all_scores_best_per_group_withZ_clean.pdb', index=False)
