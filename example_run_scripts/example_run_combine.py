import pandas as pd
from glob import glob
from casp_rna_em.process_files import reduce_df, concat_all_result_files
from casp_rna_em.process_files import get_weighted_sum_z, get_zscore

###############################################################################
# Combine all
###############################################################################

metrics = {'angle_outlier': 'min', 'avg_suitness': 'max', 'clashscore': 'min',
           'bond_outlier': 'min', 'pucker_outlier': 'min',
           'suite_outlier': 'min', 'rmsd_lga': 'min', 'gdt_ts': 'max',
           'rmsd_usalign': 'min', 'tm': 'max',
           'cc_volume': 'max', 'cc_mask': 'max', 'cc_peaks': 'max',
           'ai': 'max', 'ai_bb': 'max', 'ai_base': 'max', 'do': 'max',
           'do_bb': 'max', 'do_base': 'max',
           'tempy_ccc': 'max', 'tempy_mi': 'max',
           'tempy_lsf': 'max', 'tempy_env': 'max',
           'tempy_sccc': 'max', 'tempy_smoc': 'max',
           'tempy_smoc_1': 'max'}
           
# combine all
concat_all_result_files(glob(f'{cwd}/result_files/*scores.csv'),
                        "all_scores.csv")
concat_all_result_files(glob(f'{cwd}/result_files/*per_residue.csv'),
                        "all_per_residue.csv")
concat_all_result_files(glob(f'{cwd}/result_files/*per_threshold.csv'),
                        "all_per_threshold.csv")

# reduce to one score per model per conformation
df = pd.read_csv("all_scores.csv")
columns = ['target', 'gr_code', 'model', 'emmap']
df = reduce_df(df, score_to_choice_best=None,
               static_columns=columns, metric_dict=metrics)
df.to_csv('all_scores_best_per_conformation_per_model.csv', index=False)

# reduce to one score per model
# save clean version with only scores of interest
columns = ['target', 'gr_code', 'model']
df = reduce_df(df, score_to_choice_best=None,
               static_columns=columns, metric_dict=metrics)
df.to_csv('all_scores_best_per_model.csv', index=False)
clean_columns = ['target', 'gr_code', 'model', 'cc_mask', 'cc_peaks',
                 'tempy_mi', 'tempy_smoc', 'ai', "tm", "gdt_ts", "clashscore"]
df[clean_columns].to_csv('all_scores_best_per_model_clean.csv', index=False)

# reduce to one score per group
columns = ['target', 'gr_code']
df = reduce_df(df, score_to_choice_best=None,
               static_columns=columns, metric_dict=metrics)
df.to_csv('all_scores_best_per_group.csv', index=False)

# for each metric and puzzle, calcualte Z-scores
temp_dfs = []
for target in temp_df.target.unique():
    target_df = df[df.target==target].copy()
    for metric in metrics.keys():
        target_df["Z_"+metric] = get_zscore(target_df[metric].to_numpy(copy=True),
                                            negative=(METRICS[metric]=="min"),
                                            threshold=-2)
    temp_dfs.append(target_df)
df = pd.concat(temp_dfs)

# get weighted sum Zs of interest
# save a clean version with only scores of interest
Z_em = {'cc_mask': 1 / 5, 'cc_peaks': 1 / 5,
        'tempy_mi': 1 / 5, 'tempy_smoc': 1 / 5, 'ai': 1 / 5}
Z_topo = {"tm": 1 / 2, "gdt_ts": 1 / 2}
df["Z_em"] = get_weighted_sum_z(df, Z_em, "Z_")
df["Z_topo"] = get_weighted_sum_z(df, Z_topo, "Z_")
df.to_csv('all_scores_best_per_group_withZ.csv', index=False)
clean_columns = ['target', 'gr_code', 'cc_mask', 'cc_peaks', 'tempy_mi',
                 'tempy_smoc', 'ai', "tm", "gdt_ts", "clashscore",
                 'Z_cc_mask', 'Z_cc_peaks', 'Z_tempy_mi', 'Z_tempy_smoc',
                 'Z_ai', "Z_tm", "Z_gdt_ts", "Z_clashscore",
                 'Z_em', 'Z_topo']
df[clean_columns].to_csv(
    'all_scores_best_per_group_withZ_clean.csv', index=False)

