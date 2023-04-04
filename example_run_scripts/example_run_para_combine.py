import pandas as pd
from casp_rna_em.process_files import reduce_df, get_zscore, get_weighted_sum_z

concat_all_result_files(glob('result_files/*scores.csv'), "all_scores.csv")
concat_all_result_files(glob('result_files/*per_residue.csv'), "all_per_residue.csv")
concat_all_result_files(glob('result_files/*per_threshold.csv'), "all_per_threshold.csv")

df = pd.read_csv("all_scores.csv")
columns = ['target', 'gr_code', 'model', 'emmap', ]
df = reduce_df(df, score_to_choice_best=None, static_columns=columns, metric_dict=metrics)
df.to_csv('all_scores_best_per_conformation_per_model.csv', index=False)
columns = ['target', 'gr_code', 'model']
df = reduce_df(df, score_to_choice_best=None, static_columns=columns, metric_dict=metrics)
df.to_csv('all_scores_best_per_model.csv', index=False)
clean_columns = ['target', 'gr_code', 'model', 'cc_mask', 'cc_peaks', 'tempy_mi', 'tempy_smoc', 'ai', "tm", "gdt_ts"]
df[clean_columns].to_csv('all_scores_best_per_model_clean.csv', index=False)

columns = ['target', 'gr_code']
df = reduce_df(df, score_to_choice_best=None, static_columns=columns, metric_dict=metrics)
df.to_csv('all_scores_best_per_group.csv', index=False)

# for each metric and puzzle, calcualte Z-scores
for metric in metrics.keys():
    df["Z_" + metric] = get_zscore(df[metric].to_numpy(copy=True), negative=(metrics[metric] == "min"), threshold=-2)

# get weighted sum Zs of interest
Z_em = {'cc_mask': 1 / 5, 'cc_peaks': 1 / 5, 'tempy_mi': 1 / 5, 'tempy_smoc': 1 / 5, 'ai': 1 / 5}
Z_topo = {"tm": 1 / 2, "gdt_ts": 1 / 2}
df["Z_em"] = get_weighted_sum_z(df, Z_em, "Z_")
df["Z_topo"] = get_weighted_sum_z(df, Z_topo, "Z_")
df.to_csv('all_scores_best_per_group_withZ.csv', index=False)
clean_columns = ['target', 'gr_code', 'cc_mask', 'cc_peaks', 'tempy_mi', 'tempy_smoc', 'ai', "tm", "gdt_ts",
                 'Z_cc_mask', 'Z_cc_peaks', 'Z_tempy_mi', 'Z_tempy_smoc', 'Z_ai', "Z_tm", "Z_gdt_ts",
                 'Z_em', 'Z_topo']
df[clean_columns].to_csv('all_scores_best_per_group_withZ_clean.csv', index=False)

