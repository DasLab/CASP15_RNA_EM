import pandas as pd
import seaborn as sns
import numpy as np

METRICS = {"clashscore":"min","global_rmsd":"min","cc_mask":"max","cc_peaks":"max",
          "tempy_mi":"max","tempy_smoc":"max","ai":"max",
          "gdt_ts":"max","tm_score":"max","lddt":"max","inf_all":"max"}

def reduce_df(df,score_to_choice_best=None,static_columns=["target","gr_code","model"]):
    ''' 
    Reduce df to "best" score per X 
        (vs multiple scores per X if there are multiple Y)
    Eg for conformations:
    If score_to_choice_best is None then, for each score, the best value over all coformations is selected
    If it is specified then first the best conformation according to that score is selected,
    then all scores are taken from the conformation.
    '''
    columns_combining = [x for x in df.columns if x not in list(METRICS.keys())+static_columns]
    print("combing the following columns:", columns_combining)
    
    if score_to_choice_best is not None:
        if METRICS[score_to_choice_best] == "max":
            best_index = df.groupby(static_columns)[score_to_choice_best].idxmax().to_numpy(copy=True)
        elif METRICS[score_to_choice_best] == "min":
            best_index = df.groupby(static_columns)[score_to_choice_best].idxmin().to_numpy(copy=True)
        return df.iloc[best_index].drop(columns_combining,axis=1)
    else:
        return df.groupby(static_columns).agg(METRICS).reset_index()
    
def calc_z(arr,raw):
    '''
    Calculate Z of second array with mean and std of the first
    # arr is values to calculate mean and std of
    # raw is values to calculate z with calculated mean and std
    '''
    u = np.nanmean(arr)
    s = np.nanstd(arr)
    z = (raw-u)/s
    return z

def get_zscore(arr,negative=False,threshold=-2,mask=None):
    # from CASP14
    # 1. Calculate zscores from the raw scores for all "first" models (corresponding values from the main result table);
    # 2. Remove outliers - models with zscores below the tolerance threshold (set to -2.0);
    # 3. Recalculate zscores on the reduced dataset;
    # 4. Assign z-scores below the penalty threshold (either -2.0 or 0.0) to the value of this threshold.
    
    # copy raw scores to calculate final z later
    raw = arr.astype(float)
    old_raw = np.copy(raw)

    # mask out values for mean and std calcualtion
    if mask is not None:
        raw[~mask] = np.nan

    # set multiple if large or small is good
    if negative: mult = -1.
    else: mult = 1.

    # calculate z score
    zscore = calc_z(raw*mult,raw*mult)
    # ignore negative outliers and recalculate
    # mean and std, calculate z over all raw values
    raw[zscore<threshold] = np.nan
    zscore = calc_z(raw*mult,old_raw*mult)
    # cap z score at threshold
    zscore[zscore<threshold] = threshold
    return zscore


def get_weighted_sum_z(data,weights,prefix):
    final_z = None
    for score,weight in weights.items():
        if final_z is None:
            final_z = weight*data[prefix+score].to_numpy()
        else:
            final_z += weight*data[prefix+score].to_numpy()
    return final_z
    
def get_group_score(df,agg="sum>0",score="Z_rna"):
    if agg=="sum>0":
        return df.groupby("gr_code")[score].apply(lambda col: col[col > 0].sum()).reset_index()
    elif agg=="mean":
        return df.groupby("gr_code")[score].mean().reset_index()
    
