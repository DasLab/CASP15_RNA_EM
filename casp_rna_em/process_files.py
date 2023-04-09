from glob import glob
import subprocess as sp
import mrcfile
from os import system, path, remove
from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np
from casp_rna_em.run_metric_programs import run_usalign


def prepare_pdbs(pdbs, ignore_no_extension=False, fill_occupancy=True,
                 fill_bfactor=True, remove_protein=True):
    '''
    Prepare all pdbs to ready for scoring

    Args:
        pdbs (str): string of pdb locations, globbed
        ignore_no_extension (bool): whether to assume a file without extension is a pdb and use it or not (default False)
        fill_occupancy (bool): whether to fill all occupancy to 1, some scores do weigh by occupancy (default True)
        fill_bfactor (bool): whether to fill bfactor to 0 (default True)
        remove_protein (bool): whether to remove all protein residues, in fact any none A,C,G,U residue (default True)

    Returns:
        list of pdbs which have now been prepared
    '''
    pdbs = glob(pdbs)

    # only look at file with pdb extension
    # or if have no extension, assume pdb, .pdb
    if ignore_no_extension:
        pdbs = [x for x in pdbs if pdbs.rsplit(".", 1)[-1] == "pdb"]
    else:
        pdbs_a = [x for x in pdbs if x.rsplit(".", 1)[-1] == "pdb"]
        pdbs_to_rename = [x for x in pdbs if "." not in x]
        pdbs_b = [f'{x}.pdb' for x in pdbs_to_rename]
        for pdb in pdbs_to_rename:
            system(f'mv {pdb} {pdb}.pdb')
        pdbs = pdbs_a + pdbs_b
    # copy to keep a unpreparred version
    for pdb in pdbs:
        system(f'cp {pdb} {pdb[:-4]}_unprepared.pdb')

    # clean each pdb
    for pdb in pdbs:
        clean_pdb(pdb, fill_occupancy=fill_occupancy,
                  fill_bfactor=fill_bfactor, remove_protein=remove_protein)
    return pdbs


def clean_pdb(pdb, fill_occupancy=True, fill_bfactor=True,
              remove_protein=True):
    '''
    Helper function of prepare+pdbs, refer above for arguments
    '''
    if fill_occupancy or fill_bfactor or remove_protein:
        ppdb = PandasPdb().read_pdb(pdb)
        # make all atoms have full occupanncy
        if fill_occupancy:
            ppdb.df['ATOM'].occupancy = 1
        # remove bfactor information
        if fill_bfactor:
            ppdb.df['ATOM'].b_factor = 0
        if remove_protein:
            RNA_nucs = ["A", "C", "G", "U"]
            at_df = ppdb.df['ATOM']
            all_residues = at_df.residue_name.unique()
            removed_residues = [x for x in all_residues if x not in RNA_nucs]
            if removed_residues != []:
                print(f'Removed {removed_residues} from {pdb}')
            ppdb.df['ATOM'] = at_df[at_df.residue_name.isin(
                RNA_nucs)]
        ppdb.to_pdb(path=pdb)

    # use rna tools to prepare pdbs, action completed:
    # remove any modeled H so they are not scored
    # renames chains
    command = ['rna_pdb_tools.py', '--get-rnapuzzle-ready',
               '--inplace', "--renumber-residues", '--delete-anisou', pdb]
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print(f'ERROR: rna tools failed: on {pdb}\n{out.decode()}\n{err.decode()}')

    '''
    old code that seems covered by RNA tools
        drop any duplicated atoms
        ppdb.df['ATOM'].drop_duplicates(subset=ppdb.df['ATOM'].columns[:-1],inplace=True)

        147, 227 has first residue atom an "O" which is not recognized atom
        for f in ${puzzle}/${puzzle}TS???_?.pdb; do sed -i '/^ATOM      1  O     G/d' $f; done
        some group 035, 054, 125 have cysteines N4 as a O in the element column
        this crashes programs, so correctly change to a N
    '''


def center_map_and_native(emmap, pdb):
    '''
    Make origin of emmap 0,0,0 and move pdb along as well
    Note if center is already 0,0,0 no translation occur

    Args:
        emmap (str): location of cryoEM map to be centered
        pdb (str): location of pdb to be moved same as emmap
    '''
    system(f'cp {emmap} {emmap.rsplit(".")[0]}_uncentered.mrc')
    system(f'cp {pdb} {pdb.rsplit(".")[0]}_uncentered.mrc')
    with mrcfile.open(emmap, 'r+') as mrc:
        x_, y_, z_ = mrc.header.origin.x.item(
        ), mrc.header.origin.y.item(), mrc.header.origin.z.item()
        mrc.header.origin.x = 0
        mrc.header.origin.y = 0
        mrc.header.origin.z = 0
    ppdb = PandasPdb().read_pdb(pdb)
    ppdb.df['ATOM'].x_coord = ppdb.df['ATOM'].x_coord - x_
    ppdb.df['ATOM'].y_coord = ppdb.df['ATOM'].y_coord - y_
    ppdb.df['ATOM'].z_coord = ppdb.df['ATOM'].z_coord - z_
    ppdb.to_pdb(path=pdb)


def write_and_sbatch_scoring(pdbs, N, prefix, name, base_sbatch,
                             out_file_prefix, native, usalign_location='', EM=False,
                             phenix_location='', chimerax_location='',
                             emmap=None, resolution=None, threshold=None,
                             sbatch=True):
    '''
    Split pdbs up into N groups and sbatch off the scoring of those groups

    Args:
        pdbs (list of str): a list of pdbs
        N (int): number of groups to split into
        prefix (str): prefix of text files to write pdb list and sbatch file to
        name (str): target name, will save outfiles and sbatchfiles accordingly
        base_sbatch (str): example/header sbatch file
        out_file_prefix (str): prefix for all result files
        native (str): location of native pdb
        usalign_location (str): if empty then USalign should be in path, location of usalign, this folder should contain USalign (default='')
        EM (bool): if should score EM matrics (default False)
        phenix_location (str): if empty then phenix executable should be in path, location of phenix, this folder should contain for example phenix.rna_validate (default='') 
        chimerax_location (str): if empty then ChimeraX should be in path, location of ChimeraX, this folder should contain ChimeraX executable (default '')
        emmap (str): cryo-EM map location (default None)
        resolution (float): only if EM, resolution of cryo-EM map (default None)
        threshold (float): only if EM, threshold to use for atomic inclusion (default None)
        sbatch (bool): if true sbatch the jobs, otherwise just run them all (default True)
    '''
    seperate_pdbs_to_run_groups(pdbs, N, prefix)
    for i in range(N):
        system(f'cp {base_sbatch} {prefix}_{i}.sbatch')
        system(f'sed -i s/XXX/{name}_{i}/g {prefix}_{i}.sbatch')
        with open(f'{prefix}_{i}.sbatch', 'a') as f:
            #f.write('\nconda activate casp_rna_em \n')
            if EM:
                f.write(f'python {path.dirname(__file__)}/casp_rna_score_all_from_file.py -f {prefix}_{i}.txt --out_file_prefix {out_file_prefix}_{i} --native {native} --usalign_location {usalign_location} --chimerax_location {chimerax_location} --EM --emmap {emmap} --resolution {resolution} --threshold {threshold} --phenix_location {phenix_location} \n')
            else:
                f.write(f'python {path.dirname(__file__)}/casp_rna_score_all_from_file.py -f {prefix}_{i}.txt --out_file_prefix {out_file_prefix}_{i} --native {native} --usalign_location {usalign_location} --phenix_location {phenix_location} \n')
        if sbatch:
            system(f'sbatch {prefix}_{i}.sbatch')
        else:
            system(f'chmod +x {prefix}_{i}.sbatch')
            system(f'./{prefix}_{i}.sbatch &')


def seperate_pdbs_to_run_groups(pdbs, N, prefix):
    '''
    split pdbs into N groups and save lists in files prefix_i.txt

    Args:
        pdbs (list of str): a list of pdbs
        N (int): number of groups to split into
        prefix (str): prefix of text files to write pdb list to
    '''
    num = (len(pdbs) + 1) // N
    for i in range(N):
        with open(f'{prefix}_{i}.txt', 'w') as f:
            for pdb in pdbs[i * num:min((i + 1) * num, len(pdbs))]:
                f.write(f'{pdb}\n')


def concat_all_result_files(files, out_file):
    '''
    saves a csv of all the files combines

    Args:
        files (list of str): list of csv to combine
        out_file (str): name of file to save combined data to
    '''
    dfs = []
    for f in files:
        dfs.append(pd.read_csv(f))
    dfs = pd.concat(dfs)
    dfs.to_csv(out_file, index=False)


# A default metric_dict for CASP
METRICS = {"clashscore": "min", "global_rmsd": "min", "cc_mask": "max",
           "cc_peaks": "max", "tempy_mi": "max", "tempy_smoc": "max",
           "ai": "max", "gdt_ts": "max", "tm_score": "max", "lddt": "max",
           "inf_all": "max"}


def reduce_df(df, score_to_choice_best=None,
              static_columns=["target", "gr_code", "model"],
              metric_dict=METRICS):
    ''' 
    Reduce df to "best" score per X 
        (vs multiple scores per X if there are multiple Y)
    Eg for conformations:
    If score_to_choice_best is None then, for each score, the best value over all conformations is selected
    If it is specified then first the best conformation according to that score is selected,
    then all scores are taken from the conformation.

    Args:
        df (DataFrame): DataFrame with all columns in stat_columns and metric_dict, will combine all rows for any other columns
        score_to_choice_best (str): if None then for each score in metric_dict take best value according to agg function, 
            if it is specified then take all values from the same row which is the best score according to this column (default None)
        static_columns (list of str): The columns to group by
        metric_dict (dict str(metric):str(agg function)): the dataframe column names and the function to aggregate scores (eg 'clashscore':'min')

    Returns:
        DataFrame reduced to one value per static_columns by the aggregation in metric_dict
    '''
    columns_combining = [x for x in df.columns if x not in list(
        metric_dict.keys()) + static_columns]
    print("combing the following columns:", columns_combining)

    if score_to_choice_best is not None:
        if metric_dict[score_to_choice_best] == "max":
            best_index = df.groupby(static_columns)[
                score_to_choice_best].idxmax().to_numpy(copy=True)
        elif metric_dict[score_to_choice_best] == "min":
            best_index = df.groupby(static_columns)[
                score_to_choice_best].idxmin().to_numpy(copy=True)
        return df.iloc[best_index].drop(columns_combining, axis=1)
    else:
        return df.groupby(static_columns).agg(metric_dict).reset_index()


def calc_z(arr, raw):
    '''
    Calculate Z of second array with mean and std of the first

    Args:
        arr (array): is values to calculate mean and std of
        raw (array): is values to calculate z with calculated mean and std

    Returns:
        array of zscore
    '''
    u = np.nanmean(arr)
    s = np.nanstd(arr)
    z = (raw - u) / s
    return z


def get_zscore(arr, negative=False, threshold=-2, mask=None):
    '''
    Calculate Zscore of an array

    from CASP14
    1. Calculate zscores from the raw scores for all "first" models (corresponding values from the main result table);
    2. Remove outliers - models with zscores below the tolerance threshold (set to -2.0);
    3. Recalculate zscores on the reduced dataset;
    4. Assign z-scores below the penalty threshold (either -2.0 or 0.0) to the value of this threshold.

    Args:
        arr (array): the array of scores to get z-score of, will ignore nans
        negative (bool): true if the minimum score is the best score (default False)
        threshold (float): after initial mean and std calculations, any Zscore below -2 is ignore for the final mean and
            threshold calculation, the final Zscores are also capped at threshold (default -2)
        mask (array): will ignore any scores in arr where mask is False (default None)

    Returns:
        Array of the Z-scores
    '''

    # copy raw scores to calculate final z later
    raw = arr.astype(float)
    old_raw = np.copy(raw)

    # mask out values for mean and std calcualtion
    if mask is not None:
        raw[~mask] = np.nan

    # set multiple if large or small is good
    if negative:
        mult = -1.
    else:
        mult = 1.

    # calculate z score
    zscore = calc_z(raw * mult, raw * mult)
    # ignore negative outliers and recalculate
    # mean and std, calculate z over all raw values
    raw[zscore < threshold] = np.nan
    zscore = calc_z(raw * mult, old_raw * mult)
    # cap z score at threshold
    zscore[zscore < threshold] = threshold
    return zscore


def get_weighted_sum_z(data, weights, prefix=''):
    '''
    Obtain the weight sum of scores
    Weights should sum to 1
    Note does to explicitly require weights to sum to 1

    Args:
        data (DataFrame): contains columns for each prefix+score in weights
        weights (dict str(score):float(weight)): the weights for each score
        prefix (str): if column names have a prefix (eg 'Z_') infront of the 
            scores specified in weights, specify this here (default '')

    Returns:
        Series of the summed scores
    '''
    final_z = None
    for score, weight in weights.items():
        if final_z is None:
            final_z = weight * data[prefix + score].to_numpy()
        else:
            final_z += weight * data[prefix + score].to_numpy()
    return final_z


def get_group_score(df, agg="sum>0", score="Z_rna"):
    '''
    Aggregate the scores of each group in gr_code

    Args:
        df (DataFrame): contain gr_code and score (argument specified) column
        agg (str): how to aggregate scores options: 'sum>0','mean','sum' (default 'sum>0')
        score (str): what column to aggregate

    Returns:
        DataFrame of groups and their summed score
    '''
    vals = df.groupby("gr_code")[score]
    if agg == "sum>0":
        return vals.apply(lambda col: col[col > 0].sum()).reset_index()
    elif agg == "mean":
        return vals.mean().reset_index()
    elif agg == "sum":
        return vals.sum().reset_index()
