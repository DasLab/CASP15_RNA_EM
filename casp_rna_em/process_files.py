from glob import glob
import subprocess as sp
import mrcfile
from os import system, path, remove
from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np
from casp_rna_em.run_metric_programs import run_usalign


def prepare_pdbs(pdbs, ignore_no_extension=False, fill_occupancy=True, fill_bfactor=True, remove_protein=True):
    pdbs = glob(pdbs)

    # only look at file with pdb extension or no extension in which case add pdb
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
        clean_pdb(pdb, fill_occupancy=fill_occupancy, fill_bfactor=fill_bfactor, remove_protein=remove_protein)
    return pdbs


def clean_pdb(pdb, fill_occupancy=True, fill_bfactor=True, remove_protein=True):

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
            all_residues = ppdb.df['ATOM'].residue_name.unique()
            removed_residues = [x for x in all_residues if x not in RNA_nucs]
            if removed_residues != []:
                print(f'Removed {removed_residues} from {pdb}')
            ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM'].residue_name.isin(RNA_nucs)]
        ppdb.to_pdb(path=pdb)

    # use rna tools to prepare pdbs, action completed:
    # remove any modeled H so they are not scored
    # renames chains
    command = ['rna_pdb_tools.py', '--get-rnapuzzle-ready', '--inplace', "--renumber-residues", '--delete-anisou', pdb]
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print(f'ERROR: rna tools failed: on {pdb}\n{out.decode()}\n{err.decode()}')

    # drop any duplicated atoms
    # ppdb.df['ATOM'].drop_duplicates(subset=ppdb.df['ATOM'].columns[:-1],inplace=True)

    # 147, 227 has first residue atom an "O" which is not recognized atom, for RNP rna-puzzle
    # is not correcting
    # for f in ${puzzle}/${puzzle}TS???_?.pdb; do sed -i '/^ATOM      1  O     G/d' $f; done
    # problems it corrects in same way:
    # some group 035, 054, 125 have cysteines N4 as a O in the element column
    # this crashes programs, so correctly change to a N


def center_map_and_native(emmap, pdb):
    system(f'cp {emmap} {emmap.rsplit(".")[0]}_uncentered.mrc')
    system(f'cp {pdb} {pdb.rsplit(".")[0]}_uncentered.mrc')
    with mrcfile.open(emmap, 'r+') as mrc:
        x_, y_, z_ = mrc.header.origin.x.item(), mrc.header.origin.y.item(), mrc.header.origin.z.item()
        mrc.header.origin.x = 0
        mrc.header.origin.y = 0
        mrc.header.origin.z = 0
    ppdb = PandasPdb().read_pdb(pdb)
    ppdb.df['ATOM'].x_coord = ppdb.df['ATOM'].x_coord - x_
    ppdb.df['ATOM'].y_coord = ppdb.df['ATOM'].y_coord - y_
    ppdb.df['ATOM'].z_coord = ppdb.df['ATOM'].z_coord - z_
    ppdb.to_pdb(path=pdb)


def write_and_sbatch_scoring(pdbs, N, prefix, name, base_sbatch, out_file_prefix, native, usalign_exec, EM, phenix_location, chimerax_exec=None, emmap=None, resolution=None, threshold=None,
                             sbatch=True):
    seperate_pdbs_to_run_groups(pdbs, N, prefix)
    for i in range(N):
        system(f'cp {base_sbatch} {prefix}_{i}.sbatch')
        system(f'sed -i s/XXX/{name}_{i}/g {prefix}_{i}.sbatch')
        with open(f'{prefix}_{i}.sbatch', 'a') as f:
            #f.write('\nconda activate casp_rna_em \n')
            if EM:
                f.write(f'python {path.dirname(__file__)}/casp_rna_score_all_from_file.py -f {prefix}_{i}.txt --out_file_prefix {out_file_prefix}_{i} --native {native} --usalign_exec {usalign_exec} --chimerax_exec {chimerax_exec} --EM --emmap {emmap} --resolution {resolution} --threshold {threshold} --phenix_location {phenix_location} \n')
            else:
                f.write(f'python {path.dirname(__file__)}/casp_rna_score_all_from_file.py -f {prefix}_{i}.txt --out_file_prefix {out_file_prefix}_{i} --native {native} --usalign_exec {usalign_exec} --phenix_location {phenix_location} \n')
        if sbatch:
            system(f'sbatch {prefix}_{i}.sbatch')
        else:
            system(f'chmod +x {prefix}_{i}.sbatch')
            system(f'./{prefix}_{i}.sbatch &')


def seperate_pdbs_to_run_groups(pdbs, N, prefix):
    num = (len(pdbs) + 1) // N
    for i in range(N):
        with open(f'{prefix}_{i}.txt', 'w') as f:
            for pdb in pdbs[i * num:min((i + 1) * num, len(pdbs))]:
                f.write(f'{pdb}\n')


def concat_all_result_files(files, out_file):
    dfs = []
    for f in files:
        dfs.append(pd.read_csv(f))
    dfs = pd.concat(dfs)
    dfs.to_csv(out_file, index=False)


METRICS = {"clashscore": "min", "global_rmsd": "min", "cc_mask": "max", "cc_peaks": "max",
           "tempy_mi": "max", "tempy_smoc": "max", "ai": "max",
           "gdt_ts": "max", "tm_score": "max", "lddt": "max", "inf_all": "max"}


def reduce_df(df, score_to_choice_best=None, static_columns=["target", "gr_code", "model"], metric_dict=METRICS):
    ''' 
    Reduce df to "best" score per X 
        (vs multiple scores per X if there are multiple Y)
    Eg for conformations:
    If score_to_choice_best is None then, for each score, the best value over all coformations is selected
    If it is specified then first the best conformation according to that score is selected,
    then all scores are taken from the conformation.
    '''
    columns_combining = [x for x in df.columns if x not in list(metric_dict.keys()) + static_columns]
    print("combing the following columns:", columns_combining)

    if score_to_choice_best is not None:
        if metric_dict[score_to_choice_best] == "max":
            best_index = df.groupby(static_columns)[score_to_choice_best].idxmax().to_numpy(copy=True)
        elif metric_dict[score_to_choice_best] == "min":
            best_index = df.groupby(static_columns)[score_to_choice_best].idxmin().to_numpy(copy=True)
        return df.iloc[best_index].drop(columns_combining, axis=1)
    else:
        return df.groupby(static_columns).agg(metric_dict).reset_index()


def calc_z(arr, raw):
    '''
    Calculate Z of second array with mean and std of the first
    # arr is values to calculate mean and std of
    # raw is values to calculate z with calculated mean and std
    '''
    u = np.nanmean(arr)
    s = np.nanstd(arr)
    z = (raw - u) / s
    return z


def get_zscore(arr, negative=False, threshold=-2, mask=None):
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


def get_weighted_sum_z(data, weights, prefix):
    final_z = None
    for score, weight in weights.items():
        if final_z is None:
            final_z = weight * data[prefix + score].to_numpy()
        else:
            final_z += weight * data[prefix + score].to_numpy()
    return final_z


def get_group_score(df, agg="sum>0", score="Z_rna"):
    if agg == "sum>0":
        return df.groupby("gr_code")[score].apply(lambda col: col[col > 0].sum()).reset_index()
    elif agg == "mean":
        return df.groupby("gr_code")[score].mean().reset_index()
    elif agg == "sum":
        return df.groupby("gr_code")[score].sum().reset_index()
