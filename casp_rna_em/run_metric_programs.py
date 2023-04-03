import subprocess as sp
import numpy as np
import pandas as pd
from os import system, path
from tqdm import tqdm


def score_all_from_file(file, out_file_prefix, native, usalign_location='', chimerax_location='', EM=True, emmap=None, resolution=None, threshold=None, phenix_location=""):
    pdbs = []
    with open(file, 'r') as f:
        for line in f.readlines():
            pdbs.append(line.strip())
    score_all(pdbs, out_file_prefix=out_file_prefix, native=native, usalign_location=usalign_location,
              chimerax_location=chimerax_location, EM=EM, emmap=emmap, resolution=resolution, threshold=threshold, phenix_location=phenix_location)


def score_all(pdbs, out_file_prefix, native, usalign_location='', chimerax_location='', EM=True, emmap=None, resolution=None, threshold=None, phenix_location=""):

    single_scores = {}
    per_residue_scores = {}
    per_threshold_scores = {}

    for pdb in tqdm(pdbs):
        # print('scoring:', pdb)
        score_dict = {}
        score_dict.update(run_phenix_clashscore(pdb, phenix_location=phenix_location))
        score_dict.update(run_phenix_rna_validate(pdb, phenix_location=phenix_location))
        score_dict.update(run_lga(pdb, native, atom='C4,'))
        score_dict.update(run_usalign(pdb, native, usalign_location=usalign_location))

        if EM:
            # TODO add skip USalign if that pdb already exists.
            # TODO phenix option for docking
            # TODO add all scores (fsc, Q), maybe option of score to not do?
            score_dict.update(dock_pdb_usalign_fitmap(pdb, emmap, native, threshold, usalign_location, chimerax_location))
            docked_pdb = f'{pdb.rsplit(".",1)[0]}_{emmap.rsplit(".",1)[0].rsplit("/",1)[-1]}_usalignfitmapDOCKED.pdb'
            score_dict.update(run_phenix_cc(docked_pdb, emmap, resolution=resolution, phenix_location=phenix_location))
            score_dict.update(run_atomic_inclusion(docked_pdb, emmap, threshold, chimerax_location))
            score_dict.update(run_tempy(docked_pdb, emmap, resolution))

        single_score = {key: score_dict[key] for key in score_dict if "per_residue" not in key and "per_threshold" not in key}
        per_residue_score = {key: score_dict[key] for key in score_dict if "per_residue" in key}
        per_threshold_score = {key: score_dict[key] for key in score_dict if "per_threshold" in key}
        model_name = pdb.rsplit(".", 1)[0].rsplit('/', 1)[-1]
        single_scores[model_name] = single_score
        per_residue_scores[model_name] = per_residue_score
        per_threshold_scores[model_name] = per_threshold_score

    df = pd.DataFrame(single_scores).transpose()
    df.index.name = 'pdb'
    df = df.reset_index()
    df['native'] = native.rsplit('/', 1)[-1]
    if EM:
        df['emmap'] = emmap.rsplit('/', 1)[-1]
    df['target'] = df.pdb.apply(lambda x: x.split('TS')[0])
    df['gr_code'] = df.pdb.apply(lambda x: x.split('TS')[1].split("_")[0])
    df['model'] = df.pdb.str[-1]
    df.to_csv(f"{out_file_prefix}_scores.csv", index=False)

    per_res = pd.DataFrame(per_residue_scores).transpose()
    per_res.index.name = 'model'
    # TODO need check num residues correct?
    if len(per_res.columns) > 0:
        per_res = per_res.explode(list(per_res.columns))
    per_res.to_csv(f"{out_file_prefix}_per_residue.csv")

    per_thr = pd.DataFrame(per_threshold_scores).transpose()
    per_thr.index.name = 'model'
    # TODO need write threhsolds
    if len(per_thr.columns) > 0:
        per_thr = per_thr.explode(list(per_thr.columns))
    per_thr.to_csv(f"{out_file_prefix}_per_threshold.csv")

###############################################################################
# Individual score runs
###############################################################################


def run_phenix_clashscore(pdb, phenix_location='', nuclear=True, keep_hydrogens=True, result_file=None):
    # nuclear=False for xray
    # keep_hydrogens=False to remove any hydrogens and add using reduce
    if result_file is None:
        result_file = f'{pdb.rsplit(".")[0]}_CLASHSCORE.out'
    command = [f'{phenix_location}phenix.clashscore', pdb, f'nuclear={nuclear}', f'keep_hydrogens={keep_hydrogens}']
    out_file = open(result_file, 'w')
    p = sp.Popen(command, stdout=out_file, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print(f'ERROR: phenix.clashscore failed: on {pdb}\n{err.decode()}')
        result_file = None
    out_file.close()
    return parse_phenix_clashscore(result_file)


def run_phenix_rna_validate(pdb, phenix_location='', result_file=None):
    if result_file is None:
        result_file = f'{pdb.rsplit(".")[0]}_RNAVAL.out'
    command = [f'{phenix_location}phenix.rna_validate', pdb]
    out_file = open(result_file, 'w')
    p = sp.Popen(command, stdout=out_file, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print(f'ERROR: phenix.rna_validate failed: on {pdb}\n{err.decode()}')
        result_file = None
    out_file.close()
    return parse_rna_validate(result_file)


def run_dssr():
    # TODO
    # ${dssr_exec} -i=${f} --more --non-pair --prefix=${f} > ${f}_dssr.out
    return None


def run_phenix_cc(pdb, emmap, resolution, phenix_location='', result_file=None):
    if result_file is None:
        result_file = f'{pdb.rsplit(".",1)[0]}_phenixCC.out'
    command = [f'{phenix_location}phenix.model_vs_map', pdb, emmap, f'resolution={resolution}']
    out_file = open(result_file, 'w')
    p = sp.Popen(command, stdout=out_file, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print(f'ERROR: phenix.model_vs_map failed: on {pdb}\n{err.decode()}')
        result_file = None
    out_file.close()
    return parse_phenix_cc(result_file)


def run_phenix_fsc(pdb, emmap, phenix_location='', result_file=None):
    if result_file is None:
        result_file = f'{pdb.rsplit(".",1)[0]}_phenixFSC.out'
    command = [f'{phenix_location}phenix.mtriage', emmap, pdb]
    out_file = open(result_file, 'w')
    p = sp.Popen(command, stdout=out_file, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print(f'ERROR: phenix.mtriage failed: on {pdb}\n{err.decode()}')
        result_file = None
    out_file.close()
    return parse_phenix_fsc(result_file)


def run_atomic_inclusion(pdb, emmap, threshold, chimerax_location='', result_file=None):
    prefix = f'{pdb.rsplit(".",1)[0]}'
    if result_file is None:
        result_file = f'{prefix}_AI.out'

    system(f'{chimerax_location}ChimeraX --nogui --script "{path.dirname(__file__)}/chimerax_atom_inclusion.py {emmap} {pdb} {threshold} {prefix}" > {result_file}')

    return parse_atomic_inclusion(result_file, prefix, threshold)


def run_Qscore(pdb, emmap, mapq_script, chimera_location, resolution, result_file=None, timeout=3600, threads=16):
    # TODO likely including new script?
    score_dict = {"q": np.nan, "q_bb": np.nan, "q_base": np.nan,
                  'q_by_residue': [], 'g_bb_by_residue': [], 'q_base_by_residue': []}
    command = ['timeout', timeout, 'python', mapq_script, chimera_location, emmap, pdb, f'res={resolution}', f'np={threads}']
    ''' timeout 3600 $python_exec ${mapq_location}mapq_cmd.py $chimera_location ${cryoem_map}_CENTERED.mrc ${f}_DOCKED.pdb res=$resolution
    BACKBONE_ATOMS = ['P','OP1','OP2',\"O5'\",\"C5'\",\"C4'\",\"O4'\",\"C3'\",\"O3'\",\"C2'\",\"O2'\",\"C1'\"]
    df_heavy_backbone = df_heavy[df_heavy.atom_name.isin(BACKBONE_ATOMS)]
    df_heavy_base = df_heavy[~df_heavy.atom_name.isin(BACKBONE_ATOMS)]
    df_heavy.b_factor[df_heavy.b_factor<0] = 0
    df_heavy_backbone.b_factor[df_heavy_backbone.b_factor<0] = 0
    df_heavy_base.b_factor[df_heavy_base.b_factor<0] = 0
    print('Average Q pos:',df_heavy.b_factor.mean())
    per_nuc_Q = df_heavy.rename(columns={'b_factor':'Q'}).groupby('residue_number').mean().Q
    Note this is a fork of the main mapq repository enabeling multiple q-score calculation in parralel and tabular output in the command-line tool.
    `git clone https://github.com/rkretsch/mapq.git`

    Qscore requires chimera. Download [chimera](https://www.cgl.ucsf.edu/chimera/download.html), agreeing to non-commerical liscence, install according to instructions, and keep note of location. For scoring Q-score will need to specify `chimera_location='chimera-location'`
    '''
    return score_dict


def run_tempy(pdb, emmap, resolution):
    score_dict = {'tempy_ccc': np.nan, 'tempy_mi': np.nan,
                  'tempy_lsf': np.nan, 'tempy_env': np.nan,
                  'tempy_sccc': np.nan, 'tempy_smoc': np.nan,
                  'tempy_smoc_1': np.nan, 'tempy_smoc_per_residue': [],
                  'tempy_smoc_1_per_residue': []}

    from TEMPy.protein.structure_parser import PDBParser  # TEMPy class for parsing pdb files
    from TEMPy.maps.map_parser import MapParser           # TEMPy class for parsing mrc files
    from TEMPy.protein.structure_blurrer import StructureBlurrer
    import TEMPy.protein.scoring_functions as scf         # TEMPy class containing all scoring functions

    blurrer = StructureBlurrer()
    sc = scf.ScoringFunctions()
    em_map = MapParser.readMRC(emmap)

    model = PDBParser.read_PDB_file(pdb, pdb)

    sim_map = blurrer.gaussian_blur_real_space(model, resolution, em_map)

    score_dict['tempy_ccc'] = sc.CCC(sim_map, em_map)

    score_dict['tempy_mi'] = sc.MI(sim_map, em_map)
    score_dict['tempy_lsf'] = sc.LSF(sim_map, em_map)

    mass = model.get_prot_mass_from_atoms()
    prim_bound = em_map.get_primary_boundary(molWeight=mass, low=em_map.min(), high=em_map.max(), vol_factor=1.21)
    score_dict['tempy_env'] = sc.envelope_score(em_map, prim_bound, model)
    # print(f"tempy primary boundary {prim_bound}")
    score_dict['tempy_sccc'] = sc.SCCC(em_map, resolution, 0.356, model, model)
    SMOC = sc.SMOC(em_map, resolution, model, win=11)  # default
    SMOC_1 = sc.SMOC(em_map, resolution, model, win=1)  # default
    N, total, total1 = 0, 0, 0
    chains = list(SMOC[0].keys())
    for chain in chains:
        for res, score in SMOC[0][chain].items():
            score_dict['tempy_smoc_per_residue'].append(score)
            score_dict['tempy_smoc_1_per_residue'].append(SMOC_1[0][chain][res])
    score_dict['tempy_smoc'] = sum(score_dict['tempy_smoc_per_residue']) / len(score_dict['tempy_smoc_per_residue'])
    score_dict['tempy_smoc1'] = sum(score_dict['tempy_smoc_1_per_residue']) / len(score_dict['tempy_smoc_1_per_residue'])
    return score_dict


def run_lga(pdb, native, atom='C4,', result_file=None):
    if result_file is None:
        result_file = f'{pdb.rsplit(".",1)[0]}_{native.rsplit(".",1)[0].rsplit("/",1)[-1]}_LGA.out'
    # system(f'ulimit -s unlimited; lga')
    command = [f'runlga.mol_mol.pl', pdb, native, '-3', f'-atom:{atom}', '-gdc', '-stral', '-o2']
    out_file = open(result_file, 'w')
    p = sp.Popen(command, stdout=out_file, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print(f'ERROR: LGA failed: on {pdb}\n{err.decode()}')
        result_file = None
    out_file.close()
    result_file = f'TMP/{pdb.rsplit("/",1)[-1]}.{native.rsplit("/",1)[-1]}.stral'
    return parse_lga(result_file)


def run_usalign(pdb, native, usalign_location='', output_pdb=None, result_file=None):
    if result_file is None:
        result_file = f'{pdb.rsplit(".",1)[0]}_{native.rsplit(".",1)[0].rsplit("/",1)[-1]}_USALIGN.out'
    if output_pdb is None:
        output_pdb = f'{pdb.rsplit(".",1)[0]}_{native.rsplit(".",1)[0].rsplit("/",1)[-1]}_USALIGN'
    command = [f'{usalign_location}USalign', pdb, native, '-outfmt', '2', '-o', output_pdb]
    out_file = open(result_file, 'w')
    p = sp.Popen(command, stdout=out_file, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print(f'ERROR: usalign failed: on {pdb}\n{err.decode()}')
        result_file = None
    out_file.close()
    return parse_usalign(result_file)


###############################################################################
# Parsing of output files
###############################################################################

def parse_phenix_clashscore(result_file):
    score_dict = {'clashscore': np.nan}
    if result_file is not None:
        with open(result_file, 'r') as f:
            lines = f.readlines()
            for row in lines:
                if 'clashscore' == row[:10]:
                    score_dict['clashscore'] = row.split()[2]
    return score_dict


def parse_rna_validate(result_file):
    score_dict = {'bond_outlier': np.nan, 'angle_outlier': np.nan,
                  'pucker_outlier': np.nan, 'suite_outlier': np.nan,
                  'avg_suitness': np.nan}
    if result_file is not None:
        with open(result_file, 'r') as f:
            lines = f.readlines()
            for row in lines:
                bond_outlier = row.find('bond outliers present')
                perfect_bond = row.find('All bonds within 4.0 sigma of ideal values')
                angle_outlier = row.find('angle outliers present')
                perfect_angle = row.find('All angles within 4.0 sigma of ideal values')
                pucker_outlier = row.find('pucker outliers present')
                perfect_pucker = row.find('All puckers have reasonable geometry')
                suite_outlier = row.find('suite outliers present')
                perfect_suite = row.find('All RNA torsion suites are reasonable')
                avg_suitness = row.find('Average suiteness:')
                if bond_outlier != -1:
                    bond_outlier = row.split()[0]
                    score_dict['bond_outlier'] = float(bond_outlier.split('/')[0]) / float(bond_outlier.split('/')[1])
                elif perfect_bond != -1:
                    score_dict['bond_outlier'] = 0
                if angle_outlier != -1:
                    angle_outlier = row.split()[0]
                    score_dict['angle_outlier'] = float(angle_outlier.split('/')[0]) / float(angle_outlier.split('/')[1])
                elif perfect_angle != -1:
                    score_dict['angle_outlier'] = 0
                if pucker_outlier != -1:
                    pucker_outlier = row.split()[0]
                    score_dict['pucker_outlier'] = float(pucker_outlier.split('/')[0]) / float(pucker_outlier.split('/')[1])
                elif perfect_pucker != -1:
                    score_dict['pucker_outlier'] = 0
                if suite_outlier != -1:
                    suite_outlier = row.split()[0]
                    score_dict['suite_outlier'] = float(suite_outlier.split('/')[0]) / float(suite_outlier.split('/')[1])
                elif perfect_suite != -1:
                    score_dict['suite_outlier'] = 0
                if avg_suitness != -1:
                    score_dict['avg_suitness'] = row.split()[2]
    return score_dict


def parse_phenix_fsc(result_file):
    score_dict = {'fsc_0_mask': np.nan, 'fsc_0_unmask': np.nan,
                  'fsc_143_mask': np.nan, 'fsc_143_unmask': np.nan,
                  'fsc_50_mask': np.nan, 'fsc_50_unmask': np.nan, }
    if result_file is not None:
        with open(result_file, 'r') as f:
            lines = f.readlines()
            for row in lines:
                fsc_0_present = row.find('FSC(map,model map)=0 ')
                fsc_143_present = row.find('FSC(map,model map)=0.143')
                fsc_50_present = row.find('FSC(map,model map)=0.5')
                if fsc_0_present != -1:
                    score_dict['fsc_0_mask'] = row.split()[3]
                    score_dict['fsc_0_unmask'] = row.split()[4]
                if fsc_143_present != -1:
                    score_dict['fsc_143_mask'] = row.split()[3]
                    score_dict['fsc_143_unmask'] = row.split()[4]
                if fsc_50_present != -1:
                    score_dict['fsc_50_mask'] = row.split()[3]
                    score_dict['fsc_50_unmask'] = row.split()[4]
    return score_dict


def parse_phenix_cc(result_file):
    score_dict = {'cc_volume': np.nan, 'cc_mask': np.nan,
                  'cc_peaks': np.nan, 'cc_per_residue': [], 'per_residue_cc_residue': []}
    if result_file is not None:
        with open(result_file, 'r') as f:
            lines = f.readlines()
            for row in lines:
                if row[:11] == '  CC_volume':
                    score_dict['cc_volume'] = row.split()[1]
                if row[:9] == '  CC_mask':
                    score_dict['cc_mask'] = row.split()[2]
                if row[:10] == '  CC_peaks':
                    score_dict['cc_peaks'] = row.split()[2]
                if row.find('chain id') != -1:
                    score_dict['cc_per_residue'].append(row.split()[6])
                    score_dict['per_residue_cc_residue'].append(row.split()[4])
    return score_dict


def parse_lga(result_file):
    score_dict = {'gdt_ts': np.nan, 'rmsd_lga': np.nan}
    if result_file is not None:
        with open(result_file, 'r') as f:
            lines = f.readlines()
            for row in lines:
                if "SUMMARY(GDT)" == row[:12]:
                    score_dict['gdt_ts'] = row.split()[6]
                    score_dict['rmsd_lga'] = row.split()[5]
    return score_dict


def parse_usalign(result_file):
    score_dict = {'tm': np.nan, 'rmsd_usalign': np.nan}
    # TM1 TM2 RMSD    ID1 ID2 IDali   L1  L2  Lali
    if result_file is not None:
        with open(result_file, 'r') as f:
            row = f.readlines()[-1]
            score_dict['tm'] = row.split()[2]
            score_dict['rmsd_usalign'] = row.split()[4]
    return score_dict


def parse_atomic_inclusion(result_file, prefix, threshold):
    score_dict = {"ai": np.nan, "ai_bb": np.nan, "ai_base": np.nan,
                  "do": np.nan, "do_bb": np.nan, "do_base": np.nan,
                  'ai_per_threshold': [], 'ai_bb_per_threshold': [], 'ai_base_per_threshold': [],
                  'do_per_threshold': [], 'do_bb_per_threshold': [], 'do_base_per_threshold': [],
                  'per_threshold_threshold': [],
                  'ai_per_residue': [], 'ai_bb_per_residue': [], 'ai_base_per_residue': [],
                  'per_residue_ai_residue': []}

    with open(result_file, 'r') as f:
        lines = f.readlines()
        for row in lines:
            AI_present = row.find('AI at ')
            AI_bb_present = row.find('AI_backbone at ')
            AI_base_present = row.find('AI_base at ')
            DO_present = row.find('DO at ')
            DO_bb_present = row.find('DO_backbone at ')
            DO_base_present = row.find('DO_base at ')
            if AI_present != -1:
                score_dict['ai'] = row.split()[3]
            if AI_bb_present != -1:
                score_dict['ai_bb'] = row.split()[3]
            if AI_base_present != -1:
                score_dict['ai_base'] = row.split()[3]
            if DO_present != -1:
                score_dict['do'] = row.split()[3]
            if DO_bb_present != -1:
                score_dict['do_bb'] = row.split()[3]
            if DO_base_present != -1:
                score_dict['do_base'] = row.split()[3]
    per_res = pd.read_csv(f'{prefix}_threshold{threshold}_AI_per_residue.csv', delim_whitespace=True)
    score_dict['ai_per_residue'] = per_res.atom_inclusion.to_list()
    score_dict['ai_bb_per_residue'] = per_res.atom_inclusion_backbone.to_list()
    score_dict['ai_base_per_residue'] = per_res.atom_inclusion_base.to_list()
    score_dict['per_residue_ai_residue'] = per_res.residue.to_list()

    per_thr = pd.read_csv(f'{prefix}_DO_AI_per_threshold.csv', delim_whitespace=True)
    score_dict['ai_per_threshold'] = per_thr.atom_inclusion.to_list()
    score_dict['ai_bb_per_threshold'] = per_thr.atom_inclusion_backbone.to_list()
    score_dict['ai_base_per_threshold'] = per_thr.atom_inclusion_base.to_list()
    score_dict['do_per_threshold'] = per_thr.density_occupancy.to_list()
    score_dict['do_bb_per_threshold'] = per_thr.density_occupancy_backbone.to_list()
    score_dict['do_base_per_threshold'] = per_thr.density_occupancy_base.to_list()
    score_dict['per_threshold_threshold'] = per_thr.threshold.to_list()
    return score_dict


###############################################################################
# Docking models
###############################################################################

def dock_pdb_phenix_iterative(pdb, emmap, resolution, phenix_location='', output_pdb=None):
    if output_pdb is None:
        output_pdb = f'{pdb.rsplit(".",1)[0]}_{emmap.rsplit(".",1)[0].rsplit("/",1)[-1]}_phenixDOCKED.pdb'
    if path.isfile(output_pdb):
        remove(output_pdb)
    dock_pdb_phenix(pdb, emmap, resolution, phenix_location=phenix_location, output_pdb=output_pdb, skip_if_low_cc=True, min_cc=0.4)
    if path.isfile(output_pdb):
        return {'docked': 'phenix'}
    dock_pdb_phenix(pdb, emmap, resolution, phenix_location=phenix_location, output_pdb=output_pdb, skip_if_low_cc=False, min_cc=0.4)
    if path.isfile(output_pdb):
        return {'docked': 'phenix'}
    dock_pdb_phenix(pdb, emmap, resolution, phenix_location=phenix_location, output_pdb=output_pdb, skip_if_low_cc=False, min_cc=0.2)
    if path.isfile(output_pdb):
        return {'docked': 'phenix'}
    dock_pdb_phenix(pdb, emmap, resolution, phenix_location=phenix_location, output_pdb=output_pdb, skip_if_low_cc=False, min_cc=0.1)
    if path.isfile(output_pdb):
        return {'docked': 'phenix'}
    dock_pdb_phenix(pdb, emmap, resolution, phenix_location=phenix_location, output_pdb=output_pdb, skip_if_low_cc=False, min_cc=0.05)
    if path.isfile(output_pdb):
        return {'docked': 'phenix'}
    dock_pdb_phenix(pdb, emmap, resolution, phenix_location=phenix_location, output_pdb=output_pdb, skip_if_low_cc=False, min_cc=0.02)
    if path.isfile(output_pdb):
        return {'docked': 'phenix'}
    dock_pdb_phenix(pdb, emmap, resolution, phenix_location=phenix_location, output_pdb=output_pdb, skip_if_low_cc=False, min_cc=0.01)
    if path.isfile(output_pdb):
        return {'docked': 'phenix'}
    dock_pdb_phenix(pdb, emmap, resolution, phenix_location=phenix_location, output_pdb=output_pdb, skip_if_low_cc=False, min_cc=0.01)
    if path.isfile(output_pdb):
        return {'docked': 'failed-phenix'}


def dock_pdb_phenix(pdb, emmap, resolution, phenix_location='', output_pdb=None, result_file=None, skip_if_low_cc=True, min_cc=0.4):
    if output_pdb is None:
        output_pdb = f'{pdb.rsplit(".",1)[0]}_{emmap.rsplit(".",1)[0].rsplit("/",1)[-1]}_phenixDOCKED.pdb'
    if result_file is None:
        result_file = f'{pdb.rsplit(".",1)[0]}_{emmap.rsplit(".",1)[0].rsplit("/",1)[-1]}_phenixDOCKED.out'
    command = [f'{phenix_location}phenix.dock_in_map', emmap, pdb, f'resolution={resolution}', f'pdb_out={output_pdb}', f'skip_if_low_cc={skip_if_low_cc}', f'min_cc={min_cc}']
    out_file = open(result_file, 'w')
    p = sp.Popen(command, stdout=out_file, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print(f'ERROR: phenix dock failed: on {pdb}\n{err.decode()}')
        result_file = None
    out_file.close()


def dock_pdb_usalign_fitmap(pdb, emmap, native, threshold, usalign_location='', chimerax_location='', output_pdb=None, result_file=None):
    if output_pdb is None:
        output_pdb = f'{pdb.rsplit(".",1)[0]}_{emmap.rsplit(".",1)[0].rsplit("/",1)[-1]}_usalignfitmapDOCKED.pdb'
    if result_file is None:
        result_file = f'{pdb.rsplit(".",1)[0]}_{emmap.rsplit(".",1)[0].rsplit("/",1)[-1]}_usalignfitmapDOCKED.out'
    run_usalign(pdb, native, usalign_location)
    usalign_output_pdb = f'{pdb.rsplit(".",1)[0]}_{native.rsplit(".",1)[0].rsplit("/",1)[-1]}_USALIGN.pdb'
    system(f'{chimerax_location}ChimeraX --nogui --script "{path.dirname(__file__)}/chimerax_fit.py {usalign_output_pdb} {emmap} {output_pdb} {threshold}" > {result_file}')
    return {'docked': 'usalign-fitmap'}

# run_phenix_rna_validate('R1107TS147_1.pdb',phenix_location='/home/rachael/phenix_1.20/phenix-1.20.1-4487/build/bin/')
# run_phenix_rna_validate('R1108TS248_5.pdb',phenix_location='/home/rachael/phenix_1.20/phenix-1.20.1-4487/build/bin/')
# parse_rna_validate('R1107TS232_1_seg0.pdb_RNAVAL.out')
# print(run_phenix_clashscore('R1107TS147_1.pdb',phenix_location='/home/rachael/phenix_1.20/phenix-1.20.1-4487/build/bin/'))
# print(run_phenix_fsc('R1149TS128_1.pdb_DOCKED.pdb','../results/models/R11149.mrc_CENTERED.mrc',phenix_location='/home/rachael/phenix_1.20/phenix-1.20.1-4487/build/bin/'))
# print(run_phenix_cc('R1149TS128_1.pdb_DOCKED.pdb','../results/models/R11149.mrc_CENTERED.mrc',resolution=5,phenix_location='/home/rachael/phenix_1.20/phenix-1.20.1-4487/build/bin/'))
# print(run_usalign('R1107TS147_1.pdb','../R1107-CPEB3_ribozyme_human/D_1292119758_model-annotate_P1human.pdb','/home/groups/rhiju/philpham/casp15/us-align/USalign'))
# print(run_lga('R1107TS147_1.pdb','../R1107-CPEB3_ribozyme_human/D_1292119758_model-annotate_P1human.pdb',lga_location='/home/groups/rhiju/rkretsch/LGA/LGA_package_src/'))


# print(run_tempy('R1149TS128_1.pdb_DOCKED.pdb','../results/models/R11149.mrc_CENTERED.mrc',5))
# atomic inclusion
# qscore
