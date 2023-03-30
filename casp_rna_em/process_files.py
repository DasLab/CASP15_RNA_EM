from glob import glob
import subprocess as sp
from os import system,path,remove
from biopandas.pdb import PandasPdb
from casp_rna_em.run_metric_programs import run_usalign


def prepare_pdbs(pdbs,ignore_no_extension=False,fill_occupancy=True,fill_bfactor=True):
    pdbs = glob(pdbs)

    # only look at file with pdb extension or no extension in which case add pdb
    if ignore_no_extension:
        pdbs = [x for x in pdbs if pdbs.rsplit(".",1)[-1]=="pdb"]
    else:
        pdbs_a = [x for x in pdbs if x.rsplit(".",1)[-1]=="pdb"]
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
        clean_pdb(pdb,fill_occupancy=fill_occupancy,fill_bfactor=fill_bfactor)
    return pdbs

def clean_pdb(pdb,fill_occupancy=True,fill_bfactor=True): 
    # use rna tools to prepare pdbs, action completed:
    # remove any modeled H so they are not scored
    # renames chains
    command = ['rna_pdb_tools.py','--get-rnapuzzle-ready','--inplace',"--renumber-residues",pdb]
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = p.communicate()
    if p.returncode:
        print(f'ERROR: rna tools failed: on {pdb}\n{out.decode()}\n{err.decode()}')

    if fill_occupancy or fill_bfactor:
        ppdb = PandasPdb().read_pdb(pdb)
        # make all atoms have full occupanncy
        if fill_occupancy:
            ppdb.df['ATOM'].occupancy = 1
        # remove bfactor information
        if fill_bfactor:
            ppdb.df['ATOM'].b_factor = 0
        ppdb.to_pdb(path=pdb)

    # drop any duplicated atoms
    #ppdb.df['ATOM'].drop_duplicates(subset=ppdb.df['ATOM'].columns[:-1],inplace=True)

    # 147, 227 has first residue atom an "O" which is not recognized atom, for RNP rna-puzzle
    # is not correcting
    #for f in ${puzzle}/${puzzle}TS???_?.pdb; do sed -i '/^ATOM      1  O     G/d' $f; done
    # problems it corrects in same way:
    # some group 035, 054, 125 have cysteines N4 as a O in the element column
    # this crashes programs, so correctly change to a N



# TODO center mrc?
