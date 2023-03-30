from casp_rna_em.process_files import prepare_pdbs
from casp_rna_em.run_metric_programs import score_all
from os import system
import pandas as pd


###############################################################################
# TODO
###############################################################################
# sbatch off proccesses, run 1 group, wait until all done, and combine output files
# Qscore function and update on changes to mapq code for parralization + sig figs
# version of phenix used
# list phenix exectibles needed
# ddsr
# rna tools install 
# phenix random edits
# /home/rachael/phenix_1.20/phenix-1.20.1-4487/modules/cctbx_project/mmtbx/monomer_library/pdb_interpretation.py
# max_reasonable_bond_distance change to 100 in master_params_str
# pip install -e git+http://github.com/mmagnus/rna-tools.git#egg=rna-tools
# test refactoring on all targets from CASP15
# export PYTHONPATH=$PYTHONPATH:/home/rachael/Desktop/Das_Lab/CASP/CASP15_RNA_EM
# finish score all
# map preprocessing?
# add in segmenetation code
# add instructions for all installs
# links to data sources
###############################################################################


usalign_exec = '/home/rachael/USalign/USalign'
chimerax_exec = '/home/rachael/chimerax-1.5/bin/ChimeraX'
phenix_location = '/home/rachael/phenix_1.20/phenix-1.20.1-4487/build/bin/'

tar_file = 'R1136.tar.gz'
models_folder = 'R1136'
native = '../results/models/R1136v1_Apta_FRET-Bound_rsr15-withAt2OP1.pdb_DOCKED.pdb'
emmap = '../../fits/R1136-Apta_FRET_EM/bound/native/R1136v1_structure_22.mrc'
threshold = 0.019
resolution = 4

system(f'cp {native} R1136_native.pdb')
native = 'R1136_native.pdb'
prepare_pdbs(native)

#system(f'rm -r {models_folder}')
#system(f'tar -xf {tar_file}')
#pdbs = prepare_pdbs(f'{models_folder}/*')
# temp
from glob import glob
pdbs = glob('R1136/R1136TS???_?.pdb')

score_all(pdbs[:2],"test",emmap,native,resolution,threshold,usalign_exec,chimerax_exec,EM=True,phenix_location=phenix_location)
'''
from process_files import prepare_pdbs
from run_metric_programs import run_atomic_inclusion,dock_pdb_usalign_fitmap
from os import system

native = '../results/models/R1136v1_Apta_FRET-Bound_rsr15-withAt2OP1.pdb_DOCKED.pdb'
models_folder = 'R1136'
phenix_location = '/home/rachael/phenix_1.20/phenix-1.20.1-4487/build/bin/'
resolution = 4
emmap = '../../fits/R1136-Apta_FRET_EM/bound/native/R1136v1_structure_22.mrc'
usalign_exec = '/home/rachael/USalign/USalign'
chimerax_exec = '/home/rachael/chimerax-1.5/bin/ChimeraX'
threshold = 0.019
pdb = 'R1136/R1136TS232_2.pdb'

#system(f'cp {native} R1136_native.pdb')
#prepare_pdbs('R1136_native.pdb')

#dock_pdb_usalign_fitmap(pdb,emmap,'R1136_native.pdb',threshold,usalign_exec,chimerax_exec)
#dock_pdb_phenix_iterative(pdb,emmap,resolution,phenix_location=phenix_location)
docked_pdb = 'R1136/R1136TS232_2_R1136v1_structure_22_usalignfitmapDOCKED.pdb'
print(run_atomic_inclusion(docked_pdb,emmap,threshold,chimerax_exec))
'''