# chimerax --nogui --script 'atom_inclusion_chimerax.py CENTERED.mrc FIT.pdb rec_thresholds(can input as many as desired)'
# saves CENTERED.mrc_FIT.pdb_DO_AI.csv
# calculate density occupancy and atomic inclusion with heavy atoms and vdw radius
# for a range of threshold including ant reccomened threhsold
# will also calcualte per residue AI for thresholds inputted

import sys
import numpy as np
from itertools import product
from chimerax.core.commands import run

# load mrc and set range of thresholds to check
mrc = run(session, f'open {sys.argv[1]}')
max_value = mrc[0].data.matrix().max()
min_value = mrc[0].data.matrix().min()
std_value = np.std(mrc[0].data.matrix())
thresholds = [x/1000 for x in range(round(min_value*1000),round(max_value*1000),round(std_value*1000/10))]
specific_threshold = float(sys.argv[3])
thresholds.append(specific_threshold)
thresholds.sort()

# load pdb and get atom positions
pdb = run(session, f'open {sys.argv[2]}')
atom_coords = pdb[0].atoms.coords
elements = pdb[0].atoms.element_names
BACKBONE_ATOMS = ['P','OP1','OP2',"O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"]
backbone = np.isin(pdb[0].atoms.names,BACKBONE_ATOMS)
resnumber = pdb[0].atoms.residues.numbers

# get map values at atomic positions
atom_values = run(session, f'measure mapvalues #1 atoms #2')

# remove hydrogens and get backbone index
H_index = np.where(elements == 'H')
heavy_elements = np.delete(elements,H_index)
heavy_atom_coords = np.delete(atom_coords,H_index,axis=0)
heavy_atom_values = np.delete(atom_values[0],H_index)
heavy_resnumber = np.delete(resnumber,H_index)

# seperate base and backbone
heavy_backbone = np.delete(backbone,H_index)
num_atom = heavy_elements.shape[0]
num_backbone = heavy_backbone.sum()
num_base = num_atom-num_backbone
backbone_index = np.where(heavy_backbone)[0]
heavy_backbone_values = heavy_atom_values[heavy_backbone]
heavy_base_values = heavy_atom_values[~heavy_backbone]
heavy_backbone_resnumber = heavy_resnumber[heavy_backbone]
heavy_base_resnumber = heavy_resnumber[~heavy_backbone]

# for inputted thresholds find the AI
residues = list(set(heavy_resnumber))
AI_res = []
AI_res_backbone = []
AI_res_base = []
for res in residues:
    res_values = heavy_atom_values[heavy_resnumber==res]
    res_values_backbone = heavy_backbone_values[heavy_backbone_resnumber==res]
    res_values_base = heavy_base_values[heavy_base_resnumber==res]

    num_res_threshold = np.count_nonzero(res_values > specific_threshold)
    num_res_backbone_threshold = np.count_nonzero(res_values_backbone > specific_threshold)
    num_res_base_threshold = np.count_nonzero(res_values_base > specific_threshold)

    total_res = res_values.shape[0]
    total_backbone = res_values_backbone.shape[0]
    total_base = res_values_base.shape[0]

    AI_res.append(num_res_threshold/total_res)
    AI_res_backbone.append(num_res_backbone_threshold/total_backbone)
    AI_res_base.append(num_res_base_threshold/total_base)
np.savetxt(f'{sys.argv[4]}_threshold{specific_threshold}_AI_per_residue.csv',np.array([residues,AI_res,AI_res_backbone,AI_res_base]).T,delimiter=' ',header='residue atom_inclusion atom_inclusion_backbone atom_inclusion_base',comments='',fmt='%15.15f')

# get VDW radius
vdw_radii = {'C':1.7,'N':1.55,'O':1.52,'P':1.8}
vdw_dists = np.vectorize(vdw_radii.get)(heavy_elements)

# get vertex values and identify which within vdw of atoms
all_vertex_values = mrc[0].data.matrix().flatten()
mrc_data_shape = mrc[0].data.matrix().shape
near_atom_vertex_values = []
near_atom_vertex_values_backbone = []
near_atom_vertex_values_base = []
index_already_added = []
for i,a_coord in enumerate(heavy_atom_coords):
    min_search_coord = a_coord-6 # 6A arbitrary far enough search
    max_search_coord = a_coord+6
    min_search_ijk = mrc[0].data.xyz_to_ijk(min_search_coord).round().astype(int)
    max_search_ijk = mrc[0].data.xyz_to_ijk(max_search_coord).round().astype(int)
    search_ijk = np.array(list(product(range(min_search_ijk[0],max_search_ijk[0]+1),range(min_search_ijk[1],max_search_ijk[1]+1),range(min_search_ijk[2],max_search_ijk[2]+1))))
    search_xyz = mrc[0].data.ijk_to_xyz(search_ijk)
    dist_vector = a_coord - search_xyz
    dist = np.linalg.norm(dist_vector,axis=1)
    dist_min_vdw = dist-vdw_dists[i]
    dist_within_atom = np.where(dist_min_vdw<0)[0]
    for index in dist_within_atom:
        # save index to not repeat
        if tuple(search_ijk[index]) not in index_already_added:
            maps_value = mrc[0].data.matrix()[max(min(mrc_data_shape[0]-1,search_ijk[index][0]),0),max(min(mrc_data_shape[1]-1,search_ijk[index][1]),0),max(0,min(mrc_data_shape[2]-1,search_ijk[index][2]))]
            near_atom_vertex_values.append(maps_value)
            index_already_added.append(tuple(search_ijk[index]))
            if i in backbone_index:
                near_atom_vertex_values_backbone.append(maps_value)
            else:
                near_atom_vertex_values_base.append(maps_value)

# calculate the density occupancy and atomic inclusion at each threshold
DOs = []
DOs_backbone = []
DOs_base = []
AIs  = []
AIs_backbone = []
AIs_base = []
for threshold in thresholds:
    num_total_vertex_threshold = np.count_nonzero(all_vertex_values > threshold)
    num_near_atom_vertex_threshold = np.count_nonzero(np.array(near_atom_vertex_values) > threshold)
    num_near_atom_vertex_threshold_backbone = np.count_nonzero(np.array(near_atom_vertex_values_backbone) > threshold)
    num_near_atom_vertex_threshold_base = np.count_nonzero(np.array(near_atom_vertex_values_base) > threshold)
    num_atom_threshold = np.count_nonzero(heavy_atom_values > threshold)
    num_backbone_threshold = np.count_nonzero(heavy_backbone_values > threshold)
    num_base_threshold = np.count_nonzero(heavy_base_values > threshold)
    DOs.append(num_near_atom_vertex_threshold/num_total_vertex_threshold)
    DOs_backbone.append(num_near_atom_vertex_threshold_backbone/num_total_vertex_threshold)
    DOs_base.append(num_near_atom_vertex_threshold_base/num_total_vertex_threshold)
    AIs.append(num_atom_threshold/num_atom)
    AIs_backbone.append(num_backbone_threshold/num_backbone)
    AIs_base.append(num_base_threshold/num_base)
    if threshold == specific_threshold:
        print(f'AI at {threshold}: {AIs[-1]}')
        print(f'AI_backbone at {threshold}: {AIs_backbone[-1]}')
        print(f'AI_base at {threshold}: {AIs_base[-1]}')
        print(f'DO at {threshold}: {DOs[-1]}')
        print(f'DO_backbone at {threshold}: {DOs_backbone[-1]}')
        print(f'DO_base at {threshold}: {DOs_base[-1]}')

# save and exit
np.savetxt(f'{sys.argv[4]}_DO_AI_per_threshold.csv',np.array([thresholds,DOs,DOs_backbone,DOs_base,AIs,AIs_backbone,AIs_base]).T,delimiter=' ',header='threshold density_occupancy density_occupancy_backbone density_occupancy_base atom_inclusion atom_inclusion_backbone atom_inclusion_base',comments='',fmt='%15.15f')
run(session, 'exit')
