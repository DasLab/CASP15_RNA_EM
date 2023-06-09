# CASP15 RNA assesment - cryoEM
For the CASP15 RNA category, we aimed to assess predictions by directly comparing to expiremental data. In the respository is all scripts and code used to obtain and analyse scores. All map-to-model scores used are molecule agnostic and thus can be applied to RNA as here or other molecules, however, some of the pre-processing and other scores herein are RNA specific.

## Requirements

### Programs used

#### Python packages

To install all necessary python packages, in this repository run:

```
conda env create -f environment.yml
conda activate casp_rna_em
pip install .
```

This will also install [__TEMPy__](https://doi.org/10.1107/s2059798320014928) and [__rna-tools__](https://rna-tools.readthedocs.io/en/latest/).

To run this code always activate the appropriate conda environment before running.

`conda activate casp_rna_em`

#### Phenix

Obtain an [academic liscense](https://phenix-online.org/phenix_request/index.cgi) and [download](https://phenix-online.org/download/) phenix. Install according to instruction, for example for linux:
```
tar -xf phenix-installer-*.tar.gz
./install --prefix="desired_location"
```

You can either add `source desired_location/phenix-?.??.?-????/phenix_env.sh` to your `~/.bashrc` or in applicable function calls in this repository specify `phenix_location='desired_location/phenix-?.??.?-????/build/bin/'`.

The binaries used are `phenix.clashscore`, `phenix.rna_validate`, `phenix.model_vs_map`, `phenix.mtriage`, and `phenix.dock_in_map`.

Versions 1.18.2 and 1.20.1 (linux) were tested.

There are models whose geometries surpass the checks of phenix. They will crash phenix and not recieve phenix scores. This can be manually overwritten if desired by changing `max_reasonable_bond_distance` in `master_params_str` to `100` in `phenix_location/modules/cctbx_project/mmtbx/monomer_library/pdb_interpretation.py`

#### ChimeraX

[Download](https://www.cgl.ucsf.edu/chimerax/download.html), agreeing to non-commerical liscence, and install according to the instructions for your opperating system ChimeraX. You can either add `alias ChimeraX='chimerax-location/bin/Chimerax'` or in applicable function calls in this repository specify `chimerax_location='chimerax-location/bin/'`

#### Qscore

Qscores were calculated but not used in the final ranking because of the low resolution of the maps. 

_Section under construction_.

#### US-align

Simply install according to:

```
git clone https://github.com/pylelab/USalign.git
cd USalign
make
```
You can either add `export PATH=$PATH:/location/USalign` to your `~/.bashrc` or in applicable function calls in this repository specify `usalign_location='/location/USalign'`.


#### LGA

Request [LGA download](http://as2ts.proteinmodel.org/AS2TS/Download_code/).
```
tar -xf LGA_package_src.*.tar.gz
cd LGA_package_src/
```
Add to you `~/.bashrc`:
```
export PATH=$PATH:lga-location
ulimit -s unlimited; lga
``` 

## Pipeline

Examples of running the pipeline can be found in `example_run_scripts/` with a full run in `example_run.py` and then parralizable options compatible with a slurm computer cluster run with `example_run_para.py`. To combine results from these runs, use the file followed by `example_run_combine.py` once all the calculations are done.

## Analysis

Code can be found in `process_files.py`. With examples on how these steps were run `results/rna_results.ipynb`, `example_run_combine.py`.

### Choosing best score

For CASP15 we choice to compare all submitted models from each group (5) to the target and take the best score. If there was multiple conformation we compared to all conformations and took the best score. See `reduce_df`.

### Z score based ranking
For CASP15, we ranked based on Zscore. For each target and metric a Zscore was calculated using the best score, as described above. Outliers (Z<-2) were removed for the final mean and standard deviation calculations. See `get_zscore`. The weighted Zscore was then calculated for each target and group, see `get_weighted_sum_z`. The final ranking was decided, for each group, by summing Z-scores for every target, ignoring any Z<0, see `get_group_score`.

### Visualizations
Visualizations used in the manuscript can be found in `results/rna_results.ipynb`. Further, chimerax sessions to create the figures can be found in this folder, and the maps and models related to these files can be found [here](https://drive.google.com/file/d/1b6ZZXznF2zvvdVO1PuSEmxmOLnh2yJQY/view?usp=share_link).

### Analysis of R1149 R1156

For the CASP15 RNA expirementalist paper, R1149 and R1556 were analysed in detail including all models junction geometry and angle. This analysis is in the `R1149_R1156_analysis` folder.

## References
CASP15 papers are in preperation
