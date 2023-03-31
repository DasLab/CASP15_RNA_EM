# CASP15 RNA assesment - cryoEM
For the CASP15 RNA category, we aimed to assess predictions by directly comparing to expiremental data. In the respository is all scripts and code used to obtain and analyse scores. All map-to-model scores used are molecule agnostic and thus can be applied to RNA as here or other molecules, however, some of the pre-processing and other scores herein are RNA specific.

## Requirements

### Programs used

#### Phenix

Install phenix from https://phenix-online.org/download

#### TEMPy

`pip install BioTEMPy==2.0.0`

#### Qscore

`git clone https://github.com/gregdp/mapq.git`

#### Chimera

## Pipeline

## Analysis

### Choosing best score
Code can be found in `results/rna_results_utils.py`

For CASP15 we choice to compare all submitted models from each group (5) to the target and take the best score. If there was multiple conformation we compared to all conformations and took the best score. See `reduce_df`.

### Z score based ranking
For CASP15, we ranked based on Zscore. For each target and metric a Zscore was calculated using the best score, as described above. Outliers (Z<2) were removed for the final mean and standard deviation calculations. See `get_zscore`. The weighted Zscore was then calculated for each target and group, see `get_weighted_sum_z`. The final ranking was decided after, for each group, summing Z-scores for every target, ignoring any Z<0, see `get_group_score`.

### Visualizations
Visualizations used in the manuscript can be found in `results/rna_results.ipynb`. Further, chimerax sessions to create the figures can be found in this folder, and the maps and models related to these files can be found [here](https://drive.google.com/file/d/1b6ZZXznF2zvvdVO1PuSEmxmOLnh2yJQY/view?usp=share_link).


## References
CASP15 papers are in preperation
