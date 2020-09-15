#!/usr/bin/python

####################################
# Neighborhood Components Analysis #
####################################

# This scripts runs NCA for a new t-SNE/UMAP visualisation.
# You first need to generate a loom format file

OBJECT=/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/clustCells20PCs_30Ks_0.06667JD

/share/apps/R/3.5/bin/Rscript /home/ciro/scripts/clustering/nca_reduction_seurat2loom.R -m ${OBJECT}.RData

conda activate clustpy

python /home/ciro/scripts/clustering/nca_reduction.py --mycellsf=${OBJECT}.loom --num_NCA=15 --clusters=RNA_snn_res_0_4
