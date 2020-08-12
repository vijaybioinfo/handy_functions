#!/usr/bin/python

####################################
# Neighborhood Components Analysis #
####################################

# ## This scripts runs NCA for a new t-SNE/UMAP visualisation. Before running this
# ## you need Loom data. To do it so with a Seurat object:
# setwd('/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd8_15p_sng_nomacro_x8/clustering/zetInfo/')
# mycellsf = "clustCells25PCs_30Ks_0.06667JD.RData"
# setwd('/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x6/clustering/zetInfo/')
# mycellsf = "clustCells19PCs_30Ks_0.06667JD.RData"
# setwd('/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x6/subcluster0n1/clustering/zetInfo/')
# mycellsf = "clustCells22PCs_30Ks_0.06667JD.RData"
# setwd('/home/ciro/large/allergy_u19/results/clust_seurat/015_Asp_15_H9_30p/clustering/zetInfo/')
# mycellsf = "clustCells18PCs_30Ks_0.06667JD.RData"
#
# .libPaths('~/R/newer_packs_library/3.5/')
# load(mycellsf)
# mycells@reductions <- list()
# mycells@graphs <- list()
# mycellsloom <- Seurat::as.loom(mycells, filename = sub(".rdata", ".loom", mycellsf, ignore.case = TRUE))
# mycellsloom$close_all()

from sklearn.neighbors import NeighborhoodComponentsAnalysis
from sklearn.manifold import TSNE
from umap import UMAP
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import os
import numpy as np
import datetime
import inspect

state = 27
metric = "euclidean"
num_TSNE = 2
num_NCA = 10
suffix = ''

mycellsf = '/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd8_15p_sng_nomacro_x8/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.loom'
clusters = 'RNA_snn_res_0_6'
suffix = '_scaled'

# mycellsf10x = '/home/ciro/large/hayley/raw/sever_asthma3/COUNTS_hg19/AGGR/CD4US_filt_mapped/outs/filtered_feature_bc_matrix.h5'
mycellsf = '/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x6/clustering/zetInfo/clustCells19PCs_30Ks_0.06667JD.loom'
clusters = 'RNA_snn_res_0_4'
num_NCA = 5
num_NCA = 3

mycellsf = '/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x6/subcluster0n1/clustering/zetInfo/clustCells22PCs_30Ks_0.06667JD.loom'
clusters = 'RNA_snn_res_0_4'
num_NCA = 5
num_NCA = 3
num_NCA = 15
suffix = '_scaled'

mycellsf = '/home/ciro/large/allergy_u19/results/clust_seurat/015_Asp_15_H9_30p/clustering/zetInfo/clustCells18PCs_30Ks_0.06667JD.loom'
clusters = 'RNA_snn_res_0_4'
num_NCA = 5

fname = 'nca' + str(num_NCA) + suffix
if True:
    os.chdir(os.path.dirname(mycellsf))
    os.getcwd()
    print("Saving to:", fname)

## Loading data
if not 'data' in locals():
    print("Loading data:", datetime.datetime.now())
    data = sc.read_loom(mycellsf) # data = sc.read_10x_h5(mycellsf10x)
    print("Finished:", datetime.datetime.now())

data
hvg_genes = data.var.index[data.var.vst_variable.to_numpy() > 0]

## Running NCA
# should be the scaled one: data.X = data.layers["scale"], though data.X worked better for CD8
if 'scaled' in suffix:
    print("Using scaled data")
    X = data.layers["scale_data"].toarray()
    X = X.T[~np.any(np.isnan(X.T), axis=1)].T
else:
    print("Using all data")
    X = data.X.toarray()

X.shape
y = data.obs[clusters]

compute_dim = [True, True]
# Compute NCA
if os.path.exists(fname + '.csv'):
    print("NCA precomputed")
    precomp_df = pd.read_csv(filepath_or_buffer = fname + '.csv', index_col = 0)
    print(precomp_df.columns.values)
    data.obsm["X_nca"] = precomp_df.filter(regex = '^X_nca.$', axis = 1).to_numpy()
    compute_dim[0] = not 'X_nca_tsne1' in precomp_df.columns.values
    compute_dim[1] = not 'X_nca_umap1' in precomp_df.columns.values
    # df_c = pd.concat([data.obsm.to_df(), precomp_df], axis = 1); df_c.columns
else:
    print("NCA start:", datetime.datetime.now())
    nca = NeighborhoodComponentsAnalysis(n_components = num_NCA, random_state = state)
    data.obsm["X_nca"] = nca.fit_transform(X, y)
    print("NCA endin:", datetime.datetime.now())

num_NCA_use = list(range(num_NCA)) # num_NCA
# if not len(num_NCA_use) == num_NCA:
#     fname = 'nca' + str(len(num_NCA_use)) + suffix
#     print("New file name:", fname)

# UMAP with NCA projection
if compute_dim[1]:
    print("UMAP start:", datetime.datetime.now())
    X = data.obsm["X_nca"][:, num_NCA_use]
    umap = UMAP(n_components = 2, n_neighbors = 15, metric = metric, random_state = state)
    data.obsm["X_nca_umap"] = umap.fit_transform(X)
    print("UMAP endin:", datetime.datetime.now())
else:
    print("UMAP precomputed")
    data.obsm["X_nca_umap"] = precomp_df.filter(regex = '^X_nca_umap', axis = 1).to_numpy()

# t-SNE with NCA projection
if compute_dim[0]:
    print("t-SNE start:", datetime.datetime.now())
    X = data.obsm["X_nca"][:, num_NCA_use]
    tsne = TSNE(n_components = 2, perplexity = 30, metric = metric, random_state = state)
    data.obsm["X_nca_tsne"] = tsne.fit_transform(X)
    print("t-SNE endin:", datetime.datetime.now())
else:
    print("t-SNE precomputed")
    data.obsm["X_nca_tsne"] = precomp_df.filter(regex = '^X_nca_tsne', axis = 1).to_numpy()

# Results table
df = data.obsm.to_df().filter(regex = '^X_nca', axis = 1)
print(df.columns.values)
df.to_csv(path_or_buf = fname + '.csv')

# Plotting
if True:
    for redu in ['X_nca_tsne', 'X_nca_umap']:
        print(redu)
        fig, ax = plt.subplots(figsize=(7,7))
        scatter = ax.scatter(
            x = data.obsm[redu][:,0],
            y = data.obsm[redu][:,1],
            c = data.obs[clusters].astype(int),
            cmap = 'tab20'
        )
        fontP = FontProperties()
        fontP.set_size('small')
        legend1 = ax.legend(
            *scatter.legend_elements(),
            title = clusters#,
            #prop = fontP,
            #loc = 'upper center', bbox_to_anchor=(0.5,-0.1)
        )
        ax.add_artist(legend1)
        # ax.set_axis_off()
        plt.tight_layout() # bbox_extra_artists=(legend1), bbox_inches='tight' in savefig
        plt.savefig(fname + redu.replace("X_nca", "") + '.pdf')
        plt.close()

print("Finished:", datetime.datetime.now())
print('END')
