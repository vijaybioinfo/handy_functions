#!/usr/bin/python

####################################
# Neighborhood Components Analysis #
####################################

## This scripts runs NCA for a new t-SNE/UMAP visualisation. Before running this
## you need Loom data.

import optparse
parser = optparse.OptionParser()
parser.add_option("-i", "--mycellsf", dest="mycellsf", default='none',
                  help="Object: the loom path and name.",
                  metavar="OBJECT")
parser.add_option("-c", "--clusters", dest="clusters",
                  help="Cluster column.",
                  metavar="CLUSTERS")
parser.add_option("-n", "--num_NCA", dest="num_NCA", default=10,
                  help="Number of NCs to use for tSNE/UMAP.",
                  metavar="NC")
parser.add_option("-t", "--tsne", dest="tsne", action="store_true", default=True,
                  help="Whether to perform the t-SNE reduction.",
                  metavar="TSNE")
parser.add_option("-u", "--umap", dest="umap", action="store_true", default=True,
                  help="Whether to perform the UMAP reduction.",
                  metavar="UMAP")
parser.add_option("-m", "--metric", dest="metric", default='euclidean',
                  help="Distance metric [default: %default].",
                  metavar="METRIC")
parser.add_option("-s", "--suffix", dest="suffix", default='',
                  help="Suffix to add at the end of the output files. "
                  "If it contains 'scaled', it uses the scaled data in the object.",
                  metavar="SUFFIX")
group = optparse.OptionGroup(parser, "Dangerous Options",
                    "Caution: use these options at your own risk.  "
                    "It is believed that some of them bite.")
group.add_option("-e", "--state", dest="state", default=27,
                  help="Random seed for consistency.",
                  metavar="STATE")
parser.add_option_group(group)

(options, args) = parser.parse_args()
num_NCA = int(options.num_NCA)
state = int(options.state)
# if len(args) != 1:
#     parser.error("incorrect number of arguments")

# options.mycellsf = '/home/ciro/large/asthma_pjs/results/clust_seurat/resting_cd4_15p_sng_nomacro_x7n9/clustering/zetInfo/clustCells20PCs_30Ks_0.06667JD.loom'
# options.clusters = 'RNA_snn_res_0_4'

print("\n")
print("###### NCA analysis ######")
print("Loading libraries")
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

# import inspect
# inspect.getargspec(os.listdir)

print("\n")
print("File:", options.mycellsf)
fname = 'nca' + str(num_NCA) + options.suffix
if True:
    os.chdir(os.path.dirname(options.mycellsf))
    print("Saving to:", os.getcwd())

## Loading data
if not 'data' in locals():
    print("Loading data:", datetime.datetime.now())
    data = sc.read_loom(options.mycellsf) # data = sc.read_10x_h5(mycellsf10x)
    print("Finished:", datetime.datetime.now())

data
hvg_genes = data.var.index[data.var.vst_variable.to_numpy() > 0]

## Find if an NCA with more NCs than needed is present
import re
myncas = re.findall('nca[0-9]{1,}' + options.suffix + '.csv', str(" ".join(os.listdir())))
mynca = 'none'
for i in myncas:
    print(i)
    if num_NCA < int(str(''.join(re.findall('[0-9]{1,}', str(list(i)))))):
        mynca = i
        break

## Running NCA
# should be the scaled one: data.X = data.layers["scale"], though data.X worked better for CD8
if 'scaled' in options.suffix:
    print("Using scaled data")
    X = data.layers["scale_data"].toarray()
    X = X.T[~np.any(np.isnan(X.T), axis=1)].T
else:
    print("Using all data")
    X = data.X.toarray()

X.shape
y = data.obs[options.clusters]

compute_dim = [options.umap, options.tsne]
# Compute NCA
if os.path.exists(fname + '.csv') or os.path.exists(mynca):
    print("NCA precomputed")
    precomp_df = pd.read_csv(filepath_or_buffer = fname + '.csv', index_col = 0)
    print(precomp_df.columns.values)
    data.obsm["X_nca"] = precomp_df.filter(regex = '^X_nca.$', axis = 1).to_numpy()
    if not os.path.exists(mynca):
        print("Checking for t-SNE and UMAP")
        compute_dim[0] = not 'X_nca_tsne1' in precomp_df.columns.values
        compute_dim[1] = not 'X_nca_umap1' in precomp_df.columns.values
else:
    print("NCA start:", datetime.datetime.now())
    nca = NeighborhoodComponentsAnalysis(n_components = num_NCA, random_state = state)
    data.obsm["X_nca"] = nca.fit_transform(X, y)
    print("NCA endin:", datetime.datetime.now())

num_NCA_use = list(range(num_NCA))
# if not len(num_NCA_use) == num_NCA:
#     fname = 'nca' + str(len(num_NCA_use)) + options.suffix
#     print("New file name:", fname)

# UMAP with NCA projection
if compute_dim[1]:
    print("UMAP start:", datetime.datetime.now())
    X = data.obsm["X_nca"][:, num_NCA_use]
    umap = UMAP(n_components = 2, n_neighbors = 15, metric = options.metric, random_state = state)
    data.obsm["X_nca_umap"] = umap.fit_transform(X)
    print("UMAP endin:", datetime.datetime.now())
else:
    print("UMAP precomputed")
    data.obsm["X_nca_umap"] = precomp_df.filter(regex = '^X_nca_umap', axis = 1).to_numpy()

# t-SNE with NCA projection
if compute_dim[0]:
    print("t-SNE start:", datetime.datetime.now())
    X = data.obsm["X_nca"][:, num_NCA_use]
    tsne = TSNE(n_components = 2, perplexity = 30, metric = options.metric, random_state = state)
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
            c = data.obs[options.clusters].astype(int),
            marker = ".",
            cmap = 'tab20'
        )
        fontP = FontProperties()
        fontP.set_size('small')
        legend1 = ax.legend(
            *scatter.legend_elements(),
            title = options.clusters#,
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
