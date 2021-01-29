### SPLIT SINGLE-CELL DATASET IN GENERATION AND VALIDATION SET ###

import argparse
import pickle
import random
import anndata
import scanpy as sc
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('h5ad_file', type=str,
                    help='path to h5ad file with raw single-cell data')
parser.add_argument('annotation_file', type=str,
                    help='path to csv file with cell annotations')
parser.add_argument('--annotation_col', dest='anno_col', type=str,
                    default="annotation_1",
                    help='Name of column to use in annotation file (default: annotation_1)')
parser.add_argument('--out_dir', dest='out_dir', type=str,
                    default="/nfs/team283/ed6/simulation/lowdens_synthetic_ST_fewcells/",
                    help='Output directory')
args = parser.parse_args()

# adata_file = "/nfs/team205/vk7/sanger_projects/cell2location/notebooks/data/mouse_viseum_snrna/rawdata/all_cells_20200625.h5ad"
# annotation_file = "/nfs/team205/vk7/sanger_projects/cell2location/notebooks/results/mouse_viseum_snrna/snRNA_annotation_20200229.csv"
# anno_col = "annotation_1"
# out_dir = "/nfs/team283/ed6/simulation/lowdens_synthetic_ST_fewcells/"

adata_file = args.h5ad_file
annotation_file = args.annotation_file
anno_col = args.anno_col
out_dir = args.out_dir

### Load input single-cell data and annotations ###

adata_raw = sc.read_h5ad(adata_file)

## Cell type annotations
labels = pd.read_csv(annotation_file, index_col=0)

# # Select genes used for clustering
# adata_raw = adata_raw[:, [x for x in adata_raw.var_names if x in adata_snrna.var_names]]

adata_df = pd.DataFrame(adata_raw.X.T.toarray(), columns=adata_raw.obs_names, index=adata_raw.var_names)
adata_df = adata_df.T
adata_df.index.name = "cell"

### Subset to cells with label ###
adata_df = adata_df.loc[labels.index, :]

### Split generation and validation set ###

sc_cnt = adata_df
sc_lbl = pd.DataFrame(labels[anno_col])

# match count and label data
inter = sc_cnt.index.intersection(sc_lbl.index)

sc_lbl = sc_lbl.loc[inter, :]
sc_cnt = sc_cnt.loc[inter, :]

labels = sc_lbl.iloc[:, 0].values

# get unique labels
uni_labs, uni_counts = np.unique(labels, return_counts=True)

# only keep types with more than 50 cells
keep_types = uni_counts > 40
keep_cells = np.isin(labels, uni_labs[keep_types])

labels = labels[keep_cells]
sc_cnt = sc_cnt.iloc[keep_cells, :]
sc_lbl = sc_lbl.iloc[keep_cells, :]

uni_labs, uni_counts = np.unique(labels, return_counts=True)
n_types = uni_labs.shape[0]

seeds = random.sample(range(1000), 3)

for seed in seeds:
    random.seed(seed)
    print("Seed " + str(seed))
    # get member indices for each set
    idx_generation = []
    idx_validation = []
    for z in range(n_types):
        tmp_idx = np.where(labels == uni_labs[z])[0]
        n_generation = int(round(tmp_idx.shape[0] / 2))
        smp_gen = random.sample(list(tmp_idx), k=n_generation)
        smp_val = tmp_idx[np.isin(tmp_idx, smp_gen, invert=True)]
        idx_generation += smp_gen
        idx_validation += smp_val.tolist()
    idx_generation.sort()
    idx_validation.sort()
    # make sure no members overlap between sets
    assert len(set(idx_generation).intersection(set(idx_validation))) == 0, \
        "validation and generation set are not orthogonal"
    # assemble sets from indices
    cnt_validation = sc_cnt.iloc[idx_validation, :]
    cnt_generation = sc_cnt.iloc[idx_generation, :]
    lbl_validation = sc_lbl.iloc[idx_validation, :]
    lbl_generation = sc_lbl.iloc[idx_generation, :]
    pickle.dump(lbl_generation,
                open(out_dir + "labels_generation_" + str(seed) + ".p", "wb"))
    pickle.dump(cnt_generation,
                open(out_dir + "counts_generation_" + str(seed) + ".p", "wb"))
    pickle.dump(lbl_validation,
                open(out_dir + "labels_validation_" + str(seed) + ".p", "wb"))
    pickle.dump(cnt_validation,
                open(out_dir + "counts_validation_" + str(seed) + ".p", "wb"))
