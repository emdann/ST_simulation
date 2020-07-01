### Make design of simulated ST datasets from single-cell data
import argparse
import pickle
import numpy as np 
import pandas as pd

parser = argparse.ArgumentParser()
# parser.add_argument('lbl_gen_file', type=str,
#                     help='path to label generation pickle file')
# parser.add_argument('cnt_gen_file', type=str,
#                     help='path to label generation pickle file')
parser.add_argument('seed', type=int,
                    help='random seed of split')
parser.add_argument('--tot_spots', dest='tot_spots', type=int,
                    default=1000,
                    help='Total number of spots to simulate')
parser.add_argument('--mean_high', dest='mean_high', type=float,
                    default=3.0,
                    help='Mean cell density for high-density cell types')
parser.add_argument('--mean_low', dest='mean_low', type=float,
                    default=1.0,
                    help='Mean cell density for low-density cell types')
parser.add_argument('--percent_uniform', dest='percent_uniform', type=float,
                    default=80,
                    help='Sparsity of uniform cell types (% non-zero spots of total spots)')
parser.add_argument('--percent_sparse', dest='percent_sparse', type=float,
                    default=15,
                    help='Sparsity of sparse cell types (% non-zero spots of total spots)')
parser.add_argument('--annotation_col', dest='anno_col', type=str,
                    default="annotation_1",
                    help='Name of column to use in annotation file (default: annotation_1)')
parser.add_argument('--out_dir', dest='out_dir', type=str,
                    default='/nfs/team283/ed6/simulation/lowdens_synthetic_ST_fewcells/',
                    help='Output directory')
parser.add_argument('--assemble_id', dest='assemble_id', type=int,
                    default=1,
                    help='ID of ST assembly')

args = parser.parse_args()

# lbl_gen_file = args.lbl_gen_file
# count_gen_file = args.cnt_gen_file
seed = args.seed
tot_spots = args.tot_spots
mean_high = args.mean_high
mean_low = args.mean_low
percent_uniform = args.percent_uniform
percent_sparse = args.percent_sparse
out_dir = args.out_dir
assemble_id = args.assemble_id
anno_col = args.anno_col

### Load input data ### 
lbl_gen_file = out_dir + "labels_generation_" + str(seed) + ".p"
count_gen_file = out_dir + "counts_generation_" + str(seed) + ".p"

lbl_generation = pickle.load(open(lbl_gen_file, "rb"))
cnt_generation = pickle.load(open(count_gen_file, "rb"))

uni_labels = lbl_generation[anno_col].unique()
labels = lbl_generation
cnt = cnt_generation

### Define uniform VS sparse cell types (w more sparse)
uniform_ct = np.random.choice([0, 1], size=len(uni_labels), p=[0.8, 0.2])

#### Define low VS high density cell types (w more low density)
uni_low = np.random.choice([0, 1], size=len(uni_labels[uniform_ct == 1]),  p=[0.3, 0.7])
reg_low = np.random.choice([0, 1], size=len(uni_labels[uniform_ct == 0]),  p=[0.3, 0.7])

design_df = pd.DataFrame({'uniform': uniform_ct}, index=uni_labels)

design_df['density'] = np.nan
design_df.loc[design_df.index[design_df.uniform == 1], 'density'] = uni_low
design_df.loc[design_df.index[design_df.uniform == 0], 'density'] = reg_low

### Generate no of spots per cell type 
# Uniform ~ 60% of spots, sparse ~ 5% of spots
mean_unif = round((tot_spots / 100) * percent_uniform)
mean_sparse = round((tot_spots / 100) * percent_sparse)
sigma_unif = np.sqrt(mean_unif / 0.05)
sigma_sparse = np.sqrt(mean_sparse / 0.05)

shape_unif = mean_unif ** 2 / sigma_unif ** 2
scale_unif = sigma_unif ** 2 / mean_unif
shape_sparse = mean_sparse ** 2 / sigma_sparse ** 2
scale_sparse = sigma_sparse ** 2 / mean_sparse

unif_nspots = np.round(np.random.gamma(shape=shape_unif, scale=scale_unif, size=sum(design_df.uniform == 1)))
sparse_nspots = np.round(np.random.gamma(shape=shape_sparse, scale=scale_sparse, size=sum(design_df.uniform == 0)))
# if samples n spots is greater than total number of spots trim to the total
if (unif_nspots > tot_spots).sum() >= 1:
    unif_nspots[unif_nspots > tot_spots] = tot_spots
if (sparse_nspots > tot_spots).sum() >= 1:
    sparse_nspots[sparse_nspots > tot_spots] = tot_spots


design_df['nspots'] = np.nan
design_df.loc[design_df.index[design_df.uniform == 1], 'nspots'] = unif_nspots
design_df.loc[design_df.index[design_df.uniform == 0], 'nspots'] = sparse_nspots

### Generate avg density per spot per cell type
sigma_low = np.sqrt(mean_low / 2)
sigma_high = np.sqrt(mean_high / 2)

shape_low = mean_low ** 2 / sigma_low ** 2
scale_low = sigma_low ** 2 / mean_low
shape_high = mean_high ** 2 / sigma_high ** 2
scale_high = sigma_high ** 2 / mean_high

low_ncells_mean = np.random.gamma(shape=shape_low, scale=scale_low, size=sum(design_df.density == 1))
high_ncells_mean = np.random.gamma(shape=shape_high, scale=scale_high, size=sum(design_df.density == 0))

design_df['mean_ncells'] = np.nan
design_df.loc[design_df.index[design_df.density == 1], 'mean_ncells'] = low_ncells_mean
design_df.loc[design_df.index[design_df.density == 0], 'mean_ncells'] = high_ncells_mean

out_name = out_dir + "synthetic_ST_seed" + lbl_gen_file.split("_")[-1].rstrip(".p") + "_" + "design" + ".csv"
design_df.to_csv(out_name, sep=",", index=True, header=True)
