## Simulation of Spatial Transcriptomics spots from single-cell reference

<!-- Note: This folder provides a collection of scripts we used to simulated data but the scripts need to be edited to be used on other platforms (contains hard-coded paths for our HPC). -->
This repo provides a collection of scripts used to generate simulated spatial transcriptomics data as a mixture of single-cell transcriptomics profiles (adapting code and model from [Andersson et al. 2019](https://www.biorxiv.org/content/10.1101/2019.12.13.874495v1)). Parameters are chosen to simulate the characteristics of 10X Genomics Visium chips.
<!-- 
#### Contents

- `ST_simulation.py` functions to simulate spots
- `split_sc.py` script to split mouse brain snRNA-seq reference into generation dataset (single cells used to make the synthetic spots) and validation dataset (single cells used as a reference to train the location model)
- `assemble_ST.py` script to generate synthetic ST spots from the generation dataset
 -->
### Run simulation 

Initial input: 

- AnnData object of raw counts per single-cells (saved as `h5ad` file)
- A table of cell type annotations per cell that we want to deconvolve (saved as `csv` file)

**Step 1: Split single-cell dataset:** we split the cells in the single-cell dataset in a 'generation set', that will be used to simulate the ST spots, and a 'validation' set, that will be used to train the deconvolution models that we want to benchmark. From the command line:

```
python split_sc.py <counts_h5ad> <annotation_csv> --annotation_col annotation_1  --out_dir <output_directory>
```

Output: generation and validation count matrices and cell type annotations are saved as `pickle` files, with a random seed identifying the split. 

**Step 2: Build design matrix**: in this step we define which cell types are (A) low/high density and (B) Uniformly present in all the spots or localized in few spots (regional). To generate synthetic spots with ~10 cells per spot (as seen with nuclear segmentation on Visium spots) we reccommend setting the mean number of cells per spot per cell type < 5.

```
n_spots=100
seed=$(ls labels_generation* | sed 's/.*_//' | sed 's/.p//')
python ST_simulation/assemble_design.py \
    $seed \
  --tot_spots $n_spots --mean_high 3 --mean_low 1 \
  --out_dir <output_directory>
```

Output: `synthetic_ST_seed${seed}_${assemble_id}_design.csv` contains the design used for the simulation:

| **Column**  | **Data**                                                                                         |
|-------------|--------------------------------------------------------------------------------------------------|
| uniform     | is the cell type uniformly located across spots (1) or localized in a small subset of spots (0)  |
| density     | is the cell type present in a spot at low density (1) or high density (0)                        |
| nspots      | total number of spots in which the cell type is located                                          |
| mean_ncells | mean number of cells per spot                                                                    |

**Step 3: Assemble cell type composition per spot:** based on the design matrix, we define the cell type composition of each spot i.e. how many cells per cell type
are in each spot. An assemble ID is used to identify the assembly (we assemble many composition matrices with the same design).
```
id=1
python cell2location/pycell2location/ST_simulation/assemble_composition.py \
    $seed \
    --tot_spots $n_spots --assemble_id $id
```

Output: `synthetic_ST_seed${seed}_${assemble_id}_composition.csv` contains the number of cells per cell type in each spot, for benchmarking deconvolution models.

**(Step 4) Assemble simulated ST spots**
```
python assemble_st.py ${seed} --assemble_id $id
```

Output:

- `synthetic_ST_seed${seed}_${assemble_id}_counts.csv` contains the count matrix for the simulated ST spots
- `synthetic_ST_seed${seed}_${assemble_id}_umis.csv` contains the number of UMIs per cell type in each spot, for benchmarking deconvolution methods that model number of UMIs


### Speeding up the simulation

The current implementation is not optimized for speed, it takes ~ 2 minutes to assemble 100 spots. At the moment my suggestion to simulate thousands of spots is to assemble the design matrix once (step 1 above), then run steps 2 and 3 many times using wrapper 
```
run_simulation2.sh <seed> <n_spots> <id> 
```
then merge in one object
```
python merge_synthetic_ST.py . $seed
```


