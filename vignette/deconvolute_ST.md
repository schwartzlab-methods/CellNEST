## Running CellNEST on Spatial Transcriptomics data after deconvolution
We integrate scRNA-seq data with Spatial Transcriptomics data for deconvoluting spots into single cells using [CytoSPACE](https://github.com/digitalcytometry/cytospace). 
We first import the necessary Python libraries.

```
import os 
import pandas as pd
import numpy as np
import scanpy as sc 
import anndata as ad 
```

CytoSPACE needs scRNA-seq gene expression file to be provided in a *.csv format with genes (rows) by cells (columns) format. 
For details, please see the official site of [CytoSPACE](https://github.com/digitalcytometry/cytospace#input-files)

We download the scRNA-seq for human lymph node from [here](https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad) and save under 'data/scRNAseq_lymph/' directory. We write the following function to convert the downloaded format into CytoSPACE required format:

```
# Written By
# Deisha Paliwal

# Write single-cell RNA-seq atlas to a CSV file in the format expected by CytoSPACE 
def cytospaceRef(
        adata: ad.AnnData, 
        outFolder: str, 
        ct_col: str
    ) -> None:
    adata.var_names_make_unique()
    # sample 10,000 cells at equal cell type proportions 
    total_cells = 10000
    cell_types = adata.obs[ct_col].value_counts()
    num_cell_types = len(cell_types)
    cells_per_type = total_cells // num_cell_types
    sampled_indices = []
    for cell_type, count in cell_types.items():
        cell_type_indices = adata.obs[adata.obs[ct_col] == cell_type].index
        if count <= cells_per_type:
            sampled_indices.extend(cell_type_indices)
        else:
            sampled_indices.extend(np.random.choice(cell_type_indices, cells_per_type, replace = False))
    adata_subset = adata[sampled_indices, :].copy()
    # write raw counts to a CSV file 
    counts = getRawCounts(adata_subset)
    counts_df = pd.DataFrame(counts.T, index = adata_subset.var_names, columns = adata_subset.obs_names)
    adata_subset.var_names_make_unique()
    counts_df.reset_index(inplace = True)
    counts_df.rename(columns={"index": "GENES"}, inplace = True)
    counts_df.to_csv(os.path.join(outFolder, f"sc_counts.csv"), index = False)
    # write cell type labels to a CSV file 
    cell_type_df = pd.DataFrame({
        "Cell IDs": adata_subset.obs_names,
        "CellType": adata_subset.obs[ct_col]
    })
    # remove special characters 
    cell_type_df["CellType"] = cell_type_df["CellType"].str.replace(r"[^a-zA-Z0-9 ]", " ", regex = True)
    cell_type_df.to_csv(os.path.join(outFolder, f"sc_cell_types.csv"), index = False)
```

Then we read the scRNA-seq file and call this function.
```
in_dir = 'data/scRNAseq_lymph/'
out_dir = 'data/scRNAseq_lymph/CytoSPACE_format/'
adata_sc = sc.read(os.path.join(in_dir, "Moffit_Combined.h5ad"))
cytospaceRef(adata_sc, out_dir, "max_moffitt")
```

Let us assume we have the Spatial Transcriptomics data for [human lymph node](https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_filtered_feature_bc_matrix.tar.gz) in a standard Visium format in this path 'data/V1_Human_Lymph_Node/'. 
Then we call CytoSPACE as follows:
```
 cytospace \
    --scRNA-path data/scRNAseq_lymph/CytoSPACE_format/sc_counts.csv \
    --cell-type-path data/scRNAseq_lymph/CytoSPACE_format/sc_cell_types.csv \
    --spaceranger-path data/V1_Human_Lymph_Node/V1_Human_Lymph_Node_filtered_feature_bc_matrix.tar.gz
```

We will use following two outputs of CytoSPACE:
1. assigned_locations.csv (the X and Y coordinates of each single cell mapped to ST spots)
2. assigned_expression (Containes barcodes.tsv, genes.tsv, and matrix.mtx representing the gene expression of the resulting assignments)

We move these two outputs to 'data/V1_Human_Lymph_Node/deconvolution/' directory.
Now, we will execute the following commands to run CellNEST model on this data.

Preprocess:
```
cellnest preprocess --data_name='lymph_deconvolution' --data_from='data/V1_Human_Lymph_Node/deconvolution/assigned_expression/' --tissue_position_file='data/V1_Human_Lymph_Node/deconvolution/assigned_locations.csv'
```

Model run: 
```
nohup cellnest run --data_name='lymph_deconvolution'  --num_epoch 60000 --model_name='CellNEST_lymph_deconvolution' --run_id=1 > output_lymph_deconvolution_run1.log &
nohup cellnest run  --data_name='lymph_deconvolution'  --num_epoch 60000 --model_name='CellNEST_lymph_deconvolution' --run_id=2 > output_lymph_deconvolution_run2.log &
nohup cellnest run  --data_name='lymph_deconvolution'  --num_epoch 60000 --model_name='CellNEST_lymph_deconvolution' --run_id=3 > output_lymph_deconvolution_run3.log &
nohup cellnest run  --data_name='lymph_deconvolution'  --num_epoch 60000 --model_name='CellNEST_lymph_deconvolution' --run_id=4 > output_lymph_deconvolution_run4.log &
nohup cellnest run  --data_name='lymph_deconvolution'  --num_epoch 60000 --model_name='CellNEST_lymph_deconvolution' --run_id=5 > output_lymph_deconvolution_run5.log &
```

Postprocess:
```
cellnest postprocess --data_name='lymph_deconvolution' --model_name='CellNEST_lymph_deconvolution' --total_runs=5 
```

Visualize:
```
cellnest visualize --data_name='lymph_deconvolution' --model_name='CellNEST_lymph_deconvolution'
```

You will find the results under: 'output/lymph_deconvolution/'