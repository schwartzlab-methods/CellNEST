If you have data in non-Visium format, the easiest way is to convert the data format to anndata. We will use LUAD sample for explaining it. 
Let us assume we have downloaded the LUAD sample into a directory "CellNEST/data/LUAD_GSM5702473_TD1/". 
It comes as a *.mtx format, so following script can be run to convert it to anndata format (assuming that your current working directory is CellNEST):
```
python to_anndata.py --data_from=data/LUAD_GSM5702473_TD1/ --data_to_path=data/LUAD_GSM5702473_TD1/ \
--file_name=LUAD_GSM5702473_TD1 --tissue_position_file=data/LUAD_GSM5702473_TD1/GSM5702473_TD1_tissue_positions_list.csv
```
This is just a sample script and you are welcome to modify to_anndata.py according to your input format. Just make sure that 
the anndata has following rows and columns, as well as, spatial coordinates:
```
    # now create anndata object
    adata = anndata.AnnData(count_matrix)
    adata.obs_names = cell_barcode
    adata.var_names = gene_ids
    adata.obsm['spatial'] = coordinates
```


After that, we will run preprocessing step with **--data_type=anndata** parameter as follows:

```
./cellnest preprocess --data_name=LUAD_TD1 --data_type=anndata --data_from=data/LUAD/LUAD_GSM5702473_TD1/LUAD_GSM5702473_TD1.h5ad

```
