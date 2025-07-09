  GNU nano 5.6.1                                                                                                                to_anndata.py                                                                                                                          
# Written By 
# Fatema Tuz Zohora


print('package loading')
import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix
import pandas as pd
import gzip
import argparse
import os
import scanpy as sc
import anndata

print('user input reading')
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ################## Mandatory ####################################################################
    parser.add_argument( '--data_from', type=str, help='Name of the file with path', required=True)  
    parser.add_argument( '--data_to_path', type=str, help='Path to save the dataset in anndata format', required=True)
    parser.add_argument( '--file_name', type=str, help='File name to save the dataset in anndata format', required=True)
    parser.add_argument( '--tissue_position_file', type=str, default='None', help='If your data in --data_from \
    does not have coordinates, then please provide the path to tissue position file.')
   
    args = parser.parse_args()

    temp = sc.read_10x_mtx(args.data_from)  # NOTE: this may change based on your provided format
    count_matrix = sparse.csr_matrix(temp.X)
    print('file read done')
    cell_barcode = np.array(temp.obs.index)
    print('Number of barcodes: %d'%cell_barcode.shape[0])
    gene_ids = list(temp.var_names) 
    # now read the tissue position file. It has the format: barcode vs coordinates    
    df = pd.read_csv(args.tissue_position_file, sep=",", header=None)   
    tissue_position = df.values
    barcode_vs_xy = dict() # record the x and y coordinates for each spot/cell
    for i in range (0, tissue_position.shape[0]):
        barcode_vs_xy[tissue_position[i][0]] = [tissue_position[i][4], tissue_position[i][5]] # x and y coordinates
    
    coordinates = np.zeros((cell_barcode.shape[0], 2)) # insert the coordinates in the order of cell_barcodes
    for i in range (0, cell_barcode.shape[0]):
        coordinates[i,0] = barcode_vs_xy[cell_barcode[i]][0]
        coordinates[i,1] = barcode_vs_xy[cell_barcode[i]][1]

    # now create anndata object
    adata = anndata.AnnData(count_matrix)
    adata.obs_names = cell_barcode
    adata.var_names = gene_ids
    adata.obsm['spatial'] = coordinates
    # save it
    if 'h5ad' not in args.file_name:
        args.file_name = args.file_name + '.h5ad'

    adata.write(args.data_to_path + '/' + args.file_name, compression="gzip")
    print('file write done')  

