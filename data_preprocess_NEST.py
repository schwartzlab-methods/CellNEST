# Written By 
# Fatema Tuz Zohora


print('package loading')
import numpy as np
import pickle
from scipy import sparse
import numpy as np
import qnorm
from scipy.sparse import csr_matrix
from collections import defaultdict
import pandas as pd
import gzip
import argparse
import os
import scanpy as sc

print('user input reading')
#current_dir = 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ################## Mandatory ####################################################################
    parser.add_argument( '--data_name', type=str, help='Name of the dataset', required=True)  
    parser.add_argument( '--data_from', type=str, required=True, help='Path to the dataset to read from. Space Ranger outs/ folder is preferred. Otherwise, provide the *.mtx file of the gene expression matrix.')
    ################# default is set ################################################################
    parser.add_argument( '--data_to', type=str, default='input_graph/', help='Path to save the input graph (to be passed to GAT)')
    parser.add_argument( '--metadata_to', type=str, default='metadata/', help='Path to save the metadata')
    parser.add_argument( '--filter_min_cell', type=int, default=1 , help='Minimum number of cells for gene filtering') 
    parser.add_argument( '--threshold_gene_exp', type=float, default=98, help='Threshold percentile for gene expression. Genes above this percentile are considered active.')
    parser.add_argument( '--tissue_position_file', type=str, default='None', help='If your --data_from argument points to a *.mtx file instead of Space Ranger, then please provide the path to tissue position file.')
    parser.add_argument( '--spot_diameter', type=float, default=89.43, help='Spot/cell diameter for filtering ligand-receptor pairs based on cell-cell contact information. Should be provided in the same unit as spatia data (for Visium, that is pixel).')
    parser.add_argument( '--split', type=int, default=0 , help='How many split sections?') 
    parser.add_argument( '--neighborhood_threshold', type=float, default=0 , help='Set neighborhood threshold distance in terms of same unit as spot diameter') 
    parser.add_argument( '--database_path', type=str, default='database/NEST_database.csv' , help='Provide your desired ligand-receptor database path here. Default database is a combination of CellChat and NicheNet database.') 
    args = parser.parse_args()
    
    if args.neighborhood_threshold == 0:
        args.neighborhood_threshold = args.spot_diameter*4

    if args.data_to == 'input_graph/':
        args.data_to = args.data_to + args.data_name + '/'
    if not os.path.exists(args.data_to):
        os.makedirs(args.data_to)

    if args.metadata_to == 'metadata/':
        args.metadata_to = args.metadata_to + args.data_name + '/'
    if not os.path.exists(args.metadata_to):
        os.makedirs(args.metadata_to)
        
    ####### get the gene id, cell barcode, cell coordinates ######
    print('input data reading')
    if args.tissue_position_file == 'None': # Data is available in Space Ranger output format
        adata_h5 = sc.read_visium(path=args.data_from, count_file='filtered_feature_bc_matrix.h5')
        print('input data read done')
        gene_count_before = len(list(adata_h5.var_names) )    
        sc.pp.filter_genes(adata_h5, min_cells=args.filter_min_cell)
        gene_count_after = len(list(adata_h5.var_names) )  
        print('Gene filtering done. Number of genes reduced from %d to %d'%(gene_count_before, gene_count_after))
        gene_ids = list(adata_h5.var_names)
        coordinates = adata_h5.obsm['spatial']
        cell_barcode = np.array(adata_h5.obs.index)
        print('Number of barcodes: %d'%cell_barcode.shape[0])
        print('Applying quantile normalization')
        temp = qnorm.quantile_normalize(np.transpose(sparse.csr_matrix.toarray(adata_h5.X)))  #https://en.wikipedia.org/wiki/Quantile_normalization
        cell_vs_gene = np.transpose(temp)      
    
    else: # Data is not available in Space Ranger output format
        # read the mtx file
        temp = sc.read_10x_mtx(args.data_from)
        print('*.mtx file read done')
        gene_count_before = len(list(temp.var_names) )
        sc.pp.filter_genes(temp, min_cells=args.filter_min_cell)
        gene_count_after = len(list(temp.var_names) )
        print('Gene filtering done. Number of genes reduced from %d to %d'%(gene_count_before, gene_count_after))
        gene_ids = list(temp.var_names) 
        cell_barcode = np.array(temp.obs.index)
        print('Number of barcodes: %d'%cell_barcode.shape[0])
        print('Applying quantile normalization')
        temp = qnorm.quantile_normalize(np.transpose(sparse.csr_matrix.toarray(temp.X)))  #https://en.wikipedia.org/wiki/Quantile_normalization
        cell_vs_gene = np.transpose(temp)  
    
        
        # now read the tissue position file. It has the format:     
        df = pd.read_csv(args.tissue_position_file, sep=",", header=None)   
        tissue_position = df.values
        barcode_vs_xy = dict() # record the x and y coordinates for each spot/cell
        for i in range (0, tissue_position.shape[0]):
            barcode_vs_xy[tissue_position[i][0]] = [tissue_position[i][4], tissue_position[i][5]] # x and y coordinates
            #barcode_vs_xy[tissue_position[i][0]] = [tissue_position[i][5], tissue_position[i][4]] #for some weird reason, in the .h5 format for LUAD sample, the x and y are swapped
        
        coordinates = np.zeros((cell_barcode.shape[0], 2)) # insert the coordinates in the order of cell_barcodes
        for i in range (0, cell_barcode.shape[0]):
            coordinates[i,0] = barcode_vs_xy[cell_barcode[i]][0]
            coordinates[i,1] = barcode_vs_xy[cell_barcode[i]][1]
        
    
    
    
    ##################### make metadata: barcode_info ###################################
    i=0
    barcode_info=[]
    for cell_code in cell_barcode:
        barcode_info.append([cell_code, coordinates[i,0],coordinates[i,1], 0]) # last entry will hold the component number later
        i=i+1
    ################################################
    
    gene_info=dict()
    for gene in gene_ids:
        gene_info[gene]=''
    
    gene_index=dict()    
    i = 0
    for gene in gene_ids: 
        gene_index[gene] = i
        i = i+1
        
    #### needed if split data is used ##############
    if args.split>0:
        i=0
        node_id_sorted_xy=[]
        for cell_code in cell_barcode:
            node_id_sorted_xy.append([i, coordinates[i,0],coordinates[i,1]])
            i=i+1
        	
        node_id_sorted_xy = sorted(node_id_sorted_xy, key = lambda x: (x[1], x[2]))
        with gzip.open(metadata_to + args.data_name+'_'+'node_id_sorted_xy', 'wb') as fp:  #b, a:[0:5]   
        	pickle.dump(node_id_sorted_xy, fp)
    
    
    ####################################################################
    # ligand - receptor database 
    print('ligand-receptor database reading.')
    df = pd.read_csv(args.database_path, sep=",")
    
    '''
            Ligand   Receptor          Annotation           Reference
    0        TGFB1     TGFBR1  Secreted Signaling      KEGG: hsa04350
    1        TGFB1     TGFBR2  Secreted Signaling      KEGG: hsa04350
    '''
    print('ligand-receptor database reading done.')
    print('Preprocess start.')
    ligand_dict_dataset = defaultdict(list)
    cell_cell_contact = dict() 
    count_pair = 0
    for i in range (0, df["Ligand"].shape[0]):
        ligand = df["Ligand"][i]
        if ligand not in gene_info: # not found in the dataset
            continue    
            
        receptor = df["Receptor"][i]
        if receptor not in gene_info: # not found in the dataset
            continue   
            
        ligand_dict_dataset[ligand].append(receptor)
        gene_info[ligand] = 'included'
        gene_info[receptor] = 'included'
        count_pair = count_pair + 1
        
        if df["Annotation"][i] == 'Cell-Cell Contact':
            cell_cell_contact[receptor] = '' # keep track of which ccc are labeled as cell-cell-contact
    
    
    print('number of ligand-receptor pairs in this dataset %d '%count_pair) 
    print('number of ligands %d '%len(ligand_dict_dataset.keys()))
    
    included_gene=[]
    for gene in gene_info.keys(): 
        if gene_info[gene] == 'included':
            included_gene.append(gene)
            
    print('Total genes in this dataset: %d, number of genes working as ligand and/or receptor: %d '%(len(gene_ids),len(included_gene)))
    
    # assign id to each entry in the ligand-receptor database
    l_r_pair = dict()
    lr_id = 0
    for gene in list(ligand_dict_dataset.keys()): 
        ligand_dict_dataset[gene]=list(set(ligand_dict_dataset[gene]))
        l_r_pair[gene] = dict()
        for receptor_gene in ligand_dict_dataset[gene]:
            l_r_pair[gene][receptor_gene] = lr_id 
            lr_id  = lr_id  + 1
        
    
    ###################################################################################
    # build physical distance matrix
    from sklearn.metrics.pairwise import euclidean_distances
    distance_matrix = euclidean_distances(coordinates, coordinates)
    
    # assign weight to the neighborhood relations based on neighborhood distance 
    dist_X = np.zeros((distance_matrix.shape[0], distance_matrix.shape[1]))
    for j in range(0, distance_matrix.shape[1]): # look at all the incoming edges to node 'j'
        max_value=np.max(distance_matrix[:,j]) # max distance of node 'j' to all it's neighbors (incoming)
        min_value=np.min(distance_matrix[:,j]) # min distance of node 'j' to all it's neighbors (incoming)
        for i in range(distance_matrix.shape[0]):
            dist_X[i,j] = 1-(distance_matrix[i,j]-min_value)/(max_value-min_value) # scale the distance of node 'j' to all it's neighbors (incoming) and flip it so that nearest one will have maximum weight.
            	
        #list_indx = list(np.argsort(dist_X[:,j]))
        #k_higher = list_indx[len(list_indx)-k_nn:len(list_indx)]
        for i in range(0, distance_matrix.shape[0]):
            if distance_matrix[i,j] > args.neighborhood_threshold: #i not in k_higher:
                dist_X[i,j] = 0 # no ccc happening outside threshold distance
                
    #cell_rec_count = np.zeros((cell_vs_gene.shape[0]))
    #####################################################################################
    # Set threshold gene percentile
    cell_percentile = []
    for i in range (0, cell_vs_gene.shape[0]):
        y = sorted(cell_vs_gene[i]) # sort each row/cell in ascending order of gene expressions
        cell_percentile.append(np.percentile(y, args.threshold_gene_exp)) 
    
    ##############################################################################
    # some preprocessing before making the input graph
    count_total_edges = 0
    
    cells_ligand_vs_receptor = []
    for i in range (0, cell_vs_gene.shape[0]):
        cells_ligand_vs_receptor.append([])
        
    for i in range (0, cell_vs_gene.shape[0]):
        for j in range (0, cell_vs_gene.shape[0]):
            cells_ligand_vs_receptor[i].append([])
            cells_ligand_vs_receptor[i][j] = []

    ligand_list =  list(ligand_dict_dataset.keys())            
    start_index = 0 #args.slice
    end_index = len(ligand_list) #min(len(ligand_list), start_index+100)
    
    for g in range(start_index, end_index): 
        gene = ligand_list[g]
        for i in range (0, cell_vs_gene.shape[0]): # ligand
              
            if cell_vs_gene[i][gene_index[gene]] < cell_percentile[i]:
                continue
            
            for j in range (0, cell_vs_gene.shape[0]): # receptor
                if dist_X[i,j]==0: #distance_matrix[i,j] >= args.neighborhood_threshold: #spot_diameter*4
                    continue
    
                for gene_rec in ligand_dict_dataset[gene]:
                    if cell_vs_gene[j][gene_index[gene_rec]] >= cell_percentile[j]: # or cell_vs_gene[i][gene_index[gene]] >= cell_percentile[i][4] :#gene_list_percentile[gene_rec][1]: #global_percentile: #
                        if gene_rec in cell_cell_contact and distance_matrix[i,j] > args.spot_diameter:
                            continue
    
                        communication_score = cell_vs_gene[i][gene_index[gene]] * cell_vs_gene[j][gene_index[gene_rec]]
                        relation_id = l_r_pair[gene][gene_rec]
    
                        if communication_score<=0:
                            print('zero valued ccc score found. Might be a potential ERROR!! ')
                            continue	
                            
                        cells_ligand_vs_receptor[i][j].append([gene, gene_rec, communication_score, relation_id])
                        count_total_edges = count_total_edges + 1
                        
        print('%d genes done out of %d ligand genes'%(g+1, len(ligand_list)))
    
    
    #print('total number of edges in the input graph %d '%count_total_edges)
    ################################################################################
    # input graph generation
    ccc_index_dict = dict()
    row_col = [] # list of input edges, row = from node, col = to node
    edge_weight = [] # 3D edge features in the same order as row_col
    lig_rec = [] # ligand and receptors corresponding to the edges in the same order as row_col
    self_loop_found = defaultdict(dict) # to keep track of self-loops -- used later during visualization plotting
    for i in range (0, len(cells_ligand_vs_receptor)):
        #ccc_j = []
        for j in range (0, len(cells_ligand_vs_receptor)):
            if dist_X[i,j]>0: #distance_matrix[i][j] <= args.neighborhood_threshold: 
                count_local = 0
                if len(cells_ligand_vs_receptor[i][j])>0:
                    for k in range (0, len(cells_ligand_vs_receptor[i][j])):
                        gene = cells_ligand_vs_receptor[i][j][k][0]
                        gene_rec = cells_ligand_vs_receptor[i][j][k][1]
                        ligand_receptor_coexpression_score = cells_ligand_vs_receptor[i][j][k][2]
                        row_col.append([i,j])
                        edge_weight.append([dist_X[i,j], ligand_receptor_coexpression_score, cells_ligand_vs_receptor[i][j][k][3]])
                        lig_rec.append([gene, gene_rec])
                                                  
                        if i==j: # self-loop
                            self_loop_found[i][j] = ''
    

    total_num_cell = cell_vs_gene.shape[0]
    print('total number of nodes is %d, and edges is %d in the input graph'%(total_num_cell, len(row_col)))
    print('preprocess done.')
    print('writing data ...')

    ################## input graph #################################################
    with gzip.open(args.data_to + args.data_name + '_adjacency_records', 'wb') as fp:  
        pickle.dump([row_col, edge_weight, lig_rec, total_num_cell], fp)

    ################# metadata #####################################################
    with gzip.open(args.metadata_to + args.data_name +'_self_loop_record', 'wb') as fp: 
        pickle.dump(self_loop_found, fp)

    with gzip.open(args.metadata_to + args.data_name +'_barcode_info', 'wb') as fp:  
        pickle.dump(barcode_info, fp)
    
    ################## required for the nest interactive version ###################
    df = pd.DataFrame(gene_ids)
    df.to_csv(args.metadata_to + 'gene_ids_'+args.data_name+'.csv', index=False, header=False)
    df = pd.DataFrame(cell_barcode)
    df.to_csv(args.metadata_to + 'cell_barcode_'+args.data_name+'.csv', index=False, header=False)
    df = pd.DataFrame(coordinates)
    df.to_csv(args.metadata_to + 'coordinates_'+args.data_name+'.csv', index=False, header=False)
    
    
    ######### optional #############################################################           
    # we do not need this to use anywhere. But just for debug purpose we are saving this. We can skip this if we have space issue.
    with gzip.open(args.data_to + args.data_name + '_cell_vs_gene_quantile_transformed', 'wb') as fp:  
    	pickle.dump(cell_vs_gene, fp)
        
    print('write data done')
    
# nohup python -u data_preprocess_NEST.py --data_name='PDAC_64630_mincell3_th98p5' --data_from='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/' --filter_min_cell=3 --threshold_gene_exp=98.5 > output_data_preprocess_PDAC_64630_min_cell_3_th98p5.log &
