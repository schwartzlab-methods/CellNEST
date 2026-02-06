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
import json
import gc
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics.pairwise import euclidean_distances
print('user input reading')
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ################## Mandatory ####################################################################
    parser.add_argument( '--data_name', type=str, help='Name of the dataset', required=True)  
    parser.add_argument( '--data_from', type=str, required=True, help='Path to the dataset to read from. Space Ranger outs/ folder is preferred. Otherwise, provide the *.mtx file of the gene expression matrix.')
    ################# default is set ################################################################
    parser.add_argument( '--data_type', type=str, default='anndata', help='Set one of these two types [visium, anndata]. visium refers to spot based Visium file.')
    parser.add_argument( '--data_to', type=str, default='input_graph/', help='Path to save the input graph (to be passed to GAT)')
    parser.add_argument( '--metadata_to', type=str, default='metadata/', help='Path to save the metadata')
    parser.add_argument( '--filter_min_cell', type=int, default=1 , help='Minimum number of cells for gene filtering') 
    parser.add_argument( '--threshold_gene_exp', type=float, default=98, help='Threshold percentile for gene expression. Genes above this percentile are considered active.')
    parser.add_argument( '--tissue_position_file', type=str, default='None', help='If your --data_from argument points to a *.mtx file instead of Space Ranger, then please provide the path to tissue position file.')
    parser.add_argument( '--juxtacrine_distance', type=float, default=-1, help='Distance for filtering ligand-receptor pairs based on cell-cell contact information. Automatically calculated unless provided. It has the same unit as the coordinates (for Visium, that is pixel).')
    parser.add_argument( '--split', type=int, default=0 , help='Set to 1 if you plan to run training by splitting the input graph into multiple subgraphs') 
    parser.add_argument( '--distance_measure', type=str, default='fixed' , help='Set neighborhood cutoff criteria. Choose from [knn, fixed]')
    parser.add_argument( '--k', type=int, default=50 , help='Set neighborhood cutoff number. This will be used if --distance_measure=knn')    
    parser.add_argument( '--neighborhood_threshold', type=float, default=0 , help='Set neighborhood threshold distance in terms of same unit as the coordinates') 
    parser.add_argument( '--block_autocrine', type=int, default=0 , help='Set to 1 if you want to ignore autocrine signals.') 
    parser.add_argument( '--block_juxtacrine', type=int, default=0, help='Set to 1 if you want to ignore juxtacrine signals.')     
    parser.add_argument( '--database_path', type=str, default='database/CellNEST_database.csv' , help='Provide your desired ligand-receptor database path here. Default database is a combination of CellChat and NicheNet database.') 
    parser.add_argument( '--set_ROI', type=int, default=0, help='Set to 1 if you want to use ROI')
    parser.add_argument( '--x_min', type=int, default=-1, help='Set if you want to use ROI')
    parser.add_argument( '--x_max', type=int, default=-1, help='Set if you want to use ROI')
    parser.add_argument( '--y_min', type=int, default=-1, help='Set if you want to use ROI')
    parser.add_argument( '--y_max', type=int, default=-1, help='Set if you want to use ROI')
    args = parser.parse_args()
    
    if args.data_to == 'input_graph/':
        args.data_to = args.data_to + args.data_name + '/'
    if not os.path.exists(args.data_to):
        os.makedirs(args.data_to)

    if args.metadata_to == 'metadata/':
        args.metadata_to = args.metadata_to + args.data_name + '/'
    if not os.path.exists(args.metadata_to):
        os.makedirs(args.metadata_to)
        
    ####### get the gene id, cell barcode, cell coordinates ######
    print('Input data reading')
    if args.tissue_position_file == 'None': # Data is available in Space Ranger output format
        if args.data_type == 'visium':
            adata_h5 = sc.read_visium(path=args.data_from, count_file='filtered_feature_bc_matrix.h5')
            print('input data read done')
            gene_count_before = len(list(adata_h5.var_names) )    
            sc.pp.filter_genes(adata_h5, min_cells=args.filter_min_cell)
            gene_count_after = len(list(adata_h5.var_names) )  
            print('Gene filtering done. Number of genes reduced from %d to %d'%(gene_count_before, gene_count_after))
            gene_ids = list(adata_h5.var_names)
            gene_ids = [gene_id.upper() for gene_id in gene_ids]
            coordinates = adata_h5.obsm['spatial']
            cell_barcode = np.array(adata_h5.obs.index)
            print('Number of barcodes: %d'%cell_barcode.shape[0])
            print('Applying quantile normalization')
            temp = qnorm.quantile_normalize(np.transpose(sparse.csr_matrix.toarray(adata_h5.X)))  #https://en.wikipedia.org/wiki/Quantile_normalization
            cell_vs_gene = np.transpose(temp)
            ####### use the spot diameter as the --juxtacrine_distance
            if args.juxtacrine_distance == -1:
                file = open(args.data_from+'/spatial/scalefactors_json.json', 'r')
                data = json.load(file)
                spot_diameter = data["spot_diameter_fullres"]
                args.juxtacrine_distance = spot_diameter

        elif args.data_type == 'anndata':
            adata_h5 = sc.read_h5ad(args.data_from)
            print('input data read done')
            gene_count_before = len(list(adata_h5.var_names))    
            sc.pp.filter_genes(adata_h5, min_cells=args.filter_min_cell)
            gene_count_after = len(list(adata_h5.var_names) )  
            print('Gene filtering done. Number of genes reduced from %d to %d'%(gene_count_before, gene_count_after))
            gene_ids = list(adata_h5.var_names)
            gene_ids = [gene_id.upper() for gene_id in gene_ids]
            coordinates = np.array(adata_h5.obsm['spatial'])
            if args.set_ROI == 1:
                if args.x_min == -1:
                    args.x_min = np.min(coordinates[:][0])
                if args.x_max == -1:
                    args.x_max = np.max(coordinates[:][0])
                if args.y_min == -1:
                    args.y_min = np.min(coordinates[:][1])
                if args.y_max == -1:
                    args.y_max = np.max(coordinates[:][1])

                keep_cells = []
                for i in range(0, coordinates.shape[0]):
                    if args.x_min <= coordinates[i][0] and coordinates[i][0] <= args.x_max:
                        if args.y_min <= coordinates[i][1] and coordinates[i][1] <= args.y_max:
                            keep_cells.append(i)

                adata_h5 = adata_h5[keep_cells]
                print('after ROI cropping:')

            print(adata_h5)




            cell_barcode = np.array(adata_h5.obs_names) #obs.index)
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
        gene_ids = [gene_id.upper() for gene_id in gene_ids]
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
        
    
    #print(euclidean_distances(coordinates[0:1],coordinates[2:3]))
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
        


    # build physical distance matrix
    print('Build physical distance matrix')
    if args.distance_measure == 'fixed':
        # then you need to set the args.neighborhood_threshold
        if args.neighborhood_threshold == 0:
            distance_matrix = euclidean_distances(coordinates, coordinates)
            #### then automatically calculate the min distance between two nodes #########
            sorted_first_row = np.sort(distance_matrix[0,:])
            distance_a_b = sorted_first_row[1]
            args.neighborhood_threshold = distance_a_b * 4 # 4 times the distance between two nodes
            distance_matrix = 0
            gc.collect()

        nbrs = NearestNeighbors(radius=args.neighborhood_threshold, algorithm='kd_tree', n_jobs=-1)
        nbrs.fit(coordinates)
        distances, indices = nbrs.kneighbors(coordinates)
        print('Neighborhood distance is set to be %g (same unit as the coordinates)'%(args.neighborhood_threshold))
    else:
        print('Neighborhood distance is set to be %d nearest neighbors'%args.k)
        nbrs = NearestNeighbors(n_neighbors=args.k, algorithm='kd_tree', n_jobs=-1)
        nbrs.fit(coordinates)
        distances, indices = nbrs.kneighbors(coordinates)

    unique_distances = np.unique(distances)
    distance_a_b = sorted(unique_distances)[1]
    # distances: array of shape (n_cells, k) with the Euclidean distance from 
    # each cell to its k neighbors.
    # indices: array of shape (n_cells, k) with the neighbor indices (row indices of X).
    # each cell 'j' will receive signal from its neighbors 'i' in descending order of weights
    print('Assign weight to the neighborhood relations based on neighborhood distance')
    weightdict_i_to_j = defaultdict(dict)
    for cell_idx in range (0, indices.shape[0]):
        max_value = np.max(distances[cell_idx,:])
        min_value = np.min(distances[cell_idx,:])
        for neigh_idx in range (0, indices.shape[1]):
            neigh_cell_idx = indices[cell_idx][neigh_idx]
            distance_neigh_cell = distances[cell_idx][neigh_idx]
            flipped_distance_neigh_cell = 1-(distance_neigh_cell-min_value)/(max_value-min_value)
            # i = neigh_cell_idx, j = cell_idx
            weightdict_i_to_j[neigh_cell_idx][cell_idx] = flipped_distance_neigh_cell


    # assign weight to the neighborhood relations based on neighborhood distance 
    """
    dist_X = np.zeros((distance_matrix.shape[0], distance_matrix.shape[1]))
    for j in range(0, distance_matrix.shape[1]): # look at all the incoming edges to node 'j'
        max_value=np.max(distance_matrix[:,j]) # max distance of node 'j' to all it's neighbors (incoming)
        min_value=np.min(distance_matrix[:,j]) # min distance of node 'j' to all it's neighbors (incoming)
        for i in range(distance_matrix.shape[0]):
            dist_X[i,j] = 1-(distance_matrix[i,j]-min_value)/(max_value-min_value) # scale the distance of node 'j' to all it's neighbors (incoming) and flip it so that nearest one will have maximum weight.
            	
        if args.distance_measure=='knn':
            list_indx = list(np.argsort(dist_X[:,j]))
            k_higher = list_indx[len(list_indx)-args.k:len(list_indx)]
            for i in range(0, distance_matrix.shape[0]):
                if i not in k_higher:
                    dist_X[i,j] = 0 #-1
        else:
            for i in range(0, distance_matrix.shape[0]):
                if distance_matrix[i,j] > args.neighborhood_threshold: #i not in k_higher:
                    dist_X[i,j] = 0 # no ccc happening outside threshold distance
    



        #list_indx = list(np.argsort(dist_X[:,j]))
        #k_higher = list_indx[len(list_indx)-k_nn:len(list_indx)]
        for i in range(0, distance_matrix.shape[0]):
            if distance_matrix[i,j] > args.neighborhood_threshold: #i not in k_higher:
                dist_X[i,j] = 0 # no ccc happening outside threshold distance
    """      
    #cell_rec_count = np.zeros((cell_vs_gene.shape[0]))


    #### automatically set the args.juxtacrine_distance ############
    if args.juxtacrine_distance == -1: 
        args.juxtacrine_distance = distance_a_b

    if args.block_juxtacrine == 0:
        print("Auto calculated juxtacrine distance is %g. To change it use --juxtacrine_distance"%args.juxtacrine_distance)
    #### sort the nodes in the ascending order of x and y coordinates ####################
    i=0
    node_id_sorted_xy=[]
    for cell_code in cell_barcode:
        node_id_sorted_xy.append([i, coordinates[i,0],coordinates[i,1]])
        i=i+1

    node_id_sorted_xy = sorted(node_id_sorted_xy, key = lambda x: (x[1], x[2]))
    
    #### needed if split data is used ##############
    if args.split>0:
        with gzip.open(args.metadata_to + args.data_name+'_'+'node_id_sorted_xy', 'wb') as fp:  #b, a:[0:5]   
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
        
    print('number of ligand-receptor pairs in this dataset %d '%lr_id)     
    ###################################################################################

    #####################################################################################
    # Set threshold gene percentile
    cell_percentile = []
    for i in range (0, cell_vs_gene.shape[0]):
        y = sorted(cell_vs_gene[i]) # sort each row/cell in ascending order of gene expressions
        ## inter ##
        active_cutoff = np.percentile(y, args.threshold_gene_exp)
        if active_cutoff == min(cell_vs_gene[i][:]):
            times = 1
            while active_cutoff == min(cell_vs_gene[i][:]):
                new_threshold = args.threshold_gene_exp + 5 * times                    
                if new_threshold >= 100:
                    active_cutoff = max(cell_vs_gene[i][:])  
                    if active_cutoff == min(cell_vs_gene[i][:]): # still same as the min --> flat curve: something wrong
                        active_cutoff = max(cell_vs_gene[i][:]) + 1 # This cell will be skipped

                    break
                active_cutoff = np.percentile(y, new_threshold)
                times = times + 1 

        cell_percentile.append(active_cutoff)     

    print('set threshold gene percentile done')
    ##############################################################################
    # some preprocessing before making the input graph
    count_total_edges = 0
    """
    cells_ligand_vs_receptor = []
    for i in range (0, cell_vs_gene.shape[0]):
        cells_ligand_vs_receptor.append([])
        
    for i in range (0, cell_vs_gene.shape[0]):
        for j in range (0, cell_vs_gene.shape[0]):
            cells_ligand_vs_receptor[i].append([])
            cells_ligand_vs_receptor[i][j] = []
    """
    cells_ligand_vs_receptor = defaultdict(dict)

    ligand_list =  list(ligand_dict_dataset.keys())            
    start_index = 0 #args.slice
    end_index = len(ligand_list) #min(len(ligand_list), start_index+100)
    print('some preprocessing before making the input graph')    
    for g in range(start_index, end_index): 
        gene = ligand_list[g]
        for i in weightdict_i_to_j:
            if cell_vs_gene[i][gene_index[gene]] < cell_percentile[i]:
                continue            
            for j in weightdict_i_to_j[i]: # receptor
                if args.block_autocrine == 1 and i==j:
                    continue
                for gene_rec in ligand_dict_dataset[gene]:
                    if cell_vs_gene[j][gene_index[gene_rec]] >= cell_percentile[j]: # or cell_vs_gene[i][gene_index[gene]] >= cell_percentile[i][4] :#gene_list_percentile[gene_rec][1]: #global_percentile: #
                        if (gene_rec in cell_cell_contact) and (args.block_juxtacrine==1 or euclidean_distances(coordinates[i:i+1],coordinates[j:j+1]) > args.juxtacrine_distance):
                            continue
    
                        communication_score = cell_vs_gene[i][gene_index[gene]] * cell_vs_gene[j][gene_index[gene_rec]]
                        relation_id = l_r_pair[gene][gene_rec]
    
                        if communication_score<=0:
                            print('zero valued ccc score found. Might be a potential ERROR!! ')
                            continue	
                            
                        #cells_ligand_vs_receptor[i][j].append([gene, gene_rec, communication_score, relation_id])
                        if i in cells_ligand_vs_receptor:
                            if j in cells_ligand_vs_receptor[i]:
                                cells_ligand_vs_receptor[i][j].append([gene, gene_rec, communication_score, relation_id])
                            else:
                                cells_ligand_vs_receptor[i][j] = []
                                cells_ligand_vs_receptor[i][j].append([gene, gene_rec, communication_score, relation_id])
                        else:
                            cells_ligand_vs_receptor[i][j] = []
                            cells_ligand_vs_receptor[i][j].append([gene, gene_rec, communication_score, relation_id])

                        count_total_edges = count_total_edges + 1
                        
        print('%d/%d ligand genes processed'%(g+1, len(ligand_list)), end='\r')
    
    print('')    
    #print('total number of edges in the input graph %d '%count_total_edges)
    ################################################################################
    # input graph generation
    ccc_index_dict = dict()
    row_col = [] # list of input edges, row = from node, col = to node
    edge_weight = [] # 3D edge features in the same order as row_col
    lig_rec = [] # ligand and receptors corresponding to the edges in the same order as row_col
    self_loop_found = defaultdict(dict) # to keep track of self-loops -- used later during visualization plotting
    for i in cells_ligand_vs_receptor.keys():
        #ccc_j = []
        for j in cells_ligand_vs_receptor[i].keys():
            if i in weightdict_i_to_j and j in weightdict_i_to_j[i]: #dist_X[i,j]>0: #distance_matrix[i][j] <= args.neighborhood_threshold: 
                count_local = 0
                if len(cells_ligand_vs_receptor[i][j])>0:
                    for k in range (0, len(cells_ligand_vs_receptor[i][j])):
                        gene = cells_ligand_vs_receptor[i][j][k][0]
                        gene_rec = cells_ligand_vs_receptor[i][j][k][1]
                        ligand_receptor_coexpression_score = cells_ligand_vs_receptor[i][j][k][2]
                        row_col.append([i,j])
                        edge_weight.append([weightdict_i_to_j[i][j], ligand_receptor_coexpression_score, cells_ligand_vs_receptor[i][j][k][3]])
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
    
    ################## required for the CellNEST interactive version ###################
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
    