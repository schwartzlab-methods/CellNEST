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
import anndata
#import altair as alt
import datetime

print('user input reading')
#current_dir = 
if __name__ == "__main__":
    start_time = datetime.datetime.now()
    parser = argparse.ArgumentParser()
    ################## Mandatory ####################################################################
    parser.add_argument( '--data_name', type=str, help='Name of the dataset', default="Visium_HD_Human_Colon_Cancer_square_002um_outputs_full")  
    parser.add_argument( '--data_from', type=str, default='data/segmented_visium_hd_CRC/', help='Path to the dataset to read from. Space Ranger outs/ folder is preferred. Otherwise, provide the *.mtx file of the gene expression matrix.')
    ################# default is set ################################################################
    parser.add_argument( '--data_to', type=str, default='input_graph/', help='Path to save the input graph (to be passed to GAT)')
    parser.add_argument( '--metadata_to', type=str, default='metadata/', help='Path to save the metadata')
    parser.add_argument( '--filter_min_cell', type=int, default=1000 , help='Minimum number of cells for gene filtering') 
    parser.add_argument( '--threshold_gene_exp', type=float, default=98.5, help='Threshold percentile for gene expression. Genes above this percentile are considered active.')
    parser.add_argument( '--tissue_position_file', type=str, default='data/Visium_HD_Human_Colon_Cancer_square_002um_outputs/spatial/tissue_positions.parquet', help='If your --data_from argument points to a *.mtx file instead of Space Ranger, then please provide the path to tissue position file.')
    parser.add_argument( '--spot_diameter', type=float, default=37.04, help='Spot/cell diameter for filtering ligand-receptor pairs based on cell-cell contact information. Should be provided in the same unit as spatia data (for Visium, that is pixel).')
    parser.add_argument( '--split', type=int, default=1, help='Will you apply split options. Set 1 if yes.') 
    parser.add_argument( '--distance_measure', type=str, default='knn' , help='Set neighborhood cutoff criteria')
    parser.add_argument( '--k', type=int, default=50 , help='Set neighborhood cutoff number')    
    parser.add_argument( '--neighborhood_threshold', type=float, default=0, help='Set neighborhood threshold distance in terms of same unit as spot diameter') 
    parser.add_argument( '--database_path', type=str, default='database/NEST_database.csv' , help='Provide your desired ligand-receptor database path here. Default database is a combination of CellChat and NicheNet database.') 
    #######################################
    parser.add_argument( '--ROI', type=int, default=0 , help='Do you have a ROI')
    parser.add_argument( '--x_min', type=int, default=0 , help='Set x_min')  
    parser.add_argument( '--x_max', type=int, default=54000 , help='Set x_max')  
    parser.add_argument( '--y_min', type=int, default=0 , help='Set y_min')  
    parser.add_argument( '--y_max', type=int, default=24000 , help='Set y_max')  
    
    
    args = parser.parse_args()
    k_nn = args.k
    
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

    adata_h5 = anndata.read_h5ad(args.data_from + 'count_area_filtered_adata_p75.h5ad')
    print('input data read done')
    gene_count_before = len(list(adata_h5.var_names) )    
    sc.pp.filter_genes(adata_h5, min_cells=args.filter_min_cell)
    gene_count_after = len(list(adata_h5.var_names) )  
    print('Gene filtering done. Number of genes reduced from %d to %d'%(gene_count_before, gene_count_after))
    gene_ids = list(adata_h5.var_names)
    cell_id = np.array(adata_h5.obs.index) #id
    print('Number of barcodes: %d'%cell_id.shape[0])
    print('Applying quantile normalization')
    temp = qnorm.quantile_normalize(np.transpose(sparse.csr_matrix.toarray(adata_h5.X)))  #https://en.wikipedia.org/wiki/Quantile_normalization
    cell_vs_gene = np.transpose(temp)      

    fp = gzip.open(args.data_from + args.data_name + '_coordinate_barcode', 'rb')
    coordinates, cell_barcode = pickle.load(fp)
     
    
    ##################### make metadata: barcode_info ###################################
    i=0
    barcode_info=[]
    for cell_code in cell_id:
        barcode_info.append([cell_code, coordinates[i,0],coordinates[i,1], 0]) # last entry will hold the component number later
        i=i+1

    '''
    ###################### filter it to keep only the cells that are inside the region of interest ##################
    # filter barcode info
    x_max = 54000
    #y_max = 234225
    temp_barcode_info = []
    cell_code_active = dict()
    cell_index_active = []
    
    for i in range (0, len(barcode_info)):
        if barcode_info[i][1] <= 54000 and  barcode_info[i][2] <= y_max:
            # keep it
            temp_barcode_info.append(barcode_info[i])
            cell_code_active[barcode_info[i][0]] = ''
            cell_index_active.append(i)
            
    barcode_info = temp_barcode_info
    
    # filter the coordinate
    coordinates = np.zeros((len(barcode_info), 2)) 
    for i in range (0, len(barcode_info)):
        coordinates[i,0] = barcode_info[i][1]
        coordinates[i,1] = barcode_info[i][2]
                
    
    # filter the cell_id
    temp_cell_id = []
    for cell_code in cell_id:
        if cell_code in cell_code_active:
            temp_cell_id.append(cell_code)
            
    cell_id = temp_cell_id
    
    # filter cell_vs_gene as well
    cell_vs_gene = cell_vs_gene[cell_index_active, :]
    
    # Now print to see if the dimensions are good to proceed. 
    '''
    
    ############################ Now plot it to see how does it look ###################
    '''
    data_list=dict()
    #data_list['pathology_label']=[]
    data_list['X']=[]
    data_list['Y']=[]     

    for i in range (0, len(barcode_info)):        
        data_list['X'].append(barcode_info[i][1])
        data_list['Y'].append(barcode_info[i][2])

   
    data_list_pd = pd.DataFrame(data_list)
    chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
        alt.X('X', scale=alt.Scale(zero=False)),
        alt.Y('Y', scale=alt.Scale(zero=False)),
    )
    chart.save('/cluster/home/t116508uhn/' + args.data_name +'_tissue_altair_plot.html')
    print('Altair plot generation done')    
    '''
    ################################################
    
    gene_info=dict()
    for gene in gene_ids:
        gene_info[gene]=''
    
    gene_index=dict()    
    i = 0
    for gene in gene_ids: 
        gene_index[gene] = i
        i = i+1
        
   
    ###################################### Neighborhood Cutoff ###########################################
    # build physical distance matrix
    from sklearn.metrics.pairwise import euclidean_distances
    distance_matrix = np.zeros((len(cell_id), len(cell_id)))
    print('calculating Euclidean distance')
    for partition_id in range (1, 6+1):
        distance_matrix[:,(len(cell_id)*(partition_id-1))//6 :(len(cell_id)*partition_id)//6] = euclidean_distances(coordinates, coordinates[(len(cell_id)*(partition_id-1))//6 :(len(cell_id)*partition_id)//6])
        print("%d to %d"%((len(cell_id)*(partition_id-1))//6, (len(cell_id)*partition_id)//6))
    
    
    # assign weight to the neighborhood relations based on neighborhood distance ################
    print('Assign weight to the neighborhood relations based on neighborhood distance. This is a time consuming step.')
    
    dist_X = np.zeros((distance_matrix.shape[0], distance_matrix.shape[1]))
    for j in range(0, distance_matrix.shape[1]): # look at all the incoming edges to node 'j'
        max_value=np.max(distance_matrix[:,j]) # max distance of node 'j' to all it's neighbors (incoming)
        min_value=np.min(distance_matrix[:,j]) # min distance of node 'j' to all it's neighbors (incoming)
        for i in range(distance_matrix.shape[0]):
            dist_X[i,j] = 1-(distance_matrix[i,j]-min_value)/(max_value-min_value) # scale the distance of node 'j' to all it's neighbors (incoming) and flip it so that nearest one will have maximum weight.
            	
        if args.distance_measure=='knn':
            list_indx = list(np.argsort(dist_X[:,j]))
            k_higher = list_indx[len(list_indx)-k_nn:len(list_indx)]
            for i in range(0, distance_matrix.shape[0]):
                if i not in k_higher:
                    dist_X[i,j] = 0 #-1
        else:
            for i in range(0, distance_matrix.shape[0]):
                if distance_matrix[i,j] > args.neighborhood_threshold: #i not in k_higher:
                    dist_X[i,j] = 0 # no ccc happening outside threshold distance
    

    dist_X_list = []
    for i in range (0, dist_X.shape[0]):
        for j in range (0, dist_X.shape[1]):
            if dist_X[i,j]!=0:
                dist_X_list.append([i, j, dist_X[i,j]])
                
    print('dist_X_list length is %d'%len(dist_X_list)) #len is 595,460
    #with gzip.open(args.metadata_to + args.data_name + '_distance_weight', 'wb') as fp:  
    #    pickle.dump(dist_X_list, fp)
    
    #############################################################################################
    #with gzip.open('metadata/Visium_HD_Human_Colon_Cancer_square_002um_outputs/Visium_HD_Human_Colon_Cancer_square_002um_outputs_distance_weight', 'rb') as fp:  
    #    dist_X_list = pickle.load(fp)    

    dist_X_dict = defaultdict(dict)
    for k in range (0, len(dist_X_list)):
        i = dist_X_list[k][0]
        j = dist_X_list[k][1]
        score = dist_X_list[k][2]
        if i not in dist_X_dict:
            dist_X_dict[i][j] = score
        elif j not in dist_X_dict[i]:
            dist_X_dict[i][j] = score
        else:
            print('error')

    
    #cell_rec_count = np.zeros((cell_vs_gene.shape[0]))
    #####################################################################################        

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
    ligand_gene = dict()
    receptor_gene = dict()
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

        ligand_gene[ligand] = ''
        receptor_gene[receptor] = ''
        
        count_pair = count_pair + 1
        
        if df["Annotation"][i] == 'Cell-Cell Contact':
            cell_cell_contact[receptor] = '' # keep track of which ccc are labeled as cell-cell-contact
    
    
    print('number of ligand-receptor pairs in this dataset %d '%count_pair) 
    print('number of ligands %d '%len(ligand_dict_dataset.keys()))
    print('number of ligands %d '%len(ligand_gene.keys()))
    print('number of receptor %d '%len(receptor_gene.keys()))
    
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
    # Set threshold gene percentile
    cell_percentile = []
    for i in range (0, cell_vs_gene.shape[0]):
        y = sorted(cell_vs_gene[i]) # sort each row/cell in ascending order of gene expressions
        ## inter ##
        active_cutoff = np.percentile(y, args.threshold_gene_exp)
        if active_cutoff == min(cell_vs_gene[i][:]):
            active_cutoff = max(cell_vs_gene[i][:])
            #all_deactive_count = all_deactive_count + 1
        cell_percentile.append(active_cutoff) 

   
    ##############################################################################
    active_nodes = dict()
    if args.ROI == 1:
        print('Picking the active nodes inside region of interest (ROI)')
        for i in range (0, len(barcode_info)):
            if (args.x_min <= barcode_info[i][1] and barcode_info[i][1] <= args.x_max) and  (args.y_min <= barcode_info[i][2] and barcode_info[i][2] <= args.y_max):
                active_nodes[i] = ''        

    ############################ Now plot it to see how does it look ###################
        
        data_list=dict()
        data_list['X']=[]
        data_list['Y']=[]     
    
        for i in range (0, len(barcode_info)):    
            if i in active_nodes:
                data_list['X'].append(barcode_info[i][1])
                data_list['Y'].append(barcode_info[i][2])
    
       
        data_list_pd = pd.DataFrame(data_list)
        chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
            alt.X('X', scale=alt.Scale(zero=False)),
            alt.Y('Y', scale=alt.Scale(zero=False)),
        )
        chart.save('/cluster/home/t116508uhn/' + args.data_name +'_tissue_altair_plot_ROI.html')
        print('Altair plot generation done: '+'/cluster/home/t116508uhn/' + args.data_name +'_tissue_altair_plot_ROI.html')    
        
        
    ##############################################################################
    # some preprocessing before making the input graph
    if args.ROI == 1:
        cells_ligand_vs_receptor = defaultdict(dict)
        ligand_list =  list(ligand_dict_dataset.keys())            
        start_index = 0 #args.slice
        end_index = len(ligand_list) #min(len(ligand_list), start_index+100)
        cell_contact_found = 0
        count_total_edges = 0
        for g in range(start_index, end_index): 
            gene = ligand_list[g]
            for i in range (0, cell_vs_gene.shape[0]): # ligand
                if i not in active_nodes:
                    continue
                if cell_vs_gene[i][gene_index[gene]] < cell_percentile[i]:
                    continue
                
                for j in range (0, cell_vs_gene.shape[0]): # receptor
                    if j not in active_nodes:
                        continue
                        
                    if i not in dist_X_dict or j not in dist_X_dict[i]: #dist_X[i,j]==0: 
                        continue
    
                    #print('%d, %d'%(i, j))
                    
                    for gene_rec in ligand_dict_dataset[gene]:
                        if cell_vs_gene[j][gene_index[gene_rec]] >= cell_percentile[j]: # or cell_vs_gene[i][gene_index[gene]] >= cell_percentile[i][4] :#gene_list_percentile[gene_rec][1]: #global_percentile: #
                            #print('i and j are both high')
                            if gene_rec in cell_cell_contact:
                                if distance_matrix[i,j] > args.spot_diameter:
                                    continue
                                else:
                                    cell_contact_found = cell_contact_found + 1
        
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
                            
            print('%d genes done out of %d ligand genes count_total_edges %d'%(g+1, len(ligand_list),count_total_edges))        
        print('total number of edges in the input graph %d with cell_contact_found %d'%(count_total_edges,cell_contact_found))
    else:
        # some preprocessing before making the input graph
        cells_ligand_vs_receptor = defaultdict(dict)
        ligand_list =  list(ligand_dict_dataset.keys())            
        start_index = 0 #args.slice
        end_index = len(ligand_list) #min(len(ligand_list), start_index+100)
        cell_contact_found = 0
        count_total_edges = 0
        for g in range(start_index, end_index): 
            gene = ligand_list[g]
            for i in range (0, cell_vs_gene.shape[0]): # ligand
                if cell_vs_gene[i][gene_index[gene]] < cell_percentile[i]:
                    continue
                
                for j in range (0, cell_vs_gene.shape[0]): # receptor
                    if i not in dist_X_dict or j not in dist_X_dict[i]: #dist_X[i,j]==0: 
                        continue

    
                    #print('%d, %d'%(i, j))
                    
                    for gene_rec in ligand_dict_dataset[gene]:
                        if cell_vs_gene[j][gene_index[gene_rec]] >= cell_percentile[j]: # or cell_vs_gene[i][gene_index[gene]] >= cell_percentile[i][4] :#gene_list_percentile[gene_rec][1]: #global_percentile: #
                            #print('i and j are both high')
                            if gene_rec in cell_cell_contact:
                                if distance_matrix[i,j] > args.spot_diameter:
                                    continue
                                else:
                                    cell_contact_found = cell_contact_found + 1
        
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
                            
            print('%d genes done out of %d ligand genes count_total_edges %d'%(g+1, len(ligand_list),count_total_edges))
        
        
        print('total number of edges in the input graph %d with cell_contact_found %d'%(count_total_edges,cell_contact_found))
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
            if i in dist_X_dict and j in dist_X_dict[i]: #dist_X[i,j]>0: #distance_matrix[i][j] <= args.neighborhood_threshold: 
                count_local = 0
                if len(cells_ligand_vs_receptor[i][j])>0:
                    for k in range (0, len(cells_ligand_vs_receptor[i][j])):
                        gene = cells_ligand_vs_receptor[i][j][k][0]
                        gene_rec = cells_ligand_vs_receptor[i][j][k][1]
                        ligand_receptor_coexpression_score = cells_ligand_vs_receptor[i][j][k][2]
                        row_col.append([i,j])
                        edge_weight.append([dist_X_dict[i][j], ligand_receptor_coexpression_score, cells_ligand_vs_receptor[i][j][k][3]])
                        lig_rec.append([gene, gene_rec])
                                                  
                        if i==j: # self-loop
                            self_loop_found[i][j] = ''
    
    if args.ROI == 1:
        total_num_cell = len(active_nodes.keys())
    else:
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

    #### needed if split data is used ##############
    if args.split>0:
        if args.ROI == 1:
            i=0
            node_id_sorted_xy=[]
            for j in range (0, len(barcode_info)):  
                if j in active_nodes:
                    node_id_sorted_xy.append([j, coordinates[j,0],coordinates[j,1]]) # 
                    i=i+1
            	
            node_id_sorted_xy = sorted(node_id_sorted_xy, key = lambda x: (x[1], x[2]))
            with gzip.open(args.metadata_to + args.data_name+'_'+'node_id_sorted_xy', 'wb') as fp:  #b, a:[0:5]   
            	pickle.dump(node_id_sorted_xy, fp)
        else:
            i=0
            node_id_sorted_xy=[]
            for cell_code in cell_id: 
                node_id_sorted_xy.append([i, coordinates[i,0],coordinates[i,1]])
                i=i+1
            	
            node_id_sorted_xy = sorted(node_id_sorted_xy, key = lambda x: (x[1], x[2]))
            with gzip.open(args.metadata_to + args.data_name+'_'+'node_id_sorted_xy', 'wb') as fp:  #b, a:[0:5]   
            	pickle.dump(node_id_sorted_xy, fp)


    print("active nodes %d"%len(node_id_sorted_xy))
    ################## required for the nest interactive version ###################
    df = pd.DataFrame(gene_ids)
    df.to_csv(args.metadata_to + 'gene_ids_'+args.data_name+'.csv', index=False, header=False)

    df = pd.DataFrame(cell_id) # after segmentation
    df.to_csv(args.metadata_to + 'cell_id_'+args.data_name+'.csv', index=False, header=False)

    df = pd.DataFrame(coordinates) # after segmentation
    df.to_csv(args.metadata_to + 'coordinates_'+args.data_name+'.csv', index=False, header=False)

    
    ######### optional #############################################################           
    # we do not need this to use anywhere. But just for debug purpose we are saving this. We can skip this if we have space issue.
    with gzip.open(args.data_to + args.data_name + '_cell_vs_gene_quantile_transformed', 'wb') as fp:  
    	pickle.dump(cell_vs_gene, fp)
    #with gzip.open(args.data_to + args.data_name + '_cell_vs_gene_quantile_transformed', 'rb') as fp:  
    #	cell_vs_gene = pickle.load(fp)
    
        
    print('write data done')
    
    end_time = datetime.datetime.now()
    print('preprocess time in seconds: ', (end_time-start_time).seconds)

    
# nohup python -u data_preprocess_NEST.py --data_name='PDAC_64630_mincell3_th98p5' --data_from='/cluster/projects/schwartzgroup/fatema/pancreatic_cancer_visium/210827_A00827_0396_BHJLJTDRXY_Notta_Karen/V10M25-61_D1_PDA_64630_Pa_P_Spatial10x_new/outs/' --filter_min_cell=3 --threshold_gene_exp=98.5 > output_data_preprocess_PDAC_64630_min_cell_3_th98p5.log &
