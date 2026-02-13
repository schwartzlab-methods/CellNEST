# Written By 
# Fatema Tuz Zohora

print('package loading')
import numpy as np
import csv
import pickle
from scipy import sparse
import scipy.io as sio
import scanpy as sc 
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_hex, rgb2hex
#from typing import List
import qnorm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from collections import defaultdict
import pandas as pd
import gzip
#from kneed import KneeLocator
import copy 
import argparse
import gc
import os
#from scipy.sparse import dok_matrix
import altair as alt

##########################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument( '--data_name', type=str, default='Visium_HD_Human_Colon_Cancer_square_002um_outputs', help='The name of dataset') #  required=True
    parser.add_argument( '--model_name', type=str, default='CellNEST_Visium_HD_Human_Colon_Cancer_square_002um_outputs', help='Name of the trained model') #, required=True
    parser.add_argument( '--total_runs', type=int, default=3, help='How many runs for ensemble (at least 2 are preferred)') #, required=True
    parser.add_argument( '--total_subgraphs', type=int, default=15)
    #######################################################################################################
    parser.add_argument( '--embedding_path', type=str, default='embedding_data/', help='Path to grab the attention scores from')
    parser.add_argument( '--metadata_from', type=str, default='metadata/', help='Path to grab the metadata') 
    parser.add_argument( '--data_from', type=str, default='input_graph/', help='Path to grab the input graph from (to be passed to GAT)')
    parser.add_argument( '--output_path', type=str, default='output/', help='Path to save the visualization results, e.g., histograms, graph etc.')
    parser.add_argument( '--top_percent', type=int, default=20, help='Top N percentage communications to pick')
    parser.add_argument( '--output_all', type=int, default=1, help='Set it to 1 to output all communications')
    parser.add_argument( '--integrate_layers', type=int, default=1, help='Set it to 1 to output strictly top percent number of edges.')
    args = parser.parse_args()

    args.metadata_from = args.metadata_from + args.data_name + '/'
    args.data_from = args.data_from + args.data_name + '/'
    args.embedding_path  = args.embedding_path + args.data_name + '/'
    args.output_path = args.output_path + args.data_name + '/'
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

##################### get metadata: barcode_info ###################################
 

    with gzip.open(args.metadata_from +args.data_name+'_barcode_info', 'rb') as fp:  #b, a:[0:5]   _filtered
        barcode_info = pickle.load(fp)
    
    with gzip.open(args.data_from + args.data_name + '_adjacency_records', 'rb') as fp:  #b, a:[0:5]  _filtered 
        row_col, edge_weight, lig_rec, total_num_cell = pickle.load(fp)
    
    
    
    lig_rec_dict = defaultdict(dict)    
    for index in range (0, len(row_col)):
            i = row_col[index][0]
            j = row_col[index][1]
            if i in lig_rec_dict:
                if j in lig_rec_dict[i]:
                    lig_rec_dict[i][j].append(lig_rec[index]) 
                else:
                    lig_rec_dict[i][j] = []
                    lig_rec_dict[i][j].append(lig_rec[index])
            else:
                lig_rec_dict[i][j] = []
                lig_rec_dict[i][j].append(lig_rec[index])                

        
             
    
    ################################################################################################################
    dict_cell_edge = defaultdict(list) # key = node. values = incoming edges
    dict_cell_neighbors = defaultdict(list) # key = node. value = nodes corresponding to incoming edges/neighbors
    #nodes_active = dict()
    for i in range(0, len(row_col)):
        dict_cell_edge[row_col[i][1]].append(i) # index of the edges
        dict_cell_neighbors[row_col[i][1]].append(row_col[i][0]) # neighbor id
        #nodes_active[row_col[i][1]] = '' # to
        #nodes_active[row_col[i][0]] = '' # from


    datapoint_size = total_num_cell #len(nodes_active.keys())

    for i in range (0, datapoint_size):
        neighbor_list = dict_cell_neighbors[i]
        neighbor_list = list(set(neighbor_list))
        dict_cell_neighbors[i] = neighbor_list

    
    node_id_sorted = args.metadata_from + args.data_name+'_'+'node_id_sorted_xy'
    fp = gzip.open(node_id_sorted, 'rb')
    node_id_sorted_xy = pickle.load(fp)

    #node_id_sorted_xy_temp = []
    unfiltered_index_to_filtered_serial = dict() # not the serial. We do not need to know the serial
    filtered_serial_to_unfiltered_index = dict()
    active_node_count = 0
    for i in range(0, len(node_id_sorted_xy)):
        #if node_id_sorted_xy[i][0] in nodes_active: # skip those which are not in our ROI
            #node_id_sorted_xy_temp.append(node_id_sorted_xy[i])
        unfiltered_index_to_filtered_serial[node_id_sorted_xy[i][0]] =  active_node_count
        filtered_serial_to_unfiltered_index[active_node_count] =  node_id_sorted_xy[i][0]
        active_node_count = active_node_count + 1

    #node_id_sorted_xy = node_id_sorted_xy_temp
    ##################################################################################################################
    # split it into N set of edges

    total_subgraphs = args.total_subgraphs

    edge_list = []

    start_index = []
    id_map_old_new = [] # make an index array, so that existing node ids are mapped to new ids
    id_map_new_old = []

    for i in range (0, total_subgraphs+1):
        start_index.append((datapoint_size//total_subgraphs)*i)
        id_map_old_new.append(dict())
        id_map_new_old.append(dict())

    set_id=-1
    for indx in range (0, len(start_index)-1):
        set_id = set_id + 1
        #print('graph id %d, node %d to %d'%(set_id,start_index[indx],start_index[indx+1]))
        set1_nodes = []
        set1_edges_index = []
        node_limit_set1 = start_index[indx+1]
        set1_direct_edges = []

        for i in range (start_index[indx], node_limit_set1):
            set1_nodes.append(node_id_sorted_xy[i][0])
            # add it's edges - first hop

            for edge_index in dict_cell_edge[node_id_sorted_xy[i][0]]:
                set1_edges_index.append(edge_index) # has both row_col and edge_weight
                set1_direct_edges.append(edge_index)
            # add it's neighbor's edges - second hop
            for neighbor in dict_cell_neighbors[node_id_sorted_xy[i][0]]:
                if node_id_sorted_xy[i][0] == neighbor:
                    continue
                for edge_index in dict_cell_edge[neighbor]:
                    set1_edges_index.append(edge_index) # has both row_col and edge_weight

        set1_edges_index = list(set(set1_edges_index))

        #print('len of set1_edges_index %d'%len(set1_edges_index))
        #if len(set1_edges_index)==0:
        #    break

        # old to new mapping of the nodes
        # make an index array, so that existing node ids are mapped to new ids
        new_id = 0
        #spot_list = []
        for k in set1_edges_index:
            i = row_col[k][0] # unfiltered node index
            j = row_col[k][1] # unfiltered node index
            if i not in id_map_old_new[set_id]:
                id_map_old_new[set_id][i] = new_id # old = unfiltered, new = filtered + subgraph specific
                id_map_new_old[set_id][new_id] = i
                #spot_list.append(i) #new_id)
                new_id = new_id + 1

            if j not in id_map_old_new[set_id]:
                id_map_old_new[set_id][j] = new_id
                id_map_new_old[set_id][new_id] = j
                #spot_list.append(j) #new_id)
                new_id = new_id + 1


        #print('new id: %d'%new_id)
        set1_edges = []
        for i in set1_direct_edges:  #set1_edges_index:
            set1_edges.append([[id_map_old_new[set_id][row_col[i][0]], id_map_old_new[set_id][row_col[i][1]]], edge_weight[i]])
            #set1_edges.append([row_col[i], edge_weight[i]])

        edge_list.append(set1_edges)
        num_cell = new_id

        print("subgraph %d: number of nodes %d. Total number of edges %d"%(set_id, num_cell, len(set1_edges)))

        gc.collect()

    
    ############# load output graph #################################################
    
    #filename_suffix = ["_r1_", "r2_", "r3_", "r4_", "r5_", "r6_", "r7_", "r8_", "r9_", "r10_"]
    total_runs = args.total_runs 
    start_index = 0 
    distribution_rank = []
    distribution_score = []
    all_edge_sorted_by_rank = []
    for layer in range (0, 2):
        distribution_rank.append([])
        distribution_score.append([])
        all_edge_sorted_by_rank.append([])
    
    layer = -1
    for l in [3, 2]: #, 3]: # 2 = layer 2, 3 = layer 1 
        layer = layer + 1
        print('layer %d'%layer)
        csv_record_dict = defaultdict(list)
        for run_time in range (start_index, start_index+total_runs):
            filename_suffix = '_'+ 'r'+str(run_time+1) +'_' #str(run_time+1) +'_'
            gc.collect()
            run = run_time
            print('run %d'%run)


                            
            attention_scores = defaultdict(dict)      
            distribution = []
            ##########################
            print(args.model_name) 
            #######################################################################
            for set_id in range(0, len(edge_list)):
                print('subgraph %d'%set_id)
                ##############
                set1_exist_dict = defaultdict(dict)        
                for edge in edge_list[set_id]:
                    row_col = edge[0]
                    new_i = row_col[0] 
                    new_j = row_col[1]
                    i = id_map_new_old[set_id][new_i] # unfiltered
                    j = id_map_new_old[set_id][new_j] # unfiltered
                    set1_exist_dict[i][j] = 1
                
                ############
                 
                X_attention_filename = args.embedding_path +  args.model_name + filename_suffix + 'attention' + '_subgraph'+str(set_id)
                fp = gzip.open(X_attention_filename, 'rb')  
                X_attention_bundle = pickle.load(fp)
                print(X_attention_filename)  
                edge_found = 0
                for index in range (0, X_attention_bundle[0].shape[1]):
                    new_i = X_attention_bundle[0][0][index]
                    new_j = X_attention_bundle[0][1][index] 
                    # these returned i and j are new
                    i = id_map_new_old[set_id][new_i] # unfiltered index
                    j = id_map_new_old[set_id][new_j] # unfiltered index          
                                
                    if i in set1_exist_dict and j in set1_exist_dict[i] and set1_exist_dict[i][j]==1:
                    ###################################
                        split_i = unfiltered_index_to_filtered_serial[i] 
                        split_j = unfiltered_index_to_filtered_serial[j]
                        if split_i not in attention_scores or split_j not in attention_scores[split_i]:
                            attention_scores[split_i][split_j]=[]
                        
                        attention_scores[split_i][split_j].append(X_attention_bundle[l][index][0]) 
                        distribution.append(X_attention_bundle[l][index][0])
                        edge_found = edge_found + 1
                print('Edge found %d out of %d'%(edge_found, X_attention_bundle[0].shape[1]))
                        
            gc.collect()
            #######################    
            print('All subgraph load done')
            
            ################# scaling the attention scores so that layer 1 and 2 will be comparable ##############################        
            min_attention_score = 1000
            max_value = np.max(distribution)
            min_value = np.min(distribution)
            print('attention score is between %g to %g, total edges %d'%(np.min(distribution), np.max(distribution), len(distribution)))
            distribution = []
            for i  in attention_scores:
                for j in attention_scores[i]:
                       
                    for k in range (0, len(attention_scores[i][j])):
                        attention_scores[i][j][k] = (attention_scores[i][j][k]-min_value)/(max_value-min_value)
                        scaled_score = attention_scores[i][j][k]
                        if min_attention_score > scaled_score:
                            min_attention_score = scaled_score
                            
                        distribution.append(scaled_score)

            print('attention score is scaled between %g to %g for ensemble'%(np.min(distribution), np.max(distribution)))

            # print('PASSED')
            ###############

            csv_record = []
            csv_record.append(['from_cell', 'to_cell', 'ligand', 'receptor', 'attention_score', 'component', 'from_id', 'to_id'])
            for i in lig_rec_dict:
                for j in lig_rec_dict[i]:                           
                    split_i = unfiltered_index_to_filtered_serial[i] 
                    split_j = unfiltered_index_to_filtered_serial[j]
                    if split_i not in attention_scores or split_j not in attention_scores[split_i]:
                        continue
 
                    atn_score_list = attention_scores[split_i][split_j]
                    for k in range (0, len(atn_score_list)):
                        csv_record.append([barcode_info[i][0], barcode_info[j][0], lig_rec_dict[i][j][k][0], lig_rec_dict[i][j][k][1], attention_scores[split_i][split_j][k], barcode_info[i][3], i, j])
     
            ###########	

          
            print('records found %d'%len(csv_record))
            for i in range (1, len(csv_record)): 
                key_value = str(csv_record[i][6]) +'+'+ str(csv_record[i][7]) + '+' + csv_record[i][2] + '+' + csv_record[i][3]
                csv_record_dict[key_value].append([csv_record[i][4], run])

        ##### one run completes #####
        
        all_edge_list = []
        for key_value in csv_record_dict.keys():
            edge_score_runs = []
            edge_score_runs.append(key_value)
            for runs in csv_record_dict[key_value]:
                edge_score_runs.append(runs[0]) #
            
            all_edge_list.append(edge_score_runs) # [[key_value, score_by_run1, score_by_run2, etc.],...]
    
        ## Find the rank product #####################################################################
        ## all_edge_list has all the edges along with their scores for different runs in following format: 
        ## [edge_1_info, score_by_run1, score_by_run2, etc.], [edge_2_info, score_by_run1, score_by_run2, etc.], ..., [edge_N_info, score_by_run1, score_by_run2, etc.]
        edge_rank_dictionary = defaultdict(list)
        # sort the all_edge_list by each run's rank 
        print('total runs %d'%total_runs)
        for runs in range (0, total_runs):
            sorted_list_temp = sorted(all_edge_list, key = lambda x: x[runs+1], reverse=True) # sort based on attention score by current run: large to small
            for rank in range (0, len(sorted_list_temp)):
                edge_rank_dictionary[sorted_list_temp[rank][0]].append(rank+1) # small rank being high attention, starting from 1
                
        max_weight = len(all_edge_list) + 1 # maximum possible rank 
        all_edge_vs_rank = []
        for key_val in edge_rank_dictionary.keys():
            rank_product = 1
            score_product = 1
            attention_score_list = csv_record_dict[key_val] # [[score, run_id],...]
            avg_score = 0 #[]
            total_weight = 0
            for i in range (0, len(edge_rank_dictionary[key_val])):
                rank_product = rank_product * edge_rank_dictionary[key_val][i]
                score_product = score_product * (attention_score_list[i][0]+0.01) 
                # translated by a tiny amount to avoid producing 0 during product                
                weight_by_run = max_weight - edge_rank_dictionary[key_val][i]
                avg_score = avg_score + attention_score_list[i][0] * weight_by_run
                #avg_score.append(attention_score_list[i][0])  
                total_weight = total_weight + weight_by_run
                
            avg_score = avg_score/total_weight # lower weight being higher attention np.max(avg_score) #
            all_edge_vs_rank.append([key_val, rank_product**(1/total_runs), score_product**(1/total_runs)]) #avg_score 
            # small rank being high attention
            distribution_rank[layer].append(rank_product**(1/total_runs))
            distribution_score[layer].append(score_product**(1/total_runs)) #avg_score)
            
        all_edge_sorted_by_rank[layer] = sorted(all_edge_vs_rank, key = lambda x: x[1]) # small rank being high attention 
    
    #############################################################################################################################################
    # for each layer, I scale the attention scores [0, 1] over all the edges. 
    # So that they are comparable or mergeable between layers
    # for each edge, we have two sets of (rank, score) due to 2 layers. 
    # We take union of them.
    print('Multiple runs for each layer are ensembled.') 
    for layer in range (0, 2):
        score_min = np.min(distribution_score[layer])
        score_max = np.max(distribution_score[layer])
        rank_min = np.min(distribution_rank[layer])
        rank_max = np.max(distribution_rank[layer])
        distribution_rank[layer] = []
        distribution_score[layer] = []
        for i in range (0, len(all_edge_sorted_by_rank[layer])):
            # score is scaled between 0 to 1 again for easier interpretation 
            all_edge_sorted_by_rank[layer][i][2] = (all_edge_sorted_by_rank[layer][i][2]-score_min)/(score_max-score_min)
            all_edge_sorted_by_rank[layer][i][1] = i+1 # done for easier interpretation
            #((all_edge_sorted_by_rank[layer][i][1]-rank_min)/(rank_max-rank_min))*(b-a) + a
            distribution_rank[layer].append(all_edge_sorted_by_rank[layer][i][1])
            distribution_score[layer].append(all_edge_sorted_by_rank[layer][i][2])
    
    ################################ now take the cut-off ###############################################################################################################
    if args.integrate_layers == 0: 
        # top N percent of both layers will be taken and merged by union 
        # - may result in more than 20% edges of total edges since it takes union
        percentage_value = args.top_percent #20 ##100 #20 # top 20th percentile rank, low rank means higher attention score
        csv_record_intersect_dict = defaultdict(list)
        edge_score_intersect_dict = defaultdict(list)
        for layer in range (0, 2):
            threshold_up = np.percentile(distribution_rank[layer], percentage_value) #np.round(np.percentile(distribution_rank[layer], percentage_value),2)
            for i in range (0, len(all_edge_sorted_by_rank[layer])):
                if all_edge_sorted_by_rank[layer][i][1] <= threshold_up: # because, lower rank means higher strength
                    csv_record_intersect_dict[all_edge_sorted_by_rank[layer][i][0]].append(all_edge_sorted_by_rank[layer][i][1]) # rank
                    edge_score_intersect_dict[all_edge_sorted_by_rank[layer][i][0]].append(all_edge_sorted_by_rank[layer][i][2]) # score
        ###########################################################################################################################################
        ## get the aggregated rank for all the edges ##
        distribution_temp = []
        for key_value in csv_record_intersect_dict.keys():  
            arg_index = np.argmin(csv_record_intersect_dict[key_value]) # layer 0 or 1, whose rank to use 
            csv_record_intersect_dict[key_value] = np.min(csv_record_intersect_dict[key_value]) # use that rank. smaller rank being the higher attention
            edge_score_intersect_dict[key_value] = edge_score_intersect_dict[key_value][arg_index] # use that score
            distribution_temp.append(csv_record_intersect_dict[key_value]) 
        
        #################
        
        ################################################################################
        csv_record_dict = csv_record_intersect_dict
        csv_record = []
        csv_record.append(['from_cell', 'to_cell', 'ligand', 'receptor', 'edge_rank', 'component', 'from_id', 'to_id', 'attention_score'])
        for key_value in csv_record_dict.keys():
            item = key_value.split('+')
            i = int(item[0])
            j = int(item[1])
            ligand = item[2]
            receptor = item[3]        
            edge_rank = csv_record_dict[key_value]        
            score = edge_score_intersect_dict[key_value] # weighted average attention score, where weight is the rank, lower rank being higher attention score
            label = -1 
            csv_record.append([barcode_info[i][0], barcode_info[j][0], ligand, receptor, edge_rank, label, i, j, score])

        print('Strong LR count %d'%len(csv_record))    
        ##### save the file for downstream analysis ########        
        df = pd.DataFrame(csv_record) # output 4
        df.to_csv(args.output_path + args.model_name+'_top' + str(args.top_percent) + 'percent.csv', index=False, header=False)
        print('Result saved: '+args.output_path + args.model_name+'_top' + str(args.top_percent) + 'percent.csv')

    #############################################################################
    percentage_value = 100 #low rank means higher attention score
    csv_record_intersect_dict = defaultdict(list)
    edge_score_intersect_dict = defaultdict(list)
    # union of layers
    for layer in range (0, 2):
        for i in range (0, len(all_edge_sorted_by_rank[layer])):
            csv_record_intersect_dict[all_edge_sorted_by_rank[layer][i][0]].append(all_edge_sorted_by_rank[layer][i][1]) # rank 
            edge_score_intersect_dict[all_edge_sorted_by_rank[layer][i][0]].append(all_edge_sorted_by_rank[layer][i][2]) # score
    ###########################################################################################################################################
    ## get one aggregated rank and score for all the edges ##
    key_vs_rank = []
    for key_value in csv_record_intersect_dict.keys():  
        #arg_index = np.argmin(csv_record_intersect_dict[key_value]) # layer 0 or 1, whose rank to use  
        csv_record_intersect_dict[key_value] = np.mean(csv_record_intersect_dict[key_value]) #np.min(csv_record_intersect_dict[key_value]) # use that rank. smaller rank being the higher attention
        key_vs_rank.append([key_value, csv_record_intersect_dict[key_value]])
        edge_score_intersect_dict[key_value] = np.mean(edge_score_intersect_dict[key_value]) #edge_score_intersect_dict[key_value][arg_index] # use that score
         
    ###############################################################################################################
    ## Ranks are now floating point. Make them integer.  
    key_vs_rank = sorted(key_vs_rank, key = lambda x: x[1]) # small rank being high attention 
    csv_record_intersect_dict = dict()
    for i in range (0, len(key_vs_rank)):
        csv_record_intersect_dict[key_vs_rank[i][0]] = i + 1 # ranking starts from index 1

    ################################################################################
    csv_record_dict = copy.deepcopy(csv_record_intersect_dict)
    
    ################################################################################
    score_distribution = []
    csv_record = []
    csv_record.append(['from_cell', 'to_cell', 'ligand', 'receptor', 'edge_rank', 'component', 'from_id', 'to_id', 'attention_score']) #, 'deviation_from_median'
    for key_value in csv_record_dict.keys():
        item = key_value.split('+')
        i = int(item[0])
        j = int(item[1])
        ligand = item[2]
        receptor = item[3]        
        edge_rank = csv_record_dict[key_value]        
        score = edge_score_intersect_dict[key_value] # weighted average attention score, where weight is the rank, lower rank being higher attention score
        label = -1 
        csv_record.append([barcode_info[i][0], barcode_info[j][0], ligand, receptor, edge_rank, label, i, j, score])
        score_distribution.append(score)
    
            

    print('Total LR edge count %d, top %g percent LR edge count %d'%(len(csv_record), args.top_percent, (len(csv_record)*args.top_percent)//100))    
    print('Result saved: ')
    data_list=dict()
    data_list['attention_score']=[]
    for score in score_distribution:
        data_list['attention_score'].append(score)
        
    df = pd.DataFrame(data_list)    
    chart = alt.Chart(df).transform_density(
            'attention_score',
            as_=['attention_score', 'density'],
        ).mark_area().encode(
            x="attention_score:Q",
            y='density:Q',
        )

    chart.save(args.output_path + args.model_name+'_attention_score_distribution.html')  
    print(args.output_path + args.model_name+'_attention_score_distribution.html')
    #skewness_distribution = skew(score_distribution)
    #print('skewness of the distribution is %g'%skewness_distribution)
    if args.output_all == 1:
        df = pd.DataFrame(csv_record) # 
        df.to_csv(args.output_path + args.model_name+'_allCCC.csv', index=False, header=False)
        print(args.output_path + args.model_name+'_allCCC.csv')

    if args.integrate_layers == 1:
        # take top 20% highly ranked edges 
        csv_record[1:] = sorted(csv_record[1:], key = lambda x: x[4]) # small rank being high attention
        csv_record = csv_record[0: (len(csv_record)*args.top_percent)//100]
        ##### save the file for downstream analysis ########        
        df = pd.DataFrame(csv_record) # output 4
        df.to_csv(args.output_path + args.model_name+'_top' + str(args.top_percent) + 'percent.csv', index=False, header=False)
        print(args.output_path + args.model_name+'_top' + str(args.top_percent) + 'percent.csv')
