# Written By 
# Fatema Tuz Zohora


print('package loading')
import numpy as np
import csv
import pickle
import statistics
from scipy import sparse
import scipy.io as sio
import scanpy as sc 
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
import numpy as np
#from matplotlib.colors import LinearSegmentedColormap, to_hex, rgb2hex
#from typing import List
import qnorm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import median_abs_deviation
from scipy.stats import skew
from collections import defaultdict
import pandas as pd
import gzip
#from kneed import KneeLocator
import copy 
import argparse
import gc
import os
import altair as alt
import altairThemes # assuming you have altairThemes.py at your current directoy or your system knows the path of this altairThemes.py.


##########################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument( '--data_name', type=str, default='V1_Human_Lymph_Node_spatial', help='The name of dataset', required=True) # default='PDAC_64630',
    parser.add_argument( '--model_name', type=str, default='NEST_V1_Human_Lymph_Node_spatial', help='Name of the trained model', required=True)
    parser.add_argument( '--total_runs', type=int, default=5, help='How many runs for ensemble (at least 2 are preferred)', required=True)
    #######################################################################################################
    parser.add_argument( '--embedding_path', type=str, default='embedding_data/', help='Path to grab the attention scores from')
    parser.add_argument( '--metadata_from', type=str, default='metadata/', help='Path to grab the metadata') 
    parser.add_argument( '--data_from', type=str, default='input_graph/', help='Path to grab the input graph from (to be passed to GAT)')
    parser.add_argument( '--output_path', type=str, default='output/', help='Path to save the visualization results, e.g., histograms, graph etc.')
    parser.add_argument( '--top_percent', type=int, default=20, help='Top N percentage communications to pick')
    parser.add_argument( '--cutoff_MAD', type=int, default=-1, help='Set it to 1 to filter out communications having deviation higher than MAD')
    parser.add_argument( '--cutoff_z_score', type=float, default=-1, help='Set it to 1 to filter out communications having z_score less than 1.97 value')
    parser.add_argument( '--output_all', type=int, default=1, help='Set it to 1 to output all communications')
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
    
    
    datapoint_size = total_num_cell
    lig_rec_dict = []
    for i in range (0, datapoint_size):
        lig_rec_dict.append([])  
        for j in range (0, datapoint_size):	
            lig_rec_dict[i].append([])   
            lig_rec_dict[i][j] = []
    
    for index in range (0, len(row_col)):
            i = row_col[index][0]
            j = row_col[index][1]
            lig_rec_dict[i][j].append(lig_rec[index])  
    
    row_col = 0
    edge_weight = 0
    lig_rec = 0
    total_num_cell = 0
    
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
    for l in [3,2]: #, 3]: # 2 = layer 2, 3 = layer 1 
        layer = layer + 1
        print('layer %d'%layer)
        csv_record_dict = defaultdict(list)
        for run_time in range (start_index, start_index+total_runs):
            filename_suffix = '_'+ 'r'+str(run_time+1) +'_'
            gc.collect()
            run = run_time
            print('run %d'%run)
    
            attention_scores = []
            for i in range (0, datapoint_size):
                attention_scores.append([])   
                for j in range (0, datapoint_size):	
                    attention_scores[i].append([])   
                    attention_scores[i][j] = []
    
            distribution = []
            ##############################################
            print(args.model_name)     
            X_attention_filename = args.embedding_path +  args.model_name + filename_suffix + 'attention' #.npy
            print(X_attention_filename)
            #X_attention_bundle = np.load(X_attention_filename, allow_pickle=True) # this is deprecated
            fp = gzip.open(X_attention_filename, 'rb')  
            X_attention_bundle = pickle.load(fp)
            
            for index in range (0, X_attention_bundle[0].shape[1]):
                i = X_attention_bundle[0][0][index]
                j = X_attention_bundle[0][1][index]
                distribution.append(X_attention_bundle[l][index][0])
    
            ################# scaling the attention scores so that layer 1 and 2 will be comparable ##############################        
            min_attention_score = 1000
            max_value = np.max(distribution)
            min_value = np.min(distribution)
            print('attention score is between %g to %g, total edges %d'%(np.min(distribution), np.max(distribution), len(distribution)))
            distribution = []
            for index in range (0, X_attention_bundle[0].shape[1]):
                i = X_attention_bundle[0][0][index]
                j = X_attention_bundle[0][1][index]
                scaled_score = (X_attention_bundle[l][index][0]-min_value)/(max_value-min_value)
                attention_scores[i][j].append(scaled_score) 
                if min_attention_score > scaled_score:
                    min_attention_score = scaled_score
                distribution.append(scaled_score)
                
            print('attention score is scaled between %g to %g for ensemble'%(np.min(distribution), np.max(distribution)))
    
            #print('min attention score %g, total edges %d'%(min_attention_score, len(distribution))) 
            # should always print 0 for min attention score
            
            ccc_index_dict = dict()
            threshold_down =  np.percentile(sorted(distribution), 0)
            threshold_up =  np.percentile(sorted(distribution), 100)
            connecting_edges = np.zeros((len(barcode_info),len(barcode_info)))
            for j in range (0, datapoint_size):
                #threshold =  np.percentile(sorted(attention_scores[:,j]), 97) #
                for i in range (0, datapoint_size):
                    atn_score_list = attention_scores[i][j]
                    #print(len(atn_score_list))
                    #s = min(0,len(atn_score_list)-1)
                    for k in range (0, len(atn_score_list)):
                        if attention_scores[i][j][k] >= threshold_down and attention_scores[i][j][k] <= threshold_up: #np.percentile(sorted(distribution), 50):
                            connecting_edges[i][j] = 1
                            ccc_index_dict[i] = ''
                            ccc_index_dict[j] = ''
        
        
        
            graph = csr_matrix(connecting_edges)
            n_components, labels = connected_components(csgraph=graph,directed=True, connection = 'weak',  return_labels=True) #
            #print('number of component %d'%n_components)
        
            count_points_component = np.zeros((n_components))
            for i in range (0, len(labels)):
                 count_points_component[labels[i]] = count_points_component[labels[i]] + 1
        
            #print(count_points_component)
        
            id_label = 2 # initially all are zero. =1 those who have self edge but above threshold. >= 2 who belong to some component
            index_dict = dict()
            for i in range (0, count_points_component.shape[0]):
                if count_points_component[i]>1:
                    index_dict[i] = id_label
                    id_label = id_label+1
        
            #print('number of components with multiple datapoints is %d'%id_label)
        
        
            for i in range (0, len(barcode_info)):
            #    if barcode_info[i][0] in barcode_label:
                if count_points_component[labels[i]] > 1:
                    barcode_info[i][3] = index_dict[labels[i]] #2
                elif connecting_edges[i][i] == 1 and len(lig_rec_dict[i][i])>0: 
                    barcode_info[i][3] = 1
                else:
                    barcode_info[i][3] = 0
        
            #######################
        
        
     
            ###############
            csv_record = []
            csv_record.append(['from_cell', 'to_cell', 'ligand', 'receptor', 'attention_score', 'component', 'from_id', 'to_id'])
            for j in range (0, len(barcode_info)):
                for i in range (0, len(barcode_info)):
                    
                    if i==j:
                        if len(lig_rec_dict[i][j])==0:
                            continue # because the model generates some score for selfloop even if there was no user input
                                    # to prevent losing own information
                     
                    atn_score_list = attention_scores[i][j]
                    for k in range (0, len(atn_score_list)):
                        if attention_scores[i][j][k] >= threshold_down and attention_scores[i][j][k] <= threshold_up: 
                            if barcode_info[i][3]==0:
                                print('error')
                            elif barcode_info[i][3]==1: # selfloop = autocrine
                                csv_record.append([barcode_info[i][0], barcode_info[j][0], lig_rec_dict[i][j][k][0], lig_rec_dict[i][j][k][1], attention_scores[i][j][k], '0-single', i, j])
                            else: # paracrine or juxtacrine
                                csv_record.append([barcode_info[i][0], barcode_info[j][0], lig_rec_dict[i][j][k][0], lig_rec_dict[i][j][k][1], attention_scores[i][j][k], barcode_info[i][3], i, j])
     
            ###########	
          
            #print('records found %d'%len(csv_record))
            for i in range (1, len(csv_record)): 
                key_value = str(csv_record[i][6]) +'-'+ str(csv_record[i][7]) + '-' + csv_record[i][2] + '-' + csv_record[i][3]
                # i-j-ligandGene-receptorGene
                csv_record_dict[key_value].append([csv_record[i][4], run])
                
        '''    
        for key_value in csv_record_dict.keys():
            run_dict = defaultdict(list)
            for scores in csv_record_dict[key_value]: # entry count = total_runs 
                run_dict[scores[1]].append(scores[0]) # [run_id]=score
            
            for runs in run_dict.keys():
                run_dict[runs] = np.mean(run_dict[runs]) # taking the mean attention score
            
     
            csv_record_dict[key_value] = [] # make it blank
            for runs in run_dict.keys(): # has just one mean value for the attention score
                csv_record_dict[key_value].append([run_dict[runs],runs]) # [score, 0]
        '''
        ########## All runs combined. Now find rank product #############################
        
        all_edge_list = []
        for key_value in csv_record_dict.keys():
            edge_score_runs = []
            edge_score_runs.append(key_value)
            for runs in csv_record_dict[key_value]:
                edge_score_runs.append(runs[0]) #
            
            all_edge_list.append(edge_score_runs) # [[key_value, score_by_run1, score_by_run2, etc.],...]
    
        ## Find the rank product #####################################################################
        ## all_edge_list has all the edges along with their scores for different runs in the following format: 
        ## [edge_1_info, score_by_run1, score_by_run2, etc.], [edge_2_info, score_by_run1, score_by_run2, etc.], ..., [edge_N_info, score_by_run1, score_by_run2, etc.]
        edge_rank_dictionary = defaultdict(list)
        # sort the all_edge_list by each run's rank 
        #print('total runs %d. Ensemble them.'%total_runs)
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
            avg_score = 0 
            total_weight = 0
            for i in range (0, len(edge_rank_dictionary[key_val])):
                rank_product = rank_product * edge_rank_dictionary[key_val][i]
                score_product = score_product * (attention_score_list[i][0]+0.01) 
                # translated by a tiny amount to avoid producing 0 during product 
                weight_by_run = max_weight - edge_rank_dictionary[key_val][i]
                avg_score = avg_score + attention_score_list[i][0] * weight_by_run
                total_weight = total_weight + weight_by_run
                
            avg_score = avg_score/total_weight # lower weight being higher attention np.max(avg_score) #
            all_edge_vs_rank.append([key_val, rank_product**(1/total_runs), score_product**(1/total_runs)])  # small rank being high attention
            # or all_edge_vs_rank.append([key_val, rank_product**(1/total_runs), avg_score])
            distribution_rank[layer].append(rank_product**(1/total_runs))
            distribution_score[layer].append(score_product**(1/total_runs)) #avg_score)

        all_edge_sorted_by_rank[layer] = sorted(all_edge_vs_rank, key = lambda x: x[1]) # small rank being high attention 
        #print('rank ranges from %g to %g'%(np.min(distribution_rank[layer]), np.max(distribution_rank[layer])))
        #print('score ranges from %g to %g'%(np.min(distribution_score[layer]), np.max(distribution_score[layer]))) 
    #############################################################################################################################################
    # for each layer, I scale the attention scores [0, 1] over all the edges. So that they are comparable or mergeable between layers
    # for each edge, we have two sets of (rank, score) due to 2 layers. We take union of them.
    print('Multiple runs for each layer are ensembled.') 
    for layer in range (0, 2):
        score_min = np.min(distribution_score[layer])
        score_max = np.max(distribution_score[layer])
        rank_min = np.min(distribution_rank[layer])
        rank_max = np.max(distribution_rank[layer])
        distribution_rank[layer] = []
        distribution_score[layer] = []
        #a = 1
        #b = len(all_edge_sorted_by_rank[layer])
        #print('b %d'%b)
        for i in range (0, len(all_edge_sorted_by_rank[layer])):
            # score is scaled between 0 to 1 again for easier interpretation 
            all_edge_sorted_by_rank[layer][i][2] = (all_edge_sorted_by_rank[layer][i][2]-score_min)/(score_max-score_min)
            all_edge_sorted_by_rank[layer][i][1] = i+1 # done for easier interpretation
            #((all_edge_sorted_by_rank[layer][i][1]-rank_min)/(rank_max-rank_min))*(b-a) + a
            distribution_rank[layer].append(all_edge_sorted_by_rank[layer][i][1])
            distribution_score[layer].append(all_edge_sorted_by_rank[layer][i][2])
    
    ################################ or ###############################################################################################################
    if args.cutoff_MAD ==-1 and args.cutoff_z_score == -1:
        percentage_value = args.top_percent #20 ##100 #20 # top 20th percentile rank, low rank means higher attention score
        csv_record_intersect_dict = defaultdict(list)
        edge_score_intersect_dict = defaultdict(list)
        for layer in range (0, 2):
            threshold_up = np.percentile(distribution_rank[layer], percentage_value) #np.round(np.percentile(distribution_rank[layer], percentage_value),2)
            for i in range (0, len(all_edge_sorted_by_rank[layer])):
                if all_edge_sorted_by_rank[layer][i][1] <= threshold_up: # because, lower rank means higher strength
                    csv_record_intersect_dict[all_edge_sorted_by_rank[layer][i][0]].append(all_edge_sorted_by_rank[layer][i][1]) #i+1) # already sorted by rank. so just use i as the rank 
                    edge_score_intersect_dict[all_edge_sorted_by_rank[layer][i][0]].append(all_edge_sorted_by_rank[layer][i][2]) # score
        ###########################################################################################################################################
        ## get the aggregated rank for all the edges ##
        distribution_temp = []
        for key_value in csv_record_intersect_dict.keys():  
            arg_index = np.argmin(csv_record_intersect_dict[key_value]) # layer 0 or 1, whose rank to use # should I take the avg rank instead, and scale the ranks (1 to count(total_edges)) later? 
            csv_record_intersect_dict[key_value] = np.min(csv_record_intersect_dict[key_value]) # use that rank. smaller rank being the higher attention
            edge_score_intersect_dict[key_value] = edge_score_intersect_dict[key_value][arg_index] # use that score
            distribution_temp.append(csv_record_intersect_dict[key_value]) 
        
        #################
        
        ################################################################################
        csv_record_dict = copy.deepcopy(csv_record_intersect_dict)
        
        ################################################################################
            
        combined_score_distribution = []
        csv_record = []
        csv_record.append(['from_cell', 'to_cell', 'ligand', 'receptor', 'edge_rank', 'component', 'from_id', 'to_id', 'attention_score'])
        for key_value in csv_record_dict.keys():
            item = key_value.split('-')
            i = int(item[0])
            j = int(item[1])
            ligand = item[2]
            receptor = item[3]        
            edge_rank = csv_record_dict[key_value]        
            score = edge_score_intersect_dict[key_value] # weighted average attention score, where weight is the rank, lower rank being higher attention score
            label = -1 
            csv_record.append([barcode_info[i][0], barcode_info[j][0], ligand, receptor, edge_rank, label, i, j, score])
            combined_score_distribution.append(score)
        
                
        print('common LR count %d'%len(csv_record))
        
        
        ##### save the file for downstream analysis ########
        csv_record_final = []
        csv_record_final.append(csv_record[0])
        for k in range (1, len(csv_record)):
            ligand = csv_record[k][2]
            receptor = csv_record[k][3]
            #if ligand =='CCL19' and receptor == 'CCR7':
            csv_record_final.append(csv_record[k])
        
        
            
        df = pd.DataFrame(csv_record_final) # output 4
        df.to_csv(args.output_path + args.model_name+'_top' + str(args.top_percent) + 'percent.csv', index=False, header=False)
        print('result is saved at: '+args.output_path + args.model_name+'_top' + str(args.top_percent) + 'percent.csv')
############### skewness plot ##############
    percentage_value = 100 #20 ##100 #20 # top 20th percentile rank, low rank means higher attention score
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
    csv_record_dict = copy.deepcopy(csv_record_intersect_dict)
    
    ################################################################################
    score_distribution = []
    csv_record = []
    csv_record.append(['from_cell', 'to_cell', 'ligand', 'receptor', 'edge_rank', 'component', 'from_id', 'to_id', 'attention_score']) #, 'deviation_from_median'
    for key_value in csv_record_dict.keys():
        item = key_value.split('-')
        i = int(item[0])
        j = int(item[1])
        ligand = item[2]
        receptor = item[3]        
        edge_rank = csv_record_dict[key_value]        
        score = edge_score_intersect_dict[key_value] # weighted average attention score, where weight is the rank, lower rank being higher attention score
        label = -1 
        csv_record.append([barcode_info[i][0], barcode_info[j][0], ligand, receptor, edge_rank, label, i, j, score])
        score_distribution.append(score)
    
            
#    print('common LR count %d'%len(csv_record))
    
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
    skewness_distribution = skew(score_distribution)
    
    if args.output_all == 1:
        df = pd.DataFrame(csv_record) # 
        df.to_csv(args.output_path + args.model_name+'_allCCC.csv', index=False, header=False)
        print(args.output_path + args.model_name+'_allCCC.csv')
    ###########
    if args.cutoff_MAD !=-1:
        MAD = median_abs_deviation(score_distribution)
        median_distribution = statistics.median(score_distribution)
        csv_record_final = []
        csv_record_final.append(csv_record[0])
        csv_record_final[0].append('deviation_from_median')
        for k in range (1, len(csv_record)):
            deviation_from_median = median_distribution-csv_record[k][8]
            if deviation_from_median <= MAD:   
                temp_record = csv_record[k]
                temp_record.append(deviation_from_median)
                csv_record_final.append(temp_record)
                
    
        df = pd.DataFrame(csv_record_final) # 
        df.to_csv(args.output_path + args.model_name+'_MAD_cutoff.csv', index=False, header=False)
        
    ##### save the file for downstream analysis ########
    if args.cutoff_z_score !=-1:
        z_score_distribution = stats.zscore(score_distribution)
        csv_record_final = []
        csv_record_final.append(csv_record[0])
        csv_record_final[0].append('z-score')
        for k in range (1, len(csv_record)):
            if z_score_distribution[k-1] >= 1.97: #args.cutoff_z_score:  
                temp_record = csv_record[k]
                temp_record.append(z_score_distribution[k-1])
                csv_record_final.append(temp_record)
    
        df = pd.DataFrame(csv_record_final) # output 4
        df.to_csv(args.output_path + args.model_name+'_z_score_cutoff.csv', index=False, header=False)
        
    ###########################################################################################################################################
    
    # plot the distribution    

    # write the skewness and MAD value in a text file
    # Opening a file
    '''
    file1 = open(args.output_path + args.model_name+'_statistics.txt', 'w')
    L = ["Median Absolute Deviation (MAD):"+str(MAD)+"\n", "Skewness: "+str(skewness_distribution) +" \n"]    
    # Writing multiple strings
    file1.writelines(L)
    # Closing file
    file1.close()
    ''' 
    

 
