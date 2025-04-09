# Written By
# Fatema Tuz Zohora

print('package loading')
import numpy as np
import csv
import pickle
import statistics
import numpy as np
from collections import defaultdict
import pandas as pd
import gzip
import copy 
import argparse
import gc
import os
import subprocess

# load the CellNEST detected results
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_name', type=str, help='The name of dataset', required=True) # 
    parser.add_argument( '--model_name', type=str, help='Name of the trained model', required=True)
    parser.add_argument( '--top_percent', type=int, default=20, help='Top N percentage communications to pick')    
    parser.add_argument( '--metadata_from', type=str, default='metadata/', help='Path to grab the metadata') 
    parser.add_argument( '--output_path', type=str, default='output/', help='Path to save the visualization results, e.g., histograms, graph etc.')
    parser.add_argument( '--top_ccc_file', type=str, default='', help='Path to load the selected top CCC file produced during data postprocessing step')
    parser.add_argument( '--output_name', type=str, default='', help='Output file name prefix according to user\'s choice')
    parser.add_argument( '--N', type=int, default=10000, help='Number of times to run the model for P-value generation')  
    parser.add_argument( '--p_value_cutoff', type=float, default=0.05, help='P value cutoff for filtering the ccc')  
    args = parser.parse_args()

    ######################### read the CellNEST output in csv format ####################################################
    args.metadata_from = args.metadata_from + args.data_name + '/'
    args.data_from = args.data_from + args.data_name + '/'
    args.embedding_path  = args.embedding_path + args.data_name + '/'
    args.output_path = args.output_path + args.data_name + '/'
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

##################### get metadata: input_graph ################################## 

    
    with gzip.open(args.data_from + args.data_name + '_adjacency_records', 'rb') as fp:  #b, a:[0:5]  _filtered 
        row_col, edge_weight, lig_rec, total_num_cell = pickle.load(fp)
    
    lig_rec_db = defaultdict(dict)
    for i in range (0, len(edge_weight)):
        lig_rec_db[lig_rec[i][0]][lig_rec[i][1]] =  edge_weight[i][2]   

    ##########################################
    if args.top_ccc_file == '':
        inFile = args.output_path + args.model_name+'_top' + str(args.top_percent) + 'percent.csv'
        df_org = pd.read_csv(inFile, sep=",")
    else: 
        inFile = args.top_ccc_file
        df_org = pd.read_csv(inFile, sep=",")

    csv_record = df_org
    # columns are: from_cell, to_cell, ligand_gene, receptor_gene, rank, component, from_id, to_id,  attention_score 
    cell_cell_lr_score = defaultdict(dict)
    for record in range (1, len(csv_record)-1):
        i = csv_record[record][6]
        j = csv_record[record][7]
        ligand_gene = csv_record[record][2]
        receptor_gene = csv_record[record][3]
        lr_pair_id = lig_rec_db[ligand_gene][receptor_gene]
        if i in cell_cell_lr_score:
            if j in cell_cell_lr_score[i]: 
                cell_cell_lr_score[i][j][lr_pair_id] = csv_record[record][8]
            else:
                cell_cell_lr_score[i][j] = dict()
                cell_cell_lr_score[i][j][lr_pair_id] = csv_record[record][8]
        else:
            cell_cell_lr_score[i][j] = dict()
            cell_cell_lr_score[i][j][lr_pair_id] = csv_record[record][8]
            
    ################## N time ##########################################
    cell_cell_lr_score_shuffled = defaultdict(dict)
    for shuffle_time in range (0, args.N):
        ## permutate edge feature vector
        edge_weight_temp = copy.deepcopy(edge_weight)
        random.shuffle(edge_weight_temp)
        # reassign the first 2 dimensions = distance and coexpression. Not changing the 3rd dimension (lr pair) because we may not get any value
        # if that is not found between these two cells during edge shuffling
        for i in range (0, len(row_col)):
            edge_weight[i][0] = edge_weight_temp[i][0]
            edge_weight[i][1] = edge_weight_temp[i][1]
        # save edge_weight as a temp
         with gzip.open(args.data_from + args.data_name+'_shuffled' + '_adjacency_records', 'wb') as fp:  #b, a:[0:5]  _filtered 
            pickle.dump([row_col, edge_weight, lig_rec, total_num_cell],fp)
             
        # run the model 3 times, first two asynchronous and 3rd one synchronous
        subprocess.Popen(["bash", "cellnest", "run", "--data_name="+args.data_name+'_shuffled', "--num_epoch=60000", "--model_name="+args.model_name+'_shuffled', "--run_id=1"], stdout="output_"+args.data_name+'_shuffled'+"_run1.log" )
        subprocess.Popen(["bash", "cellnest", "run", "--data_name="+args.data_name+'_shuffled', "--num_epoch=60000", "--model_name="+args.model_name+'_shuffled', "--run_id=2"], stdout="output_"+args.data_name+'_shuffled'+"_run2.log" )
        subprocess.run(["bash", "cellnest", "run", "--data_name="+args.data_name+'_shuffled', "--num_epoch=60000", "--model_name="+args.model_name+'_shuffled', "--run_id=3"], stdout="output_"+args.data_name+'_shuffled'+"_run3.log" )
       
        # postprocess results
        subprocess.run(["bash", "cellnest", "postprocess", "--data_name="+args.data_name+'_shuffled', "--model_name="+args.model_name+'_shuffled', "--total_runs=3"], stdout="output_"+args.data_name+'_shuffled'+"_postproc.log" )

        # read it and get the values
        inFile = args.output_path + args.model_name+'_shuffled'+'_top' + str(args.top_percent) + 'percent_temporary.csv'
        df = pd.read_csv(inFile, sep=",")    
        csv_record = df
        for record in range (1, len(csv_record)-1):
            i = csv_record[record][6]
            j = csv_record[record][7]
            ligand_gene = csv_record[record][2]
            receptor_gene = csv_record[record][3]
            lr_pair_id = lig_rec_db[ligand_gene][receptor_gene]
            if i in cell_cell_lr_score_shuffled:
                if j in cell_cell_lr_score_shuffled[i]: 
                    cell_cell_lr_score_shuffled[i][j][lr_pair_id] = csv_record[record][8]
                else:
                    cell_cell_lr_score_shuffled[i][j] = deafultdict(list)
                    cell_cell_lr_score_shuffled[i][j][lr_pair_id].append(csv_record[record][8])
            else:
                cell_cell_lr_score_shuffled[i][j] = defeaultdict(list)
                cell_cell_lr_score_shuffled[i][j][lr_pair_id].append(csv_record[record][8])
           
    ######################## N times done. Now assign P values ##############################
    # for each i and j cells, for each k lr_pair, find how many times the attention score was 
    # above the original attention score recorded in cell_cell_lr_score
    for i in cell_cell_lr_score:
        for j in cell_cell_lr_score[i]:
            for lr_pair in cell_cell_lr_score[i][j]:
                original_score = cell_cell_lr_score[i][j][lr_pair]
                # how many times higher
                count_higher = 0
                for atn_score in cell_cell_lr_score_shuffled[i][j][lr_pair]:
                    if atn_score > original_score:
                        count_higher = count_higher + 1
                        
                cell_cell_lr_score[i][j][lr_pair] = count_higher/N  # p-value
        
    #########################################################################  
    csv_record = df_org
    # columns are: from_cell, to_cell, ligand_gene, receptor_gene, rank, component, from_id, to_id,  attention_score 
    cell_cell_lr_score = defaultdict(dict)
    csv_record[0].append('p-value')
    for record in range (1, len(csv_record)-1):
        i = csv_record[record][6]
        j = csv_record[record][7]
        ligand_gene = csv_record[record][2]
        receptor_gene = csv_record[record][3]
        lr_pair_id = lig_rec_db[ligand_gene][receptor_gene]
        csv_record[0].append(cell_cell_lr_score[i][j][lr_pair_id])

    csv_record_final = []
    csv_record_final.append(csv_record[0])
    for record in range (1, len(csv_record)-1):
        if csv_record[record][9] <= args.p_value_cutoff:
            csv_record_final.append(csv_record[record])

    
    df = pd.DataFrame(csv_record_final) 
    df.to_csv(args.output_path + args.model_name+'_ccc_pvalue_filtered.csv', index=False, header=False)
