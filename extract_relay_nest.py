# Written By 
# Fatema Tuz Zohora

import numpy as np
import csv
import pickle
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.colors import  rgb2hex # LinearSegmentedColormap, to_hex,
from scipy.sparse import csr_matrix
from collections import defaultdict
import pandas as pd
import gzip
import argparse
import os
import scipy.stats
from scipy.sparse.csgraph import connected_components
from pyvis.network import Network
import networkx as nx
from networkx.drawing.nx_agraph import write_dot
import altair as alt
import altairThemes # assuming you have altairThemes.py at your current directoy or your system knows the path of this altairThemes.py.
import gc
import copy
alt.themes.register("publishTheme", altairThemes.publishTheme)
# enable the newly registered theme
alt.themes.enable("publishTheme")


#current_directory = ??

##########################################################
# preprocessDf, plot: these two functions are taken from GW's repository                                                                                                                                                                     /mnt/data0/gw/research/notta_pancreatic_cancer_visium/plots/fatema_signaling/hist.py                                                                                                                                                                                         

def preprocessDf(df):
  """Transform ligand and receptor columns."""
  df["ligand-receptor"] = df["ligand"] + '-' + df["receptor"]
  df["component"] = df["component"] #.astype(str).str.zfill(2)

  return df


def plot(df):
  set1 = altairThemes.get_colour_scheme("Set1", len(df["component"].unique()))
  set1[0] = '#000000'
  base = alt.Chart(df).mark_bar().encode(
            x=alt.X("ligand-receptor:N", axis=alt.Axis(labelAngle=45), sort='-y'),
            y=alt.Y("count()"),
            color=alt.Color("component:N", scale = alt.Scale(range=set1)),
            order=alt.Order("component:N", sort="ascending"),
            tooltip=["component"]
        )
  p = base

  return p

####################### Set the name of the sample you want to visualize ###################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( '--data_name', type=str, required = True, help='The name of dataset') # 
    parser.add_argument( '--metadata', type=str, default='metadata/', help='The name of dataset') #
    parser.add_argument( '--barcode_info_file', type=str, default='', help='Path to load the barcode information file produced during data preprocessing step')
    parser.add_argument( '--annotation_file_path', type=str, default='', help='Path to load the annotation file in csv format (if available) ') #_ayah_histology
    parser.add_argument( '--selfloop_info_file', type=str, default='', help='Path to load the selfloop information file produced during data preprocessing step')
    parser.add_argument( '--top_ccc_file', type=str, required = True, help='Path to load the selected top CCC file produced during data postprocessing step')
    parser.add_argument( '--output_path', type=str, default='NEST_figures_output/', help='Output file name prefix according to user\'s choice')
    parser.add_argument( '--target_relay', type=str, default='', help='Target relay to plot according to user\'s choice. For example: ')
    
    args = parser.parse_args()




    output_path = args.output_path
    
    ##################### make cell metadata: barcode_info ###################################
    if args.barcode_info_file == '':
        args.barcode_info_file = args.metadata + args.data_name + '/' + args.data_name + '_barcode_info'
        
    with gzip.open(args.barcode_info_file, 'rb') as fp:     
        barcode_info = pickle.load(fp)    

    ###############################  read which spots have self loops ###############################################################
    if args.selfloop_info_file == '':
    	args.selfloop_info_file = args.metadata + args.data_name + '/' + args.data_name + '_self_loop_record'
        
    with gzip.open(args.selfloop_info_file, 'rb') as fp: 
        self_loop_found = pickle.load(fp)

    ####### load annotations ##############################################
    if args.annotation_file_path != '':
        annotation_data = pd.read_csv(args.annotation_file_path, sep=",")
        pathologist_label=[]
        for i in range (0, len(annotation_data)):
            pathologist_label.append([annotation_data['Barcode'][i], annotation_data['IX_annotation'][i]])
    
        barcode_type=dict() # record the type (annotation) of each spot (barcode)
        for i in range (0, len(pathologist_label)):
            barcode_type[pathologist_label[i][0]] = pathologist_label[i][1]

    else:
        barcode_type=dict() # record the type (annotation) of each spot (barcode)
        for i in range (0, len(barcode_info)):
            barcode_type[barcode_info[i][0]] = ''
      
    #################################################################################################################
    '''
    csv_record = df.values.tolist() # barcode_info[i][0], barcode_info[j][0], ligand, receptor, edge_rank, label, i, j, score
    ## sort the edges based on their rank (column 4), low to high, low being higher attention score
    # columns are: from_cell, to_cell, ligand_gene, receptor_gene, rank, attention_score, component, from_id, to_id
    df_column_names = list(df.columns)
#    print(df_column_names)

    print(len(csv_record))

    if args.top_edge_count != -1:
        csv_record_final = [df_column_names] + csv_record [0:min(args.top_edge_count, len(csv_record))]
    '''
    #####################################
    # read csv_record_final
    inFile = args.top_ccc_file
    df = pd.read_csv(inFile, sep=",")
    csv_record_final = [list(df.columns)] + df.values.tolist()

    #####################################
    component_list = dict()
    for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
        record = csv_record_final[record_idx]
        i = record[6]
        j = record[7]
        component_label = record[5]
        barcode_info[i][3] = component_label #?
        barcode_info[j][3] = component_label #?
        component_list[component_label] = ''

    component_list[0] = ''
    unique_component_count = max(max(component_list.keys()), 2)	+ 1 # at least 2, 0 = no CCC, 1 = with CCC.
    # max(component_list.keys()) = because it will keep the colors matched with the original component plot 

    ####################### pattern finding ##########################################################################
    # make a dictionary to keep record of all the outgoing edges [to_node, ligand, receptor] for each node
    each_node_outgoing = defaultdict(list)
    for k in range (1, len(csv_record_final)-1): # last record is a dummy for histogram preparation
        i = csv_record_final[k][6]
        j = csv_record_final[k][7]
        if i == j:
            continue        
        ligand = csv_record_final[k][2]
        receptor = csv_record_final[k][3]
        each_node_outgoing[i].append([j, ligand, receptor, k]) 
    
    # all possible 2-hop pattern count
    pattern_distribution = defaultdict(list)
    pattern_distribution_cell_info = defaultdict(list)
    # pattern_distribution['ligand-receptor to ligand-receptor']=[1,1,1,1, ...]
    edge_list_2hop = []
#    target_relay = 'APOE-SDC4 to APOE-SDC1' #'PSAP-LRP1 to APOE-LRP1' #'PSAP-LRP1 to PSAP-LRP1'
    for i in each_node_outgoing:
        for tupple in each_node_outgoing[i]: # first hop
            j = tupple[0]
            lig_rec_1 = tupple[1]+'-'+tupple[2]
            record_id_1 = tupple[3]
            if j in each_node_outgoing:
                for tupple_next in each_node_outgoing[j]: # second hop
                    k = tupple_next[0]
                    if k == i or k == j:
                        continue
                    lig_rec_2 = tupple_next[1]+'-'+tupple_next[2]
                    record_id_2 = tupple_next[3]
                    pattern_distribution[lig_rec_1 + ' to ' + lig_rec_2].append(1)
                    pattern_distribution_cell_info[lig_rec_1 + ' to ' + lig_rec_2].append([barcode_info[i], barcode_info[j], barcode_info[k]]) # a_id, b_id, c_id
                    relay = lig_rec_1 + ' to ' + lig_rec_2
                    if args.target_relay != "" and relay == args.target_relay:
                        edge_list_2hop.append([record_id_1,record_id_2])
    
    

    two_hop_pattern_distribution = []
    same_count = 0
    for key in pattern_distribution:
        count = len(pattern_distribution[key])
        two_hop_pattern_distribution.append([key, count]) 
        #if lig_rec_1 == lig_rec_2:
        #    same_count = same_count + 1
    
    two_hop_pattern_distribution = sorted(two_hop_pattern_distribution, key = lambda x: x[1], reverse=True) # high to low
    
    data_list=dict()
    data_list['X']=[]
    data_list['Y']=[] 
    for i in range (0, 100): #len(two_hop_pattern_distribution)):
        data_list['X'].append(two_hop_pattern_distribution[i][0])
        data_list['Y'].append(two_hop_pattern_distribution[i][1])
        
    data_list_pd = pd.DataFrame({
        'Relay Patterns': data_list['X'],
        'Pattern Abundance (#)': data_list['Y']
    })

    chart = alt.Chart(data_list_pd).mark_bar().encode(
        x=alt.X("Relay Patterns:N", axis=alt.Axis(labelAngle=45), sort='-y'),
        y='Pattern Abundance (#)'
    )
    chart.save(output_path + 'CellNEST_' + args.data_name +'_relay_pattern_histograms.html')
    ######################### save as table format and numpy format ########
    #data_list_pd.to_csv(output_path + args.data_name +'_top20p_topEdge'+str(args.top_edge_count)+'_relay_count.csv', index=False)
    data_list_pd.to_csv(output_path + 'CellNEST_' + args.data_name +'_relay_pattern_count.csv', index=False)
    with gzip.open(output_path + 'CellNEST_' + args.data_name + '_relay_pattern_cell_info', 'wb') as fp:
        # This will be used during cell type finding for relays
        pickle.dump(pattern_distribution_cell_info, fp)

    print('Relay pattern extraction done and saved at:\n'+
          output_path + 'CellNEST_' + args.data_name +'_relay_pattern_histograms.html' + '\n' +
          output_path + 'CellNEST_' + args.data_name +'_relay_pattern_count.csv' + '\n' +
          output_path + 'CellNEST_' + args.data_name + '_relay_pattern_cell_info'
          )

    '''
    with gzip.open(output_path + 'CellNEST_' + args.data_name + '_relay_pattern_cell_info', 'rb') as fp:
            pattern_distribution_cell_info = pickle.load(fp)

    pattern_list = list(pattern_distribution_cell_info.keys()) 
    # see which cells are forming the first pattern
    group_list = pattern_distribution_cell_info[pattern_list[0]]
    for group_id in range(0, len(group_list)):
        print('group: %d'%group_id)
        print(group_list[group_id])
    '''
    ######################### Plotting the two hops #####################
    set1 = altairThemes.get_colour_scheme("Set1", unique_component_count)
    colors = set1
    colors[0] = '#000000' # black means no CCC
    ids = []
    x_index=[]
    y_index=[]
    colors_point = []
    for i in range (0, len(barcode_info)):    
        ids.append(i)
        x_index.append(barcode_info[i][1])
        y_index.append(barcode_info[i][2])    
        colors_point.append(colors[barcode_info[i][3]]) 
    
    max_x = np.max(x_index)
    max_y = np.max(y_index)



    g = nx.MultiDiGraph(directed=True) 
    for i in range (0, len(barcode_info)):        
        marker_size = 'circle'
        label_str =  str(i)+'_c:'+str(barcode_info[i][3]) #  label of the node or spot is consists of: spot id, component number
        if args.annotation_file_path != '': 
            if 'tumor' in str(barcode_type[barcode_info[i][0]]): # Tumour
                marker_size = 'box'
                label_str = label_str +'_'+ 'tumor' # also add the type of the spot to the label if annotation is available 
        
        g.add_node(int(ids[i]), x=int(x_index[i]), y=int(y_index[i]), label = label_str, pos = str(x_index[i])+","+str(-y_index[i])+" !", physics=False, shape = marker_size, color=matplotlib.colors.rgb2hex(colors_point[i]))    


    # scale the edge scores [0 to 1] to make width in the plot work
    score_list = []
    for k in range (1, len(csv_record_final)-1):
        score_list.append(csv_record_final[k][8])

    min_score = np.min(score_list)
    max_score = np.max(score_list)

    count_edges = 0
    if len(edge_list_2hop) > 0:
        for relay in edge_list_2hop:
            print(relay)
            for hop in range (0, len(relay)):
                k = relay[hop] # record id for the hops
                i = csv_record_final[k][6]
                j = csv_record_final[k][7] 
                print('%d %d'%(i, j))
                ligand = csv_record_final[k][2]
                receptor = csv_record_final[k][3]            
                edge_score = csv_record_final[k][8]
                edge_score = (edge_score-min_score)/(max_score-min_score)   
                title_str =  "L:" + ligand + ", R:" + receptor+ ", "+ str(edge_score) #+
                g.add_edge(int(i), int(j), label = title_str, color=colors_point[i], value=np.float64(edge_score)) #
                count_edges = count_edges + 1
    
        print("total edges plotted: %d"%count_edges)
    
        nt = Network( directed=True, height='1000px', width='100%') #"500px", "500px",, filter_menu=True     
        nt.from_nx(g)
        nt.save_graph(output_path + 'CellNEST_' +args.data_name +'_relay_mygraph.html')
        print('Edge graph plot generation done and saved at: '+output_path + 'CellNEST_' +args.data_name +'_relay_mygraph.html')
        ########################################################################
        # convert it to dot file to be able to convert it to pdf or svg format for inserting into the paper
        write_dot(g, output_path + 'CellNEST_' + args.data_name + "_relay_test_interactive.dot")
        print('dot file generation done and saved at:'+output_path + 'CellNEST_' + args.data_name + "_relay_test_interactive.dot")

