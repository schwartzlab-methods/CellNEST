# Written By 
# Fatema Tuz Zohora

print('package loading')
import numpy as np
import csv
import pickle
import matplotlib
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

    parser.add_argument( '--data_name', type=str, help='The name of dataset', default="Visium_HD_Human_Colon_Cancer_square_002um_outputs") # , required=True
    parser.add_argument( '--model_name', type=str, help='Name of the trained model', default='NEST_Visium_HD_Human_Colon_Cancer_square_002um_outputs') #, required=True
    parser.add_argument( '--top_edge_count', type=int, default=135000 ,help='Number of the top communications to plot. To plot all insert -1') # 
    parser.add_argument( '--top_percent', type=int, default=20, help='Top N percentage communications to pick')    
    parser.add_argument( '--metadata_from', type=str, default='metadata/', help='Path to grab the metadata') 
    parser.add_argument( '--output_path', type=str, default='output/', help='Path to save the visualization results, e.g., histograms, graph etc.')
    parser.add_argument( '--barcode_info_file', type=str, default='', help='Path to load the barcode information file produced during data preprocessing step')
    parser.add_argument( '--annotation_file_path', type=str, default='', help='Path to load the annotation file in csv format (if available) ')
    parser.add_argument( '--selfloop_info_file', type=str, default='', help='Path to load the selfloop information file produced during data preprocessing step')
    parser.add_argument( '--top_ccc_file', type=str, default='', help='Path to load the selected top CCC file produced during data postprocessing step')
    parser.add_argument( '--output_name', type=str, default='', help='Output file name prefix according to user\'s choice')
    parser.add_argument( '--filter', type=int, default=0, help='Set --filter=-1 if you want to filter the CCC')
    parser.add_argument( '--filter_by_ligand_receptor', type=str, default='', help='Set ligand-receptor pair, e.g., --filter_by_ligand_receptor="CCL19-CCR7" if you want to filter the CCC by LR pair')
    parser.add_argument( '--filter_by_annotation', type=str, default='', help='Set cell or spot type, e.g., --filter_by_annotation="T-cell" if you want to filter the CCC')
    parser.add_argument( '--filter_by_component', type=int, default=-1, help='Set component id, e.g., --filter_by_component=9 if you want to filter by component id')
    
    
    args = parser.parse_args()
    if args.metadata_from=='metadata/': # if default one is used, then concatenate the dataname. Otherwise, use the user provided path directly
        args.metadata_from = args.metadata_from + args.data_name + '/'

    if args.output_path=='output/': # if default one is used, then concatenate the dataname. Otherwise, use the user provided path directly
        args.output_path = args.output_path + args.data_name + '/'
    print('Top %d communications will be plot. To change the count use --top_edge_count parameter'%args.top_edge_count)

    if args.output_name=='':
        output_name = args.output_path + args.model_name
    else: 
        output_name = args.output_path + args.output_name
    
    ##################### make cell metadata: barcode_info ###################################
    if args.barcode_info_file=='':
        with gzip.open(args.metadata_from +args.data_name+'_barcode_info', 'rb') as fp:  #b, a:[0:5]   
            barcode_info = pickle.load(fp)
    else:
        with gzip.open(args.barcode_info_file, 'rb') as fp:  #b, a:[0:5]        
            barcode_info = pickle.load(fp)    


    ###############################  read which spots have self loops ################################################################
    if args.selfloop_info_file=='':
        with gzip.open(args.metadata_from + args.data_name +'_self_loop_record', 'rb') as fp:  #b, a:[0:5]   _filtered
            self_loop_found = pickle.load(fp)
    else:
        with gzip.open(args.selfloop_info_file, 'rb') as fp:  #b, a:[0:5]   _filtered
            self_loop_found = pickle.load(fp)

    ####### load annotations ##############################################
    if args.annotation_file_path != '':
        pathologist_label=[]
        annotation_data = pd.read_csv(args.annotation_file_path, sep=",")
        for i in range (0, len(annotation_data)):
            pathologist_label.append([annotation_data['Barcode'][i], annotation_data['Type'][i]])

        barcode_type=dict() # record the type (annotation) of each spot (barcode)
        for i in range (0, len(pathologist_label)):
            barcode_type[pathologist_label[i][0]] = pathologist_label[i][1]


    else:
        barcode_type=dict() # record the type (annotation) of each spot (barcode)
        for i in range (0, len(barcode_info)):
            barcode_type[barcode_info[i][0]] = ''

    ######################### read the NEST output in csv format ####################################################
    if args.top_ccc_file == '':
        inFile = args.output_path + args.model_name+'_top' + str(args.top_percent) + 'percent.csv'
        df = pd.read_csv(inFile, sep=",")
    else: 
        inFile = args.top_ccc_file
        df = pd.read_csv(inFile, sep=",")


    csv_record = df.values.tolist() # barcode_info[i][0], barcode_info[j][0], ligand, receptor, edge_rank, label, i, j, score

    ## sort the edges based on their rank (column 4), low to high, low being higher attention score
    csv_record = sorted(csv_record, key = lambda x: x[4])
    ## add the column names and take first top_edge_count edges
    # columns are: from_cell, to_cell, ligand_gene, receptor_gene, rank, component, from_id, to_id,  attention_score 
    df_column_names = list(df.columns)
#    print(df_column_names)

    print("Total edges in the topN%% file is: %d"%len(csv_record))

    if args.top_edge_count != -1:
        csv_record_final = [df_column_names] + csv_record[0:min(args.top_edge_count, len(csv_record))]

    ## add a dummy row at the end for the convenience of histogram preparation (to keep the color same as altair plot)
    '''
    in_region_node = -1
    for i in range (0, len(barcode_info)):
        if barcode_info[i][1] <= 54000 :
            in_region_node = i
            break
  
    i = in_region_node
    j = in_region_node
    csv_record_final.append([barcode_info[i][0], barcode_info[j][0], 'no-ligand', 'no-receptor', 0, 0, i, j, 0]) # dummy for histogram
    '''

    csv_record_final.append([csv_record_final[1][0], csv_record_final[1][0], 'no-ligand', 'no-receptor', 0, 0, csv_record_final[1][6], csv_record_final[1][6], 0]) # dummy for histogram
    csv_record = 0
    gc.collect()

    ######################## connected component finding #################################
    print('Finding connected component')
    connecting_edges = np.zeros((len(barcode_info),len(barcode_info)))  
    for k in range (1, len(csv_record_final)-1): # last record is a dummy for histogram preparation
        i = csv_record_final[k][6]
        j = csv_record_final[k][7]
        connecting_edges[i][j]=1
            
    graph = csr_matrix(connecting_edges)
    n_components, labels = connected_components(csgraph=graph,directed=True, connection = 'weak',  return_labels=True) # It assigns each SPOT to a component based on what pair it belongs to
#    print('Number of connected components %d'%n_components) 

    count_points_component = np.zeros((n_components))
    for i in range (0, len(labels)):
        count_points_component[labels[i]] = count_points_component[labels[i]] + 1

    id_label = 2 # initially all are zero. =1 those who have self edge but above threshold. >= 2 who belong to some component
    index_dict = dict()
    for i in range (0, count_points_component.shape[0]):
        if count_points_component[i]>1:
            index_dict[i] = id_label
            id_label = id_label+1

    print('component count %d'%id_label)

    for i in range (0, len(barcode_info)):
        if count_points_component[labels[i]] > 1:
            barcode_info[i][3] = index_dict[labels[i]] #2
        elif connecting_edges[i][i] == 1 and (i in self_loop_found and i in self_loop_found[i]): # that is: self_loop_found[i][i] do exist 
            barcode_info[i][3] = 1
        else: 
            barcode_info[i][3] = 0

    # update the label based on found component numbers
    #max opacity
    for record in range (1, len(csv_record_final)-1):
        i = csv_record_final[record][6]
        label = barcode_info[i][3]
        csv_record_final[record][5] = label
    
    ############################################### Optional filtering ########################################################
    if args.filter == 1:
        ## change the csv_record_final here if you want histogram for specific components/regions only. e.g., if you want to plot only stroma region, or tumor-stroma regions etc.    ##
        #region_of_interest = [...] 
        csv_record_final_temp = []
        csv_record_final_temp.append(csv_record_final[0])
        component_dictionary_dummy = dict()
        for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
            if args.filter_by_component!=-1:
                if csv_record_final[record_idx][5] == int(args.filter_by_component):
                    csv_record_final_temp.append(csv_record_final[record_idx])                
            elif args.filter_by_annotation!='': 
                if barcode_type[csv_record_final[record_idx][0]] == args.filter_by_annotation and barcode_type[csv_record_final[record_idx][1]] == args.filter_by_annotation: # if from_node == type and to_node == type
                    csv_record_final_temp.append(csv_record_final[record_idx])   
                if csv_record_final[record_idx][5] not in component_dictionary_dummy:
                    component_dictionary_dummy[csv_record_final[record_idx][5]] = csv_record_final[record_idx]                
            elif args.filter_by_ligand_receptor!='':
                ligand = (args.filter_by_ligand_receptor).split('-')[0]
                receptor = (args.filter_by_ligand_receptor).split('-')[1]
                if csv_record_final[record_idx][2] == ligand and csv_record_final[record_idx][3] == receptor:
                    csv_record_final_temp.append(csv_record_final[record_idx])                  
            else: # if no 'filter by' options are provided by mistake 
                csv_record_final_temp.append(csv_record_final[record_idx])
            
    
        # insert just one record from each other components so that the color scheme does not change in the altair scatter plot and histogram :-(
        for component_id in component_dictionary_dummy:
            csv_record_final_temp.append(component_dictionary_dummy[component_id])
        
        csv_record_final_temp.append(csv_record_final[len(csv_record_final)-1])
        csv_record_final = copy.deepcopy(csv_record_final_temp)
                    

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
    unique_component_count = max(len(component_list.keys()), id_label)

    ####################### pattern finding ##########################################################################
    '''
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
    # pattern_distribution['ligand-receptor to ligand-receptor']=[1,1,1,1, ...]
    edge_list_2hop = []
    target_relay = 'C3-CXCR4 to C3-LRP1'
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
                    relay = lig_rec_1 + ' to ' + lig_rec_2
                    if relay == target_relay:
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
    for i in range (0, 15): #len(two_hop_pattern_distribution)):
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

    chart.save(output_name + '_pattern_distribution.html')
    '''

    ##################################### Altair Plot ##################################################################
    ## dictionary of those spots who are participating in CCC ##
    active_spot = defaultdict(list)
    for record_idx in range (1, len(csv_record_final)-1): #last entry is a dummy for histograms, so ignore it.
        record = csv_record_final[record_idx]
        i = record[6]
        pathology_label = barcode_type[barcode_info[i][0]]
        component_label = record[5]
        X = barcode_info[i][1]
        Y = -barcode_info[i][2]
        opacity = float(record[8])
        active_spot[i].append([pathology_label, component_label, X, Y, opacity])
        
        j = record[7]
        pathology_label = barcode_type[barcode_info[j][0]]
        component_label = record[5]
        X = barcode_info[j][1]
        Y = -barcode_info[j][2]
        opacity = float(record[8])   
        active_spot[j].append([pathology_label, component_label, X, Y, opacity])
        ''''''
        
    ######### color the spots in the plot with opacity = attention score #################
    opacity_list = []
    for i in active_spot:
        sum_opacity = []
        for edges in active_spot[i]:
            sum_opacity.append(edges[4])
            
        avg_opacity = np.max(sum_opacity) #np.mean(sum_opacity)
        opacity_list.append(avg_opacity) 
        active_spot[i]=[active_spot[i][0][0], active_spot[i][0][1], active_spot[i][0][2], active_spot[i][0][3], avg_opacity]

    min_opacity = np.min(opacity_list)
    max_opacity = np.max(opacity_list)

    #### making dictionary for converting to pandas dataframe to draw altair plot ###########
    data_list=dict()
    data_list['pathology_label']=[]
    data_list['component_label']=[]
    data_list['X']=[]
    data_list['Y']=[]   
    data_list['opacity']=[]  

    for i in range (0, len(barcode_info)):  
        if barcode_info[i][1] > 54000:
            continue
        if i in active_spot:
            data_list['pathology_label'].append(active_spot[i][0])
            data_list['component_label'].append(active_spot[i][1])
            data_list['X'].append(active_spot[i][2])
            data_list['Y'].append(active_spot[i][3])
            data_list['opacity'].append((active_spot[i][4]-min_opacity)/(max_opacity-min_opacity))
            
        else:
            data_list['pathology_label'].append(barcode_type[barcode_info[i][0]])
            data_list['component_label'].append(0) # make it zero so it is black
            data_list['X'].append(barcode_info[i][1])
            data_list['Y'].append(-barcode_info[i][2])
            data_list['opacity'].append(0.1)
            # barcode_info[i][3] = 0



    # converting to pandas dataframe

    data_list_pd = pd.DataFrame(data_list)
    id_label = len(list(set(data_list['component_label']))) # 
    unique_component_count = id_label
    set1 = altairThemes.get_colour_scheme("Set1", id_label)
    set1[0] = '#000000'
    chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
        alt.X('X', scale=alt.Scale(zero=False)),
        alt.Y('Y', scale=alt.Scale(zero=False)),
        shape = alt.Shape('pathology_label:N'), #shape = "pathology_label",
        color=alt.Color('component_label:N', scale=alt.Scale(range=set1)),
        #opacity=alt.Opacity('opacity:N'), #"opacity", 
        tooltip=['component_label'] #,'opacity'
    )

    chart.save(output_name +'_component_plot.html')
    print('Altair plot generation done')

    ###################################  Histogram plotting #################################################################################

    df = pd.DataFrame(csv_record_final)
    df.to_csv('temp_csv.csv', index=False, header=False)
    df = pd.read_csv('temp_csv.csv', sep=",")
    os.remove('temp_csv.csv') # delete the intermediate file

    print('len of loaded csv for histogram generation is %d'%len(df))
    df = preprocessDf(df)
    p = plot(df)
    outPath = output_name +'_histogram_test.html'
    p.save(outPath)	
    print('Histogram plot generation done')

    ############################  Network/edge graph plot ######################

    set1 = altairThemes.get_colour_scheme("Set1", unique_component_count)
    colors = set1
    colors[0] = '#000000' # black means no CCC
    ids = []
    x_index=[]
    y_index=[]
    colors_point = []
    barcode_to_serial = dict()
    for i in range (0, len(barcode_info)): 
        if barcode_info[i][1] > 54000:
            continue
        ids.append(i)
        x_index.append(barcode_info[i][1])
        y_index.append(barcode_info[i][2])    
        barcode_to_serial[i] = len(colors_point)
        colors_point.append(colors[barcode_info[i][3]]) 
    
    max_x = np.max(x_index)
    max_y = np.max(y_index)



    g = nx.MultiDiGraph(directed=True) 
    list_index = 0
    for i in ids:
        marker_size = 'circle'
        label_str =  str(i)+'_c:'+str(barcode_info[i][3]) #  label of the node or spot is consists of: spot id, component number
        if args.annotation_file_path != '':
            label_str = label_str +'_'+ barcode_type[barcode_info[i][0]] # also add the type of the spot to the label if annotation is available 
        
        g.add_node(int(i), x=int(x_index[list_index]), y=int(y_index[list_index]), label = label_str, pos = str(x_index[list_index])+","+str(-y_index[list_index])+" !", physics=False, shape = marker_size, color="#808080") #=matplotlib.colors.rgb2hex(colors_point[list_index])    
        list_index = list_index + 1


    # scale the edge scores [0 to 1] to make plot work
    score_list = []
    for k in range (1, len(csv_record_final)-1):
        score_list.append(csv_record_final[k][8])

    min_score = np.min(score_list)
    max_score = np.max(score_list)
    ###################################################################
    count_edges = 0
    for k in range (1, len(csv_record_final)-1):
        i = csv_record_final[k][6]
        j = csv_record_final[k][7] 

        ligand = csv_record_final[k][2]
        receptor = csv_record_final[k][3]

        #if ligand=='CCL19' and receptor=='CCR7':
        #    print('CCL19-CCR7')

        edge_score = csv_record_final[k][8]
        edge_score = (edge_score-min_score)/(max_score-min_score)   
        title_str =  "L:" + ligand + ", R:" + receptor+ ", "+ str(edge_score) #+
        g.add_edge(int(i), int(j), label = title_str, color=colors_point[barcode_to_serial[i]], value=np.float64(edge_score)) #
        count_edges = count_edges + 1

    print("total edges plotted: %d"%count_edges)
    nt = Network( directed=True, height='1000px', width='100%') #"500px", "500px",, filter_menu=True     
    nt.from_nx(g)
    nt.save_graph(output_name +'_mygraph.html')
    print('Edge graph plot generation done')
    ########################################################################
    # convert it to dot file to be able to convert it to pdf or svg format for inserting into the paper
    write_dot(g, output_name + "_test_interactive.dot")
    print('dot file generation done')
    print('All done')
    ############################ or Relay plot ########################################

    count_edges = 0
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
            g.add_edge(int(i), int(j), label = title_str, color="#FF0000", value=np.float64(edge_score)) #
            count_edges = count_edges + 1

    print("total edges plotted: %d"%count_edges)

  
    nt = Network( directed=True, height='1000px', width='100%') #"500px", "500px",, filter_menu=True     
    nt.from_nx(g)
    nt.save_graph(output_name +'_relay_mygraph.html')
    print('Edge graph plot generation done')
    ########################################################################
    # convert it to dot file to be able to convert it to pdf or svg format for inserting into the paper
    write_dot(g, output_name + "_relay_test_interactive.dot")
    print('dot file generation done')
    print('All done')
