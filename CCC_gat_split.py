from scipy import sparse
import pickle
import numpy as np
from collections import defaultdict
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import DeepGraphInfomax #Linear, 
from torch_geometric.data import Data, DataLoader
#from torch_sparse import SparseTensor
import gzip
import gc
from GATv2Conv_NEST import GATv2Conv
from scipy.sparse import csr_matrix



def get_split_graph(training_data, node_id_sorted, total_subgraphs, expression_matrix_path=''): # use this if you don't want to save the split graph into disk due to space issue
    
    fp = gzip.open(training_data, 'rb')  
    row_col, edge_weight, lig_rec, total_num_cell = pickle.load(fp)
    
    dict_cell_edge = defaultdict(list) # key = node. values = incoming edges
    dict_cell_neighbors = defaultdict(list) # key = node. value = nodes corresponding to incoming edges/neighbors


    nodes_active = dict()
    for i in range(0, len(row_col)): 
        dict_cell_edge[row_col[i][1]].append(i) # index of the edges
        dict_cell_neighbors[row_col[i][1]].append(row_col[i][0]) # neighbor id
        nodes_active[row_col[i][1]] = '' # to 
        nodes_active[row_col[i][0]] = '' # from
    
    
    datapoint_size = len(nodes_active.keys())
    
    for i in range (0, datapoint_size):
        neighbor_list = dict_cell_neighbors[i]
        neighbor_list = list(set(neighbor_list))
        dict_cell_neighbors[i] = neighbor_list
    
    
    fp = gzip.open(node_id_sorted, 'rb')
    node_id_sorted_xy = pickle.load(fp)
    
    node_id_sorted_xy_temp = []
    for i in range(0, len(node_id_sorted_xy)):
        if node_id_sorted_xy[i][0] in nodes_active: # skip those which are not in our ROI
            node_id_sorted_xy_temp.append(node_id_sorted_xy[i])
    
    node_id_sorted_xy = node_id_sorted_xy_temp    
    ##################################################################################################################

    if expression_matrix_path == '':
        # one hot vector used as node feature vector
        #X = np.eye(datapoint_size, datapoint_size)
        #np.random.shuffle(X)
        X = np.arange(datapoint_size)        
        np.random.shuffle(X)
        num_feature = datapoint_size

    else:
        f = gzip.open(expression_matrix_path, 'rb')
        X = pickle.load(f)
        num_feature = X.shape[1]
    
    # split it into N set of edges
    
    total_subgraphs = total_subgraphs
    
    #edge_list = []
    graph_bag = []
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
        spot_list = []
        for k in set1_edges_index:
            i = row_col[k][0]
            j = row_col[k][1]
            if i not in id_map_old_new[set_id]:
                id_map_old_new[set_id][i] = new_id
                id_map_new_old[set_id][new_id] = i
                spot_list.append(new_id)
                new_id = new_id + 1
    
            if j not in id_map_old_new[set_id]:
                id_map_old_new[set_id][j] = new_id
                id_map_new_old[set_id][new_id] = j
                spot_list.append(new_id)
                new_id = new_id + 1
    
    
        #print('new id: %d'%new_id)
        set1_edges = []
        for i in set1_direct_edges:  #set1_edges_index:
            set1_edges.append([[id_map_old_new[set_id][row_col[i][0]], id_map_old_new[set_id][row_col[i][1]]], edge_weight[i]])
            #set1_edges.append([row_col[i], edge_weight[i]])
            
        #edge_list.append(set1_edges)
        
        # create new X matrix
        num_cell = new_id
        X_data = np.zeros((num_cell, datapoint_size))
        spot_id = 0
        for spot in spot_list:
            #X_data[spot_id] = X[spot,:]
            one_column_position = X[spot]
            X_data[spot_id, one_column_position] = 1
            spot_id = spot_id + 1    
        

        row_col_temp = []
        edge_weight_temp = []
        for i in range (0, len(set1_edges)):
            row_col_temp.append(set1_edges[i][0])
            edge_weight_temp.append(set1_edges[i][1])
    
        print("subgraph %d: number of nodes %d, each having feature dimension %d. Total number of edges %d"%(set_id, num_cell, num_feature, len(row_col_temp)))
        
        
        X_data = torch.tensor(X_data, dtype=torch.float)
        X_data = X_data.to_sparse()
        
        edge_index = torch.tensor(np.array(row_col_temp), dtype=torch.long).T
        edge_attr = torch.tensor(np.array(edge_weight_temp), dtype=torch.float)
        
        data = Data(x=X_data, edge_index=edge_index, edge_attr=edge_attr)
        data_loader = DataLoader([data], batch_size=1)
      
        graph_bag.append(data_loader)
        gc.collect()
        
    return graph_bag, num_feature    


def get_graph(training_data): # use this if you already have saved the split graph
    """Add Statement of Purpose
    Args:
        training_data: Path to the input graph    
    Returns:
        List of torch_geometric.data.Data type: Contains the input graph
        Integer: Dimension of node embedding
    """
    
    f = gzip.open(training_data , 'rb')
    graph_bag, num_feature = pickle.load(f)
    
    ###########
    # split it into N set of edges
    ###########

    return graph_bag, num_feature #graph_bags


class Encoder(nn.Module):
    def __init__(self, in_channels, hidden_channels, heads, dropout):
        """Add Statement of Purpose
        Args: [to be]
               
        Returns: [to be]
    
        """

        
        super(Encoder, self).__init__()
        print('incoming channel %d'%in_channels)

        heads = heads
        self.conv =  GATv2Conv(in_channels, hidden_channels, edge_dim=3, heads=heads, concat = False)#, dropout=dropout)
        self.conv_2 =  GATv2Conv(hidden_channels, hidden_channels, edge_dim=3, heads=heads, concat = False)#, dropout=0)

        self.attention_scores_mine_l1 = 'attention_l1'
        self.attention_scores_mine_unnormalized_l1 = 'attention_unnormalized_l1'

        self.attention_scores_mine = 'attention'
        self.attention_scores_mine_unnormalized = 'attention_unnormalized'

        #self.prelu = nn.Tanh(hidden_channels)
        self.prelu = nn.PReLU(hidden_channels)


    def forward(self, data):
        """Add Statement of Purpose
        Args: [to be]
               
        Returns: [to be]
    
        """

        # layer 1
        x, attention_scores, attention_scores_unnormalized = self.conv(data.x, data.edge_index, edge_attr=data.edge_attr, return_attention_weights = True)
        self.attention_scores_mine_l1 = attention_scores
        self.attention_scores_mine_unnormalized_l1 = attention_scores_unnormalized


        # layer 2
        x, attention_scores, attention_scores_unnormalized  = self.conv_2(x, data.edge_index, edge_attr=data.edge_attr, return_attention_weights = True)  # <---- ***
        self.attention_scores_mine = attention_scores #self.attention_scores_mine_l1 #attention_scores
        self.attention_scores_mine_unnormalized = attention_scores_unnormalized #self.attention_scores_mine_unnormalized_l1 #attention_scores_unnormalized

        ###############################
        x = self.prelu(x)

        return x #, attention_scores

class my_data():
    def __init__(self, x, edge_index, edge_attr):
        """Add Statement of Purpose
        Args: [to be]
               
        Returns: [to be]
    
        """
        self.x = x
        self.edge_index = edge_index
        self.edge_attr = edge_attr


def corruption(data):
    """Add Statement of Purpose
    Args: [to be]
           
    Returns: [to be]

    """
    #print('inside corruption function')
    #print(data.x)

    data.x = data.x.to_dense()
    x = data.x[torch.randperm(data.x.size(0))]
    x = x.to_sparse()
    gc.collect()
    return my_data(x, data.edge_index, data.edge_attr)


def train_NEST(args, graph_bag, in_channels):
    """Add Statement of Purpose
    Args: [to be]
           
    Returns: [to be]

    """
    loss_curve = np.zeros((args.num_epoch//500+1))
    loss_curve_counter = 0

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    DGI_model = DeepGraphInfomax(
        hidden_channels=args.hidden,
        encoder=Encoder(in_channels=in_channels, hidden_channels=args.hidden, heads=args.heads, dropout = args.dropout),
        summary=lambda z, *args, **kwargs: torch.sigmoid(z.mean(dim=0)),
        corruption=corruption).to(device)

    DGI_optimizer = torch.optim.Adam(DGI_model.parameters(), lr=args.lr_rate) #1e-5)#5 #6
    DGI_filename = args.model_path+'DGI_'+ args.model_name  +'.pth.tar'

    if args.load:
        DGI_load_path = args.model_path+'DGI_'+ args.load_model_name+'.pth.tar'
        DGI_model.load_state_dict(torch.load(DGI_load_path))
        DGI_optimizer.load_state_dict(torch.load(args.model_path+'DGI_optimizer_'+ args.load_model_name  +'.pth.tar'))

    import datetime
    start_time = datetime.datetime.now()
    min_loss=10000
    print('Saving init model state ...')
    torch.save(DGI_model.state_dict(), args.model_path+'DGI_init'+ args.model_name  + '.pth.tar')
    torch.save(DGI_optimizer.state_dict(), args.model_path+'DGI_optimizer_init'+ args.model_name  + '.pth.tar')
    print('training starts ...')
    for epoch in range(args.num_epoch):
        DGI_model.train()
        DGI_optimizer.zero_grad()
        DGI_all_loss = []

        for subgraph in graph_bag:
            for data in subgraph:
                data = data.to(device)
                pos_z, neg_z, summary = DGI_model(data=data)
                DGI_loss = DGI_model.loss(pos_z, neg_z, summary)
                DGI_loss.backward()
                DGI_all_loss.append(DGI_loss.item())
                
        DGI_optimizer.step()

        if ((epoch)%500) == 0:
            print('Epoch: {:03d}, Loss: {:.4f}'.format(epoch+1, np.mean(DGI_all_loss)))
            loss_curve[loss_curve_counter] = np.mean(DGI_all_loss)
            loss_curve_counter = loss_curve_counter + 1

            if np.mean(DGI_all_loss)<min_loss:

                min_loss=np.mean(DGI_all_loss)

                # save the current model state
                torch.save(DGI_model.state_dict(), DGI_filename)
                torch.save(DGI_optimizer.state_dict(), args.model_path+'DGI_optimizer_'+ args.model_name  +'.pth.tar')
                save_tupple=[pos_z, neg_z, summary] 
                ############################################################################################################
                subgraph_id = -1
                for subgraph in graph_bag:
                    subgraph_id = subgraph_id + 1
                    for data in subgraph:
                        data = data.to(device)
                        pos_z, neg_z, summary = DGI_model(data=data)              
               
                    # save the node embedding
                    X_embedding = pos_z
                    X_embedding = X_embedding.cpu().detach().numpy()
                    X_embedding_filename =  args.embedding_path + args.model_name + '_Embed_X' #.npy
                    with gzip.open(X_embedding_filename+'_subgraph'+str(subgraph_id), 'wb') as fp:  
                        pickle.dump(X_embedding, fp)
                        

                    # save the attention scores
    
                    X_attention_index = DGI_model.encoder.attention_scores_mine[0]
                    X_attention_index = X_attention_index.cpu().detach().numpy()
    
                    # layer 1
                    X_attention_score_normalized_l1 = DGI_model.encoder.attention_scores_mine_l1[1]
                    X_attention_score_normalized_l1 = X_attention_score_normalized_l1.cpu().detach().numpy()
                    # layer 1 unnormalized
                    X_attention_score_unnormalized_l1 = DGI_model.encoder.attention_scores_mine_unnormalized_l1
                    X_attention_score_unnormalized_l1 = X_attention_score_unnormalized_l1.cpu().detach().numpy()
    
                    # layer 2
                    X_attention_score_normalized = DGI_model.encoder.attention_scores_mine[1]
                    X_attention_score_normalized = X_attention_score_normalized.cpu().detach().numpy()
                    # layer 2 unnormalized
                    X_attention_score_unnormalized = DGI_model.encoder.attention_scores_mine_unnormalized
                    X_attention_score_unnormalized = X_attention_score_unnormalized.cpu().detach().numpy()
    
                    print('making the bundle to save')
                    X_attention_bundle = [X_attention_index, X_attention_score_normalized_l1, X_attention_score_unnormalized, X_attention_score_unnormalized_l1, X_attention_score_normalized]
                    X_attention_filename =  args.embedding_path + args.model_name + '_attention'+'_subgraph'+str(subgraph_id)
                    with gzip.open(X_attention_filename, 'wb') as fp:  
                        pickle.dump(X_attention_bundle, fp)

                print('all saved')
                ############################################################################################################################
                logfile=open(args.model_path+'DGI_'+ args.model_name+'_loss_curve.csv', 'wb')
                np.savetxt(logfile,loss_curve, delimiter=',')
                logfile.close()

                #print(DGI_model.encoder.attention_scores_mine_unnormalized_l1[0:10])

#        if ((epoch)%40000) == 0:
#            DGI_optimizer = torch.optim.Adam(DGI_model.parameters(), lr=0.00001)  #5 #6

    end_time = datetime.datetime.now()

#        torch.save(DGI_model.state_dict(), DGI_filename)
    print('Training time in seconds: ', (end_time-start_time).seconds)
    DGI_model.load_state_dict(torch.load(DGI_filename))
    print("debug loss")
    DGI_loss = DGI_model.loss(pos_z, neg_z, summary)
    print("debug loss latest tupple %g"%DGI_loss.item())
    DGI_loss = DGI_model.loss(save_tupple[0], save_tupple[1], save_tupple[2])
    print("debug loss min loss tupple %g"%DGI_loss.item())

    return DGI_model

