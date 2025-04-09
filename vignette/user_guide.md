We provide the arguments for running various steps of CellCellNEST model along with the preferred values as Default. Most of the arguments are self-explanatory. The parameters that may need changing are discussed. 

## Data Preprocess 
Sample command:
```
cellnest preprocess --data_name='V1_Human_Lymph_Node_spatial' --data_from='data/V1_Human_Lymph_Node_spatial/'
```
### Arguments
    --data_name = Name of the dataset. Type is String. Required
    --data_from = Path to the dataset to read from. Space Ranger outs/ folder is preferred. Otherwise, provide the *.mtx file of the gene expression matrix. Type is String. Required.
    --data_to = Path to save the input graph (to be passed to GAT). Type is String. Default = 'input_graph/'.
    --metadata_to = Path to save the metadata. Type is String. Default='metadata/'
    --filter_min_cell = Minimum number of cells for gene filtering. Type is Int. Default=1. 
    --threshold_gene_exp = Threshold percentile for gene expression. Genes above this percentile are considered active. Type is float. Default=98
    --tissue_position_file = If your --data_from argument points to a *.mtx file instead of Space Ranger, then please provide the path to tissue position file. Type is String. Default='None'
    --spot_diameter = Spot/cell diameter for filtering ligand-receptor pairs based on cell-cell contact information. Should be provided in the same unit as spatia data (for Visium, that is pixel). Type is float. Default=89.43
    --split = How many split sections? Type is Int. Default=0. 
    --neighborhood_threshold = Set neighborhood threshold distance in terms of same unit as spot diameter. Type is float. If not set, then it is set to four times of the spot diameter.
    --database_path = Provide your desired ligand-receptor database path here. Default database is a combination of CellChat and NicheNet database. Type is String. Default='database/CellNEST_database.csv'

Varying --filter_min_cell does not bring a big change in the output, and we recommend using a higher value for it to reduce GPU memory consumption. --threshold_gene_exp is a crucial parameter, and we recommend keeping it around 98%. Making it too low will result in too big input graph and may not fit in the GPU. On the other hand, keeping it too high causes the risk of losing important genes. --spot_diameter is set for Visium Spatial Transcriptomics data. If other data format is used, please see the vignette for details. All other parameters can be kept at the default setting unless special circumstances arise, like splitting the input graph using --split to accommodate a very large input dataset.


We have the following optional arguments for integrating intracellular pathways with the ligand-receptor coexpression and recommend keeping the parameters at their default values:
### Arguments for integrating intracellular signaling with inter-cellular signals
    --intra_database_path = Provide your desired ligand-receptor database path here. Default database is a postprocessed NicheNet database as explained in the paper. Type is String. Default='database/nichenet_pathways_CellNEST.csv'. 
    --add_intra = Set it to 1 for intracellular signaling pathway. Type is Int. Default=1
    --num_hops = Maximum number of hops for intracellular signaling pathway search. Type is Int. Default=10
    --threshold_gene_exp_intra = Threshold percentile for gene expression. Genes above this percentile are considered active. This should be kept very low to detect most of the intracellular signals. Type is float. Default=20% 

## Model Run 
Sample command:
```
nohup cellnest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 60000 --model_name='CellNEST_V1_Human_Lymph_Node_spatial' --run_id=1 > output_human_lymph_node_run1.log &
```
### Arguments
    --data_name = Name of the dataset. Type is String. Required.  
    --model_name = Provide a model name. Type is String. Required. 
    --run_id = Please provide a running ID, for example: 0, 1, 2, etc. Five runs are recommended. Type is Int.
    #=========================== default is set ======================================
    --num_epoch = Number of epochs or iterations for model training. Number of epochs or iterations for model training. Type is Int. Default=60000
    --model_path = Path to save the model state. Type is String. Default='model/'  
    --embedding_path = Path to save the node embedding and attention scores. Type is String. Default='embedding_data/'
    --hidden = Hidden layer dimension or dimension of node embedding. Type is Int. Default=512
    --training_data = Path to input graph. Type is String. Default='input_graph/'
    --heads = Number of heads in the attention model. Type is Int. Default=1
    --dropout = Dropout value for model training. Type is Float. Default=0.
    --lr_rate = Model training parameter. Type is Float. Default=0.00001
    --manual_seed = Set 'Yes' if you want to input your own seed for random number generation. Type is String. Default='no'. 
    --seed = Input the seed if --manual_seed=Yes. Type is Int. 
    --split = Set 1 to split the input graph to accommodate the memory. Type is Int. Default=0. 
    --total_subgraphs = If --split = 1, then input the number of desired subgraph. Type is Int. Default=1
    --metadata_to = Path to load the metadata. Type is String. Default='metadata/'
    #=========================== optional ======================================
    --load = Set 1 to load a previously saved model state. Type is Int. Default=0.  
    --load_model_name = Provide the model name that you want to reload. Type is String. Default='None'

Most of the parameters can be kept at default. Set --manual_seed='Yes' if you want to reproduce the results using user defined --seed. For details explanation on the required arguments please see the vignette.  

## Model Prediction Ensemble 
Sample command:
```
cellnest postprocess --data_name='V1_Human_Lymph_Node_spatial' --model_name='CellNEST_V1_Human_Lymph_Node_spatial' --total_runs=5 
```
### Arguments
    --data_name = The name of dataset. Type is String. Required.
    --model_name = Name of the trained model. Type is String. Required.
    --total_runs = How many runs for ensemble (at least 2 are preferred). Type is Int. Required.
    #######################################################################################################
    --embedding_path = Path to grab the attention scores from.  Type is String. Default='embedding_data/'
    --metadata_from =  Path to grab the metadata. Type is String. Default='metadata/' 
    --data_from = Path to grab the input graph from (to be passed to GAT).  Type is String. Default='input_graph/'
    --output_path = Path to save the visualization results, e.g., histograms, graph etc.  Type is String. Default='output/'
    --top_percent = Top N percentage communications to pick. Type is Float. Default=20.
    --cutoff_MAD = Set it 1 to keep the communication having deviation lower than MAD. Type is Int. Default=-1
    --cutoff_z_score = Set it 1 to keep the communication having z_score higher than 1.97. Type is Int. Default=-1
    --output_all = Set it to 1 to output all communications. Type is Int. Default=-1
    
By default top 20% (--top_percent=20) CCC are kept, and the rest are discarded as most of the true connections are detected within top 20% based on synthetic data. Increasing this value will result in a lot of false positives. If you want to filter based on median absolute deviation (MAD) or z scores (1.97), please use --cutoff_MAD=1 or --cutoff_z_score=1, respectively. 

## Confidence Interval
Sample command: 
```
cellnest confidence_interval --data_name='V1_Human_Lymph_Node_spatial' --model_name='CellNEST_V1_Human_Lymph_Node_spatial'
```
Please note that to run this command, you have to first run previous command (postprocess) with argument: --output_all = 1
### Arguments
    --ccc_list_path = Path to *_allCCC.csv file generated by previous command. Type is Str. Default = 'output/' 
    --model_name = Name of the trained model. Type is Str. Required = True
    --data_name = Name of the data. Type is Str. Required = True.
    --output_path = Path to save the visualization results, e.g., histograms, graph etc. Type is Str. Default = 'output/' 
    --top20 = Set to 1 to print the 95th percent confidence interval of top20 cutoff and output the CCC having the attention score within that range. Type is Int. Default = -1
    --std = Set to 1 to print the 95th percent confidence interval of standard deviation and output the CCC having the attention score within that range. Type is Int. Default = -1 

## Output Visualize 
### Arguments
    --data_name = The name of dataset. Type is String.  Required.
    --model_name = Name of the trained model. Type is String. Required.
    --top_edge_count = Number of the top communications to plot. To plot all insert -1. Type is Int. Default=1500.
    --top_percent = Top N percentage communications to pick. Type is Int. Default=20    
    --metadata_from = Path to grab the metadata. Type is String. Default='metadata/' 
    --output_path = Path to save the visualization results, e.g., histograms, graph etc. Type is String. Default='output/'
    --barcode_info_file = Path to load the barcode information file produced during data preprocessing step. Type is String. Default path points to the metadata folder.
    --annotation_file_path = Path to load the annotation file in csv format (if available). Type is String. Default path points to the metadata folder.
    --selfloop_info_file = Path to load the selfloop information file produced during data preprocessing step. Type is String. Default path points to the metadata folder.
    --top_ccc_file = Path to load the selected top CCC file produced during data postprocessing step. Type is String. Default path points to the output folder.
    --output_name = Output file name prefix according to user's choice. Type is String.  Default path points to the metadata folder.
    --filter = Set --filter=1 if you want to filter the CCC. Type is Int. Default=0
    --filter_by_ligand_receptor = Set ligand-receptor pair, e.g., --filter_by_ligand_receptor="CCL19-CCR7" if you want to filter the CCC by LR pair. Type is String. 
    --filter_by_annotation = Set cell or spot type, e.g., --filter_by_annotation="T-cell" if you want to filter the CCC. Type is String.
    --filter_by_component = Set component id, e.g., --filter_by_component=9 if you want to filter by component id. Type is String. Type is Int.  Default=-1
    --sort_by_attentionScore = Set --histogram_attention_score=1 if you want to sort the histograms of CCC by attention score. Type is String. Default=-1

--top_edge_count is a crucial parameter, and the effect of varying this parameter is explained with figures in the vignette. Additionally, --filter should be set 1 if you want to filter the CCC by user-specified ligand-receptor, spot/cell type annotation, or components, which are all explained in vignette. All other parameters can be kept at default.  

### Output:
We provide several visualization outputs to best quantify CellNEST's predictions:
1. An HTML file showing the graph of connected components on the tissue. 
2. A CSV file containing a list of CCCs with assigned connected components.
3. A NetworkX visualization of the CCC network.
4. A DOT file for plotting the graph in easy-to-share formats such as PDF and SVG with tools like Graphviz.
5. A Histogram CSV and corresponding figure output delineating abundant CCC within each connected component of tissue.
6. A scatter plot of the mean top downstream signaling gene expressions of a user-specified receptor gene.



## Downstream gene expression plot for a given receptor gene
Sample command:
```
cellnest downstream --adata_path='data/V1_Human_Lymph_Node_spatial/V1_Human_Lymph_Node_filtered_feature_bc_matrix' --positions_path='data/V1_Human_Lymph_Node_spatial/spatial/tissue_positions_list.csv' --gene='CCR7' 
```
### Arguments
    --adata_path = The path of gene expression matrix. Type is Str. Required=True 
    --positions_path = The path of position file. Type is Str. Required=True
    --gene = Gene name to plot downstream TF. Type is Str. Required=True.
    --ppi_path = Path to ppi database. Type is Str. Default='database/human_signaling_ppi.csv'  
    --tf_path = Path to TF list. Type is Str. Default='database/human_tf_target.csv'
    --output_file = Path to save the visualization. Type is Str. Default=''
 

## Relay network related commands
### Extract relay
Sample command:
```
cellnest relay_extract --data_name='V1_Human_Lymph_Node_spatial' --metadata='metadata/' --top_ccc_file='output/V1_Human_Lymph_Node_spatial/V1_Human_Lymph_Node_spatial_ccc_list_top3000.csv' --output_path='relay_validation_sample_data/lymph_node/'
```
### Arguments
    --data_name = The name of dataset. Type is Str. Required = True.  
    --metadata = The name of dataset. Type is Str. Default='metadata/'.
    --barcode_info_file = Path to load the barcode information file produced during data preprocessing step Type is Str. Default=''.
    --annotation_file_path = Path to load the annotation file in csv format (if available). Type is Str. Default=''
    --selfloop_info_file = Path to load the selfloop information file produced during data preprocessing step. Type is Str. Default='' 
    --top_ccc_file = Path to load the selected top CCC file produced during data postprocessing step. Type is Str. Required = True
    --output_path = Output file name prefix according to user\'s choice. Type is Str. Default='CellNEST_figures_output/'

Output: An HTML file showing histograms of detected relay counts and a CSV file showing the corresponding counts in a tabular format.

### Cell type identification for relays
Sample command:
```
cellnest relay_celltype --input_path='relay_validation_sample_data/lymph_node/' --output_path='CellNEST_figures_output/' --annotation_file='relay_validation_sample_data/lymph_node/fractional_abundances_by_spot.csv' --modality='spot'
```
### Arguments
    --input_path' = Directory containing CellNEST relay outputs. Type is Str. Required = True 
    --output_path' = Directory to write output plots to. Type is Str. Required = True 
    --annotation_file = Path to CSV file with cell type annotations. Type is = argparse.FileType('r'). Required = True 
    --modality' = Spatial modality with choices = ["sc", "spot"]. Type is Str. Required = True
    --additional_network = Append additional network to bar chart/create additional pie chart. Type is Str. Required = False

Output: HTML pie and bar charts to visualize cell-type proportions in relay networks. A CSV file containing cell-type proportions in a tabular format. 

### Confidence score for the CellNEST detected relays
Sample command:
```
cellnest relay_confidence --input_path='relay_validation_sample_data/lymph_node/' --output_path='CellNEST_figures_output/' --organism='human' --database_dir='database/'
```
### Arguments
    --input_path = Path to a CSV file containing relay network outputs from CellNEST. Must contain column 'Relay Patterns'. Type is  Str. Required = True 
    --database_dir = Directory containing PPI and TF-target gene databases. Type is  Str. Required = True 
    --organism = Organism profiled in spatial transcriptomics experiment with choices = ["human", "mouse"]. Type is  Str. Required = True 
    --output_path = Path to csv file to write confidence scoring output. Type is  Str. Required = True 

Output: A CSV file containing putative intracellular pathways and confidence scores associated with relay networks.
   
