# singularity run --home=nest_container  \
# nest_container/nest_image.sif \
# nohup python -u preprocess_flowsig.py \
# --model_name=NEST_PDAC_exp2_D1_64630_manualDB \
# --ccc_file=output/PDAC_exp2_D1_64630_manualDB/NEST_PDAC_exp2_D1_64630_manualDB_top20percent.csv \
# --data_from=../../data/notta_pdac_visium/spaceranger_outputs_partial_no_header/exp2_D1/outs/ > output.log &


import scanpy as sc
import pandas as pd
from collections import defaultdict
from scipy.sparse import csr_matrix
import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( '--model_name', type=str, default='NEST_PDAC_exp2_D1_64630_manualDB', 
                        help='Name of the trained model') #, required=True
    parser.add_argument( '--ccc_file', type=str, default='output/PDAC_exp2_D1_64630_manualDB/NEST_PDAC_exp2_D1_64630_manualDB_top20percent.csv', 
                        help='Path to grab the CellNEST detected intercellular comunication list.')
    parser.add_argument( '--data_from', type=str, default='../../data/notta_pdac_visium/spaceranger_outputs_partial_no_header/exp2_D1/outs/', 
                        help='Path to grab the raw count')
    parser.add_argument( '--output_path', type=str, default='output/', 
                        help='Path to save the visualization results, e.g., histograms, graph etc.')
    parser.add_argument( '--database_path', type=str, default='database/CellNEST_database_no_predictedPPI.csv', 
                        help='Provide your desired ligand-receptor database path here. ' \
                        'Default database is a combination of CellChat and NicheNet database. It should have Ligand and Receptor columns') 

    args = parser.parse_args()
    #adata = sc.read('/content/drive/MyDrive/PostDoc/flowSig/Mouse_embryo_svg_E9.5.h5ad')
    adata = sc.read_visium(path=args.data_from, count_file='filtered_feature_bc_matrix.h5')
    adata.var_names_make_unique()
    #adata.raw = adata
    adata.layers['count'] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]


    barcode_info = list(adata.obs.index)
    total_cell = adata.X.shape[0]

    cellnest_ccc_adata = pd.read_csv(args.ccc_file, sep=",")
    lr_pair_dict = defaultdict(list)
    gene_id = dict()
    for i in range(0, len(cellnest_ccc_adata)):
        from_cell = cellnest_ccc_adata['from_cell'][i]
        to_cell = cellnest_ccc_adata['to_cell'][i]
        ligand = cellnest_ccc_adata['ligand'][i]
        receptor = cellnest_ccc_adata['receptor'][i]
        gene_id[ligand] = 1
        gene_id[receptor] = 1
        attention_score = cellnest_ccc_adata['attention_score'][i]
        lr_pair_name = ligand + '-' + receptor
        lr_pair_dict[lr_pair_name].append([from_cell, to_cell, attention_score])


    ##########################################################
    df = pd.read_csv(args.database_path, sep=",")
    '''
            Ligand   Receptor          Annotation           Reference
    0        TGFB1     TGFBR1  Secreted Signaling      KEGG: hsa04350
    1        TGFB1     TGFBR2  Secreted Signaling      KEGG: hsa04350
    '''
    print('ligand-receptor database reading done.')
    uns_df = defaultdict(list)
    uns_df['ligand'] = list(df['Ligand'])
    uns_df['receptor'] = list(df['Receptor'])
    for i in range (0, len(df)):
        uns_df['pathway'].append('unknown')

    uns_df = pd.DataFrame(uns_df)
    uns_df =uns_df.drop_duplicates()
    uns_df = uns_df.reset_index(drop=True)
    adata.uns['commot-cellchat-info'] = {}
    adata.uns['commot-cellchat-info']['df_ligrec'] = uns_df
    #######################################################
    # obsm_sum_sender = defaultdict(list)
    #adata.obsm['commot-cellchat-sum-sender']
    # this is a df with cell barcode as index, and columns as s-lr_pair_entry
    # if cell A is acting as a sender for LR pair L1R1, then, 
    # obsm_sum_sender[A]['s-L1R1'] = attention score of that

    obsm_sum_sender = defaultdict(list)
    obsm_sum_receiver= defaultdict(list)
    cellname_to_index = dict()
    i = 0
    for cell_barcode in barcode_info:
        obsm_sum_sender['cellname'].append(cell_barcode)
        obsm_sum_receiver['cellname'].append(cell_barcode)
        cellname_to_index[cell_barcode] = i
        i = i+1

    for i in range(0, len(uns_df)):
        #if uns_df['ligand'][i] not in gene_id or uns_df['receptor'][i] not in gene_id:

        lr_pair_name = uns_df['ligand'][i]+ '-' + uns_df['receptor'][i]
        s_lrpair = 's-' + lr_pair_name
        r_lrpair = 'r-' + lr_pair_name
        obsm_sum_sender[s_lrpair] = []
        obsm_sum_receiver[r_lrpair] = []
        for cell_barcode in barcode_info:
            obsm_sum_sender[s_lrpair].append(0)
            obsm_sum_receiver[r_lrpair].append(0)


        # make a total_cell x total_cell matrics
        temp_x = np.zeros((total_cell, total_cell))
        for links in lr_pair_dict[lr_pair_name]:
            from_cell = links[0]
            to_cell = links[1]
            score = links[2]
            
            from_indx = cellname_to_index[from_cell]
            to_indx = cellname_to_index[to_cell]

            temp_x[from_indx][to_indx] = score   
            obsm_sum_sender[s_lrpair][from_indx] = obsm_sum_sender[s_lrpair][from_indx] + score     
            obsm_sum_receiver[r_lrpair][to_indx] = obsm_sum_receiver[r_lrpair][to_indx] + score


        temp_x = csr_matrix(temp_x, dtype=np.float64)
        adata.obsp['commot-cellchat-'+lr_pair_name] = temp_x


    obsm_sum_sender = pd.DataFrame(obsm_sum_sender)
    obsm_sum_sender = obsm_sum_sender.set_index('cellname') # change the index column of them
    obsm_sum_receiver = pd.DataFrame(obsm_sum_receiver)
    obsm_sum_receiver = obsm_sum_receiver.set_index('cellname') # change the index column of them



    adata.obsm['commot-cellchat-sum-receiver'] = obsm_sum_receiver
    adata.obsm['commot-cellchat-sum-sender'] = obsm_sum_sender

    # save it. 
    adata.write(args.output_path + '/' + args.model_name + '_flowsig_preprocessed.h5ad', compression="gzip")
    print('saved: args.output_path' + '/' + args.model_name + '_flowsig_preprocessed.h5ad')
