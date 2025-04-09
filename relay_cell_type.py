# Written By
# Deisha Paliwal
import yaml 
import pandas as pd
import os
import gzip
import glob
import pickle
from typing import Dict, Any, List
import argparse
import re
import altair as alt
import altairThemes 
alt.themes.register("publishTheme", altairThemes.publishTheme)
alt.themes.enable("publishTheme")

# determine number of legend columns 
def n_legend_cols(value):
    if 1 < value <= 10:
        return 1
    elif 11 < value <= 20:
        return 2
    elif 21 < value <= 31:
        return 3
    else:
        return 4

# split relay network name into list: [L1, R1_L2, R2]
def get_relay_pairs(
        network: str
    ):
    str = re.sub(r'\s', '', network)
    str_list = re.split(r'[to,-]', str)
    del str_list[2]
    lr_list = [str_list[0], f"{str_list[1]}_{str_list[2]}", str_list[3]] # [L1, R1_L2, R2] 

    return lr_list

# create a pie chart for each individual network
def piechart(
        relay: dict, 
        network: str, 
        ctp: pd.DataFrame,
        modality: str
    ):
    df = query_ctp(relay, network, ctp, modality)
    chart_list = [] # for horizontal concatenation 
    for col in df.columns: # df.columns: [0,1,2] - triplet of spots/cells in relay network 
        ctp_data = df[[col]].reset_index()
        ctp_data.columns = ['cell_type', 'proportion'] # index -> cell_type
        relay_pairs = get_relay_pairs(network)
        n_cell_types = len(ctp_data["cell_type"].unique())
        set1 = altairThemes.get_colour_scheme("Set1", n_cell_types)
        n_cols =  n_legend_cols(n_cell_types)
        chart = alt.Chart(ctp_data).mark_arc().encode(
            theta = 'proportion:Q',
            color = alt.Color('cell_type:N', scale = alt.Scale(range = set1), legend = alt.Legend(title="Cell type", columns = n_cols, symbolLimit = 0)),
            tooltip = ['cell_type:N', 'proportion:Q']
        ).properties(
            title = relay_pairs[col],
            width = 130, 
            height = 130
        )
        chart_list.append(chart)
    combined_chart = alt.hconcat(*chart_list).configure_concat(
        spacing = 2
    ) # concatenate horizontally

    return combined_chart, df

# create a box plot for top 10 networks 
def boxplot(
        relay: dict, 
        network_list: List[str], 
        ctp: pd.DataFrame,
        modality: str
    ):
    df_list = []
    for network in network_list:
        df = query_ctp(relay, network, ctp, modality)
        df = df.stack().reset_index(drop = False) # df.stack(): indices remain unchanged - column names become their own column
        df.rename(columns={0: 'proportion'}, inplace = True)
        network_comps = get_relay_pairs(network)
        df["network"] = network
        df_list.append(df)
    df_all = pd.concat(df_list).reset_index(drop = True)
    n_cell_types = len(df_all["cell_type"].unique())
    set1 = altairThemes.get_colour_scheme("Set1", n_cell_types)
    n_cols = n_legend_cols(n_cell_types)
    chart = alt.Chart(df_all).mark_bar().encode(
        x = alt.X("idx:N", title = "cell" if modality == "sc" else "spot"),
        y = alt.Y("proportion:Q", title = "Proportion"),
        column = alt.Column("network:N", title = None, spacing = 45),
        color = alt.Color("cell_type:N", scale = alt.Scale(range = set1), legend = alt.Legend(title = "Cell type", columns = n_cols, symbolLimit = 0)),
        tooltip = ['cell_type:N', 'proportion:Q']
    ).properties(
            width = 120, 
            height = 150
        )

    return chart 

# prepare average proportion df for spot-based and cell-based ST data 
def query_ctp(
        relay: dict, 
        network: str, 
        annot_df: pd.DataFrame, 
        modality: str
    ):
    instances = relay[network]
    result_dict = {}
    for instance in instances:
        for idx, item in enumerate(instance):
            if idx not in result_dict:
                if modality == "sc":
                    result_dict[idx] = {ct: 0 for ct in annot_df["annotation"].unique()}
                elif modality == "spot":
                    result_dict[idx] = {ct: [] for ct in annot_df.columns}
            id = item[0]
            if id in annot_df.index:
                if modality == "sc":
                    annotation = annot_df.loc[id, "annotation"]
                    result_dict[idx][annotation] += 1
                elif modality == "spot":
                    proportion = annot_df.loc[id]
                    for ct in annot_df.columns:
                        result_dict[idx][ct].append(proportion[ct])
    result = []
    for idx, vals in result_dict.items():
        if modality == "sc":
            total_counts = sum(vals.values())
            for ct, count in vals.items():
                prop = count / total_counts if total_counts > 0 else 0 
                result.append({'idx': idx, 'cell_type': ct, 'proportion': prop})
        elif modality == "spot":
            for ct, props in vals.items():
                avg_prop = sum(props) / len(props) if props else 0
                result.append({'idx': idx, 'cell_type': ct, 'proportion': avg_prop})
    result_df = pd.DataFrame(result)
    pivoted_df = result_df.pivot(index = "cell_type", columns = "idx", values = "proportion")

    return pivoted_df

# create pie and bar charts from user input 
def make_plots(
        input_dir: str,
        output_dir: str,
        annotation_file: str,
        modality: str,
        additional_network: str = None
    ):
    os.makedirs(output_dir, exist_ok = True)
    relay_ranked = pd.read_csv(glob.glob(os.path.join(input_dir, "*relay*count.csv"))[0])
    top_networks = relay_ranked.iloc[:, 0].head(5).tolist()
    if additional_network:
        top_networks.append(additional_network)
    relay_cell_info = glob.glob(os.path.join(input_dir, "*pattern_cell_info"))[0]
    with gzip.open(relay_cell_info, "rb") as fp:
        relay = pickle.load(fp)
    annot_df = pd.read_csv(annotation_file, index_col = 0)
    box_plot = boxplot(relay, top_networks, annot_df, modality)
    box_plot.save(os.path.join(output_dir, "bar.html"))
    for network in top_networks:
        scrubbed_network = network.replace(" ", "_").replace('"', '')
        pie_chart, relay_proportion_df = piechart(relay, network, annot_df, modality)
        pie_chart.save(os.path.join(output_dir, f"{scrubbed_network}_pie.html"))
        relay_proportion_df.to_csv(os.path.join(output_dir, f"{scrubbed_network}_proportions.csv"))

def main():
    parser = argparse.ArgumentParser(description="Visualize cell types participating in relay networks")
    parser.add_argument('--input_path', type = str, required = True, help = "Directory containing CellNEST relay outputs")
    parser.add_argument('--output_path', type = str, required = True, help = "Directory to write output plots to")
    parser.add_argument('--annotation_file', type = argparse.FileType('r'), required = True, help = "Path to csv file with cell type annotations")
    parser.add_argument('--modality', type = str, required = True, help = "Spatial modality", choices = ["sc", "spot"])
    parser.add_argument('--additional_network', type = str, required = False, help = "Append additional network to bar chart/create additional pie chart")
    args = parser.parse_args()

    make_plots(
        input_dir = args.input_dir,
        output_dir = args.output_dir,
        annotation_file = args.annotation_file,
        modality = args.modality,
        additional_network = args.additional_network
    )

if __name__ == "__main__":
    main()

########## sample input ##########################
'''
relay_cell_type.py --input_dir='relay_validation_sample_data/lymph_node/' --output_dir='CellNEST_figures_output/' 
    --annotation_file='relay_validation_sample_data/lymph_node/fractional_abundances_by_spot.csv' --modality='spot' 
'''

