# Written By
# Tingxiao Gao
# Edited by 
# Fatema Tuz Zohora

import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import altair as alt
from PIL import Image
import matplotlib.pyplot as plt
import sys
import altairThemes
# Set Altair theme
alt.themes.register("publishTheme", altairThemes.publishTheme)
alt.themes.enable("publishTheme")


def load_data(ppi_path, tf_path):
    """
    Load PPI and TF matrices.

    Parameters
    ----------
    ppi_path : str
        Path to the PPI matrix file.
    tf_path : str
        Path to the TF matrix file.

    Returns
    -------
    pd.DataFrame, pd.DataFrame
        Loaded PPI and TF dataframes.
    """
    ppi_matrix = pd.read_csv(ppi_path)
    tf_matrix = pd.read_csv(tf_path)
    return ppi_matrix, tf_matrix


def extract_downstream_genes(ppi_matrix, source_gene, top_percent=20):
    """
    Extract top downstream genes based on experimental scores.

    Parameters
    ----------
    ppi_matrix : pd.DataFrame
        PPI matrix containing 'source', 'target', and 'experimental_score' columns.
    source_gene : str
        Source gene to extract downstream genes from.
    top_percent : int, optional
        Percentage of top-scored genes to include, by default 20.

    Returns
    -------
    list
        List of top downstream gene names.
    """
    immediate_genes = ppi_matrix[ppi_matrix["source"] == source_gene]["target"].tolist()
    downstream_genes_with_scores = []
    for gene in immediate_genes:
        targets = ppi_matrix[ppi_matrix["source"] == gene][["target", "experimental_score"]]
        downstream_genes_with_scores.extend(targets.values.tolist())
    downstream_genes_with_scores.sort(key=lambda x: x[1], reverse=True)
    top_index = int(len(downstream_genes_with_scores) * (top_percent / 100))
    top_genes = [gene[0] for gene in downstream_genes_with_scores[:top_index]]
    return list(set(top_genes + immediate_genes))


def load_spatial_data(adata_path, positions_path):
    """
    Load spatial data and align with the AnnData object.

    Parameters
    ----------
    adata_path : str
        Path to the 10x h5 file.
    positions_path : str
        Path to the spatial positions CSV file.

    Returns
    -------
    anndata.AnnData
        AnnData object with spatial information.
    """
    adata = sc.read_10x_h5(adata_path)
    adata.var_names_make_unique()
    positions = pd.read_csv(positions_path, header=None)
    positions.columns = ["barcode", "col1", "col2", "col3", "x", "y"]
    positions.set_index("barcode", inplace=True)
    positions = positions.reindex(adata.obs_names)
    adata.obsm["spatial"] = positions[["x", "y"]].values
    adata.obsm["spatial"] = adata.obsm["spatial"][:, [1, 0]]  # Flip coordinates
    return adata


def preprocess_adata(adata, genes):
    """
    Preprocess AnnData object and calculate average expression of specified genes.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object.
    genes : list
        List of genes to include in the calculations.

    Returns
    -------
    pd.DataFrame
        DataFrame containing spatial coordinates and gene expression.
    """
    genes = [gene for gene in genes if gene in adata.var_names]
    adata.obs["CCR7_downstream_20_avg_log"] = np.log1p(adata[:, genes].X.mean(axis=1).A1)
    adata.obs["x"] = adata.obsm["spatial"][:, 0]
    adata.obs["y"] = adata.obsm["spatial"][:, 1]
    plot_data = adata.obs[["x", "y", "CCR7_downstream_20_avg_log"]].copy()
    return plot_data


def create_spatial_expression_plot(data, title, output_file):
    """
    Create a spatial expression plot using Altair.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame with x, y, and expression level columns.
    title : str
        Plot title.
    output_file : str
        Path to save the plot.
    """
    data["y"] = -data["y"]
    scatter_plot = (
        alt.Chart(data)
        .mark_circle(size=14)
        .encode(
            x=alt.X("x:Q", title=None, axis=alt.Axis(grid=False, ticks=False, labels=False, domain=False)),
            y=alt.Y("y:Q", title=None, axis=alt.Axis(grid=False, ticks=False, labels=False, domain=False)),
            color=alt.Color("CCR7_downstream_20_avg_log:Q", scale=alt.Scale(scheme="magma"), title="Expression Level"),
            tooltip=["CCR7_downstream_20_avg_log"]
        )
        .properties(width=250, height=250, title=title)
        .configure_view(stroke=None, strokeWidth=0)
    )
    scatter_plot.save(output_file)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( '--adata_path', type=str, help='The path of gene expression matrix', required=True) # 
    parser.add_argument( '--positions_path', type=str, help='The path of position file', required=True)
    parser.add_argument( '--gene', type=str, help='Gene name to plot downstream TF', required=True)
    parser.add_argument( '--ppi_path', type=str, default='database/human_signaling_ppi.csv' ,help='Path to ppi database') # 
    parser.add_argument( '--tf_path', type=str, default='database/human_tf_target.csv', help='Path to TF list') 
    parser.add_argument( '--top_percent', type=int, default=20, help='Top N percentage to plot')    
    parser.add_argument( '--output_file', type=str, default='', help='Path to save the visualization')
    args = parser.parse_args()
    if args.output_file == '':
        args.output_file = 'output/'+'downstreamGene_'+args.gene+'.html'

    print("Reading database...")
    ppi_matrix, _ = load_data(args.ppi_path, args.tf_path)
    print("Extracting downstream genes for the input gene...")
    top_genes = extract_downstream_genes(ppi_matrix, args.gene)
    print("Loading spatial data...")
    adata = load_spatial_data(args.adata_path, args.positions_path)
    print("Plotting...")
    plot_data = preprocess_adata(adata, top_genes)
    create_spatial_expression_plot(plot_data, "Mean Expression of Top 20% "+ args.gene +" Downstream Genes", args.output_file)
    print("Plot saved at %s"%args.output_file)

