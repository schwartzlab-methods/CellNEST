# author: Deisha Paliwal
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Union
import networkx as nx
import os
from collections import deque
import argparse

def query_tf(
    grn_db: pd.DataFrame, 
    target_gene: str,
    activation_only: bool,
) -> Dict[float, List[str]]:
    """
    Query the gene regulatory network (GRN) database for transcription factors (TFs) regulating a given target gene, optionally filtering for activation interactions only.

    Parameters
    ----------
    grn_db : pd.DataFrame
        DataFrame containing the GRN database with columns 
        'source', 'target', 'mode', and 'confidence_score'.
    target_gene : str
        The target gene for which to query transcription factors.
    activation_only : bool
        If True, only consider interactions where mode == 1 (activation). 
        If False, consider all interactions. 

    Returns
    -------
    Dict[float, List[str]]
        A dictionary mapping confidence scores to lists of transcription factors that regulate the target gene with that score. 
        If no transcription factors are found, returns an empty dictionary. 
    """
    if activation_only:
        target_df = grn_db[(grn_db["target"] == target_gene) & (grn_db["mode"] == 1)]
    else:
        target_df = grn_db[grn_db["target"] == target_gene]
    if target_df.empty:
        return {}
    tf_scores = {}
    for tf, score in zip(target_df["source"], target_df["confidence_score"]):
        tf_scores.setdefault(score, []).append(tf) 

    return tf_scores

def first_weighted_path(
    graph: nx.Graph, 
    start: str, 
    end_nodes: List[str],
) -> Tuple[Union[float, None], Union[List[str], None]]:
    """
    Perform a breadth-first search to find the first path from start to any of the end_nodes, calculating the cumulative product of edge weights along the path.

    Parameters
    ----------
    graph : nx.Graph
        A NetworkX graph where edges have a 'weight' attribute representing the confidence score of the
        interaction.
    start : str
        The starting node. 
    end_nodes : List[str]
        A list of target nodes. 

    Returns
    -------
    Tuple[Union[float, None], Union[List[str], None]]   
        A tuple containing the cumulative product of edge weights along the first path found from start to any of the end_nodes, and the list of nodes in that path. 
        If no path is found, returns (None, None).
    """
    if start not in graph:
        return None, None
    queue = deque([(start, [start], 1.0)])
    while queue:
        current_node, path, current_score = queue.popleft()
        if current_node in end_nodes:
            return current_score, path
        visited = set(path)
        for neighbor in graph[current_node]:
            if neighbor not in visited:
                weight = graph[current_node][neighbor]['weight']
                new_score = current_score * weight
                new_path = path + [neighbor]
                queue.append((neighbor, new_path, new_score))

    return None, None 

def query_relay_network(
    input_path: str, 
    database_dir: str,
    organism: str, 
    output_path: str,
    activation_only: bool = True
) -> None:
    """
    Report the confidence of relay patterns by integrating confidence scores reported by PPI and GRN databases. 

    Parameters
    ----------
    input_path : str
        Path to the CSV file containing relay network outputs from CellNEST. Must contain a column
        'Relay Patterns' with format 'Ligand1-Receptor1 to Ligand2-Receptor2'.
    database_dir : str
        Directory containing the PPI and GRN databases as CSV files named '{organism}_tf_target.csv' and '{organism}_signaling_ppi.csv'.
    organism : str
        The organism profiled in the spatial transcriptomics experiment. 
    output_path : str
        Path to the CSV file where the confidence scoring output will be written. 
    activation_only : bool, optional
        If True, only consider activating interactions in the GRN database when querying for TFs. 
        Default: True.
    
    Returns
    -------
    None
        Writes the confidence scoring results to the specified output CSV file. 
        Output format: columns 'relay_pattern', 'receptor_1', 'ligand_2', 'combined_score', 'tf', and 'path'. If no valid path is found, 'combined_score' and 'tf' will be None, and 'path' will contain an explanation.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok = True) 
    relay_network = pd.read_csv(input_path)
    grn_db = pd.read_csv(os.path.join(database_dir, f"{organism}_tf_target.csv"))
    ppi_db = pd.read_csv(os.path.join(database_dir, f"{organism}_signaling_ppi.csv"))
    relay_network[['ligand_1', 'receptor_1', 'ligand_2', 'receptor_2']] = relay_network["Relay Patterns"].str.extract(r'(\w+)-(\w+) to (\w+)-(\w+)')
    if organism == "mouse":
        for col in ['ligand_1', 'receptor_1', 'ligand_2', 'receptor_2']:
            relay_network[col] = relay_network[col].str.capitalize()
    grn_targets = set(grn_db["target"])
    ppi_sources = set(ppi_db["source"])
    for _, row in relay_network.iterrows():
        receptor_1 = row["receptor_1"]
        ligand_2 = row["ligand_2"]
        source_gene = receptor_1 if receptor_1 in ppi_sources else f"{row['ligand_1']}-{receptor_1}"
        target_gene = ligand_2 if ligand_2 in grn_targets else f"{ligand_2}-{row['receptor_2']}"
        result = {
            "relay_pattern": row["Relay Patterns"],
            "receptor_1": source_gene,
            "ligand_2": target_gene,
            "combined_score": None,
            "tf": None,
            "path": ""
        }
        if target_gene not in grn_targets:
            result["path"] += "ligand 2 not found; "
        if source_gene not in ppi_sources:
            result["path"] += "receptor 1 not found; "
        if target_gene in grn_targets and source_gene in ppi_sources:
            tf_dict = query_tf(grn_db, target_gene, activation_only)
            if tf_dict:
                best_score = 0.0
                best_path = None
                all_tfs = [tf for tfs in tf_dict.values() for tf in tfs]
                for weight in np.arange(0.5, 0.0, -0.1):
                    filtered_ppi_db = ppi_db[ppi_db["experimental_score"] >= weight]
                    G = nx.DiGraph()
                    for _, ppi_row in filtered_ppi_db.iterrows():
                        G.add_edge(ppi_row["source"], ppi_row["target"], weight = ppi_row["experimental_score"])                
                    ppi_score, putative_path = first_weighted_path(G, source_gene, all_tfs)  
                    if putative_path and ppi_score is not None:
                        grn_score = next((score for score, tfs in tf_dict.items() if putative_path[-1] in tfs), None)
                        total_score = ppi_score * grn_score
                        if total_score > best_score:
                            best_score = total_score
                            best_total_score, best_path = total_score, putative_path
                if best_path:
                    result.update({
                        "path": "; ".join(best_path),
                        "combined_score": best_total_score,
                        "tf": best_path[-1]
                    })
                else:
                    result["path"] += "No valid path found."
            else:
                result["path"] += "No TFs found."
        pd.DataFrame([result]).to_csv(
            output_path,
            mode = 'a',
            header = not os.path.exists(output_path),
            index = False
        )

def main():
    parser = argparse.ArgumentParser(description = "Query confidence score of a relay network")
    parser.add_argument('--input_path', type = str, required = True, help = "Path to csv file containing relay network outputs from CellNEST. Must contain column 'Relay Patterns'")
    parser.add_argument('--database_dir', type = str, required = True, help = "Directory containing PPI and TF-target gene databases")
    parser.add_argument('--organism', type = str, required = True, help = "Organism profiled in spatial transcriptomics experiment", choices = ["human", "mouse"])
    parser.add_argument('--output_path', type = str, required = True, help = "Path to CSV file to write confidence scoring output")
    args = parser.parse_args()
    query_relay_network(
        input_path = args.input_path,
        database_dir = args.database_dir,
        organism = args.organism, 
        output_path = args.output_path
    )

if __name__ == "__main__":
    main()
