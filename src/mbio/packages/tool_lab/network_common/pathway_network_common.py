import json
import pandas as pd
import os

def get_cytoscape_file(network_json,output_dir):
    a = json.load(open(network_json))
    node_infos = a["nodes_id_infos"]
    link_infos = a["link_infos"]
    node_df = pd.DataFrame(node_infos)
    link_df = pd.DataFrame(link_infos)
    map_ids = list(node_df["id"])
    node_infos_df = node_df.set_index("id")
    node_infos_dict = node_infos_df.to_dict("index")
    new_link_df = link_df
    new_link_df["id_1"] = new_link_df["source"].apply(lambda x: map_ids[x])
    new_link_df["id_1_name"] = new_link_df["id_1"].apply(lambda x: node_infos_dict[x]["name"])
    new_link_df["id_2"] = new_link_df["target"].apply(lambda x: map_ids[x])
    new_link_df["id2_name"] = new_link_df["id_2"].apply(lambda x: node_infos_dict[x]["name"])
    new_link_df = new_link_df[["id_1", "id_1_name", "id_2", "id2_name", "type"]]
    all_used_ids = set(new_link_df["id_1"]) | set(new_link_df["id_2"])
    new_nodes_df = node_df[node_df["id"].isin(all_used_ids)]
    new_nodes_df["degree"] = new_nodes_df.apply(
        lambda x: (new_link_df == x["id"]).sum().sum() if x["type"] == "actual" else "-", axis=1)
    new_nodes_df = new_nodes_df[["id","name","type","size","degree","pvalue"]]
    new_nodes_df=new_nodes_df.rename(columns={"size":"gene_num"})
    nodes_file = os.path.join(output_dir,"nodes.txt")
    new_nodes_df.to_csv(nodes_file,sep="\t",index=False)
    edges_file = os.path.join(output_dir,"edges.txt")
    new_link_df.to_csv(edges_file,sep="\t",index=False)
