from mapman_packages import *
import os
import json
from mbio.packages.project_demo.interaction_rerun.interaction_delete import linkfile
import shutil

xml_file_dir = "/mnt/lustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/mapman/test_data_one/b_xmls"
svg_file_dir = "/mnt/lustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/mapman/test_data_one/b_imgs"
target_pathway_dir = "/mnt/lustre/users/sanger-dev/app/database/Tool_lab/mapman_e44/pathways"
target_relation_dir = "/mnt/lustre/users/sanger-dev/app/database/Tool_lab/mapman_e44/common_relations"
all_xmls = os.listdir(xml_file_dir)
bind_id2pathway = dict()
pathway2bind_id = dict()
pathway_bind_id_pos = dict()
all_pathways = []

for xml in all_xmls:
    xml_file_path = os.path.join(xml_file_dir,xml)
    pathway_name,bind_ids,id2pos = extract_bind_id_from_xml(xml_file_path)
    all_pathways.append(pathway_name)
    for bind_id in bind_ids:
        bind_id2pathway[bind_id] =  pathway_name
    pathway2bind_id[pathway_name] =bind_ids
    pathway_bind_id_pos[pathway_name] = id2pos

with open(os.path.join(target_relation_dir,"pathway2bind_id"),"w") as f:
    json.dump(pathway2bind_id, f, sort_keys=True, indent=4,default=set_default)

with open(os.path.join(target_relation_dir,"bind_id2pathway"),"w") as f:
    json.dump(bind_id2pathway, f, sort_keys=True, indent=4,default=set_default)

for pathway in all_pathways:
    if os.path.exists(os.path.join(target_pathway_dir,pathway)):
        shutil.rmtree(os.path.join(target_pathway_dir,pathway))
    os.makedirs(os.path.join(target_pathway_dir,pathway))
    pathway_dir = os.path.join(target_pathway_dir,pathway)
    print("{} has done".format(pathway))
    linkfile(os.path.join(svg_file_dir,pathway+".svg"),os.path.join(pathway_dir,pathway+".svg"))
    linkfile(os.path.join(xml_file_dir, pathway + ".xml"), os.path.join(pathway_dir, pathway + ".xml"))

    with open(os.path.join(target_pathway_dir,pathway,"pathway_infos"),"w") as f:
        pathway_infos = dict()
        pathway_infos["name"] = pathway
        pathway_infos["svg"] = os.path.join(pathway_dir,pathway+".svg")
        pathway_infos["xml"] = os.path.join(pathway_dir, pathway + ".xml")
        pathway_infos["bind_ids"] = pathway2bind_id[pathway]
        pathway_infos["id_pos"] = pathway_bind_id_pos[pathway]
        json.dump(pathway_infos, f, sort_keys=True, indent=4,default=set_default)
    print("{} success".format(pathway))



