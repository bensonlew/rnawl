# -*- coding: utf-8 -*-
import argparse
from mapman_packages import *
import json
import pandas as pd
import shutil
import traceback

class Mapman(object):
    def __init__(self,plot_type = '',species="",version="",output="",database="",exp=""):
        print("开始初始化")
        self.database = database
        self.output = output
        self.plot_type = plot_type
        self.version = version
        self.species = species
        self.exp =exp
        self.exp_df =  pd.read_table(self.exp)
        self.samples =  list(self.exp_df.columns[1:])
        self.exp_df.columns = ["seq_id"]+self.samples
        self.all_gene_ids = list(self.exp_df["seq_id"])

        self.f_exp_df = self.exp_df.set_index("seq_id")
        self.species_info_file = os.path.join(database,"species",species,version,"species_infos")
        self.species_info = json.load(open( self.species_info_file))
        self.mapping_file_dir = os.path.join(database,"pathways")
        self.common_info_dir = os.path.join(database,"common_relations")
        self.bind_id2pathway =  json.load(open( os.path.join(self.common_info_dir,"bind_id2pathway")))
        self.pathway2bind_id = json.load(open(os.path.join(self.common_info_dir, "pathway2bind_id")))

        self.all_bind_ids = self.get_all_bindids()
        self.pathways =self.get_all_pathways()
        self.check_options()

    def check_options(self):
        if len(set(self.all_gene_ids) - set(self.species_info["g2b"].keys())) > 0:
            print("gene_id 和物种存在不对应的情况")
            raise Exception("gene_id 和物种存在不对应的情况")

    def get_all_bindids(self):
        g2b = self.species_info["g2b"]
        all_bind_ids = set()
        for gene_id in self.all_gene_ids:
            bind_ids = g2b.get(gene_id)
            if bind_ids:
                for bind_id in bind_ids:
                    all_bind_ids.add(bind_id)
        return all_bind_ids

    def get_all_pathways(self):
        all_pathways = set()
        for bind_id in self.all_bind_ids:
            try:
                pathway = self.bind_id2pathway[bind_id]
                all_pathways.add(pathway)
            except:
                pass
            if len(all_pathways) >= 62:
                break
        return all_pathways

    def plot_all_pathways(self):
        for pathway in self.pathways:
            pathway_dir = os.path.join(self.mapping_file_dir,pathway)
            pathway_infos = json.load(open(os.path.join(self.mapping_file_dir,pathway,"pathway_infos")))
            pathway_bind_ids = self.pathway2bind_id[pathway]
            used_bind_ids = set(pathway_bind_ids) & self.all_bind_ids
            svg_file = pathway_infos["svg"]
            id_pos = pathway_infos["id_pos"]
            if os.path.exists(os.path.join(self.output,pathway)):
                shutil.rmtree(os.path.join(self.output,pathway))
            os.makedirs(os.path.join(self.output,pathway))
            os.makedirs(os.path.join(self.output,pathway,"draw_svgs"))
            draw_svg_dir = os.path.join(self.output,pathway,"draw_svgs")
            all_draw_svgs = []
            for bind_id in used_bind_ids:
                draw_svg_infos = {}
                bind_id_gene_ids = self.species_info["b2g"][bind_id]
                used_gene_ids = set(self.species_info["b2g"][bind_id]) & set(bind_id_gene_ids)
                draw_svg_df = self.f_exp_df.loc[list(used_gene_ids),]
                draw_line_svgs(work_dir= draw_svg_dir,data_df=draw_svg_df,svg_name=bind_id+".svg" )
                draw_svg_infos["file_path"] = os.path.join(draw_svg_dir,bind_id+".svg")
                draw_svg_infos["x_coord"] = id_pos[bind_id][0]
                draw_svg_infos["y_coord"] = id_pos[bind_id][1]
                all_draw_svgs.append(draw_svg_infos)
            insert_chart_svg_to_basic(basic_svg=svg_file, draw_svg_infos=all_draw_svgs,
                                      target_file_path=os.path.join(self.output, pathway,pathway+".svg"))
            # try:
            #     insert_chart_svg_to_basic(basic_svg = svg_file,draw_svg_infos=draw_svg_infos,target_file_path=os.path.join(self.output,pathway))
            #     print("{}绘制成功".format(pathway))
            # except Exception as e:
            #     print("error :{}\n".format(e))
            #     print("{}绘制失败".format(pathway))
            #     raise Exception("{}绘制失败".format(pathway))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-exp', type=str, default=None,
                        help="path of expression value table, tab as separator."
                             " If None, the second column of count_matrix must be gene length which"
                             " will be used to calculate fpkm or tpm. NOTE: Expression table "
                             "has nothing to do with differential analysis; Only used in report.")
    parser.add_argument('-plot_type', type=str, default="line", help='line or heat. Default: line')
    parser.add_argument('-species', type=str, default="Arabidopsis_thaliana", help='species_name')
    parser.add_argument('-version', type=str, default="Ath_AFFY_ATH1_TAIR9_Jan2010", help='annot_version')
    parser.add_argument('-output', type=str, default=None, help='output directory.')
    parser.add_argument('-database', type=str, default=None, help='database directory.')


    args = parser.parse_args()
    toolbox = Mapman(exp =args.exp,
                     plot_type = args.plot_type,
                     species = args.species,
                     version = args.version,
                     output= args.output,
                     database= args.database
                     )
    toolbox.plot_all_pathways()
