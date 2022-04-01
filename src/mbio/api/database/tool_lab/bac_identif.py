# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from bson.objectid import ObjectId
import datetime
import os
import json
import types
from mbio.api.database.whole_transcriptome.api_base import ApiBase
from Bio import Phylo


class BacIdentif(ApiBase):
    def __init__(self, bind_object):
        super(BacIdentif, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_table(self, project_sn='tool_lab', main_id=None, task_id='tool_lab', params=None):
        if main_id is None:
            name = "strain_identification" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='strain_identification',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('sg_bac_id', [main_info])
        else:
            main_id = ObjectId(main_id)
        return main_id

    def add_detail(self, main_id, stat_file, taxon_file, dir, dir2, undefined=None):
        self.add_stat(stat_file, main_id, undefined)
        tree_dict = self.add_taxon(taxon_file, dir, main_id)
        table_inf = json.dumps(tree_dict, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_bac_id', main_id, main_id=main_id, tree_type=table_inf)
        self.add_anno(dir2, main_id, undefined)
        tree_16s = {"name": "16s_tree", "condition": {}}
        table_info = json.dumps(tree_16s, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_bac_id', main_id, main_id=main_id, s16_tree=table_info)
        house_tree = {"name": "housekeeping_tree", "condition": {}}
        table_info1 = json.dumps(house_tree, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_bac_id', main_id, main_id=main_id, house_tree=table_info1)
        tree_check = [{"field": "branch_len", "title": "分支枝长"}, {"field": "bootstraps", "title": "bootstraps 值"}]
        self.update_db_record('sg_bac_id', main_id, status="end", main_id=main_id, tree_check=tree_check)


    def add_stat(self, file, main_id, anl_type=None):
        data_list = []
        if anl_type:
            table_dict = {"column": [
                {"field": "sample", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                {"field": "s16", "filter": "false", "sort": "false", "title": "16S序列", "type": "int"},
                {"field": "house_keeping", "filter": "false", "sort": "false", "title": "看家基因", "type": "int"},
                {"field": "type", "filter": "false", "sort": "false", "title": "Type", "type": "string"}],
                "condition": {}}
            with open(file, "r") as f:
                lines = f.readlines()
                for line in lines:
                    lin = line.strip().split("\t")
                    type = ""
                    if len(lin) == 3:
                        type = "Target"
                    elif len(lin) == 4:
                        type = "Reference"
                    insert_data = {
                        'strain_id': main_id,
                        "sample": lin[0],
                        "s16": lin[1],
                        "house_keeping": lin[2],
                        "type": type,
                    }
                    data_list.append(insert_data)
        else:
            table_dict = {"column": [
                {"field": "sample", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
                {"field": "s16", "filter": "false", "sort": "false", "title": "16S序列", "type": "int"},
                {"field": "house_keeping", "filter": "false", "sort": "false", "title": "看家基因", "type": "int"}],
                "condition": {}}
            with open(file, "r") as f:
                lines = f.readlines()
                for line in lines:
                    lin = line.strip().split("\t")
                    insert_data = {
                        'strain_id': main_id,
                        "sample": lin[0],
                        "s16": lin[1],
                        "house_keeping": lin[2],
                    }
                    data_list.append(insert_data)
        self.create_db_table('sg_strain_stat', data_list)
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_bac_id', main_id, main_id=main_id, table_data_sum=table_info)

    def add_taxon(self, file, dir, main_id):
        dict = {}
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        data_list = []
        table_dict = {"column": [
            {"field": "sample", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
            {"field": "taxon", "filter": "false", "sort": "false", "title": "Most Similar taxonomy（GTDB)",
             "type": "string"},{"field": "s16_tree", "filter": "false", "sort": "false", "title": "s16_tree", "type": "string"},{"field": "house_tree", "filter": "false", "sort": "false", "title": "house_tree", "type": "string"}], "condition": {}}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines:
                list_sam=[]
                lin = line.strip().split("\t")
                insert_data = {
                    "strain_id": main_id,
                    "sample": lin[0],
                    "taxon": lin[1],
                    "direction": "h"
                }
                if os.path.exists(dir + "/" + lin[0] + ".16s.nwk"):
                    samples = self.get_sample_from_tree(dir + "/" + lin[0] + ".16s.nwk")
                    dict1 = [{"groupname": "all", "value": samples}]
                    insert_data['group'] = dict1
                    with open(dir + "/" + lin[0] + ".16s.nwk", "r") as f:
                        lines = f.readlines()
                        insert_data['s16_tree'] = lines[0].strip()
                        list_sam.append({"field": "s16_tree", "title": "16s"})
                if os.path.exists(dir + "/" + lin[0] + ".house_keeping.nwk"):
                    with open(dir + "/" + lin[0] + ".house_keeping.nwk", "r") as f:
                        lines = f.readlines()
                        insert_data['house_tree'] = lines[0].strip()
                        list_sam.append({"field": "house_tree", "title": "看家基因"})
                data_list.append(insert_data)
                dict[lin[0]] =list_sam
        self.create_db_table('sg_strain_classify', data_list)
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_bac_id', main_id, main_id=main_id, table_data_taxon=table_info)
        return dict

    def add_anno(self, dir, main_id, anl_type=None):
        if len(os.listdir(dir)) > 0:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
            data_list = []
            if anl_type:
                for dd in os.listdir(dir):
                    with open(dir + "/" + dd, "r") as f:
                        lines = f.readlines()
                        if len(lines[0].strip().split("\t")) == 4:
                            table_dict = {"column": [
                                {"field": "sample", "filter": "false", "sort": "false", "title": "Sample name",
                                 "type": "string"},
                                {"field": "ref_id", "filter": "false", "sort": "false", "title": "Reference ID",
                                 "type": "string"},
                                {"field": "s16_similarity", "filter": "false", "sort": "false",
                                 "title": "16SrRNA gene sequence similarity(%)", "type": "float"},
                                {"field": "ncbi_taxon", "filter": "false", "sort": "false", "title": "NCBI taxonomy",
                                 "type": "string"}],
                                "condition": {}}
                        elif len(lines[0].strip().split("\t")) == 5:
                            table_dict = {"column": [
                                {"field": "sample", "filter": "false", "sort": "false", "title": "Sample name",
                                 "type": "string"},
                                {"field": "ref_id", "filter": "false", "sort": "false", "title": "Reference ID",
                                 "type": "string"},
                                {"field": "s16_similarity", "filter": "false", "sort": "false",
                                 "title": "16SrRNA gene sequence similarity(%)", "type": "float"},
                                {"field": "ani", "filter": "false", "sort": "false", "title": "ANI(%)",
                                 "type": "float"},
                                {"field": "ncbi_taxon", "filter": "false", "sort": "false", "title": "NCBI taxonomy",
                                 "type": "string"}],
                                "condition": {}}
                        for line in lines[1:]:
                            lin = line.strip().split("\t")
                            if len(lin) == 4:
                                insert_data = {
                                    'strain_id': main_id,
                                    "sample": lin[0],
                                    "ref_id": lin[1],
                                    "s16_similarity": lin[2],
                                    "ncbi_taxon": lin[3],
                                }
                                data_list.append(insert_data)
                            elif len(lin) == 5:
                                insert_data = {
                                    'strain_id': main_id,
                                    "sample": lin[0],
                                    "ref_id": lin[1],
                                    "s16_similarity": lin[2],
                                    "ani": lin[3],
                                    "ncbi_taxon": lin[4],
                                }
                                data_list.append(insert_data)
            else:
                for dd in os.listdir(dir):
                    with open(dir + "/" + dd, "r") as f:
                        lines = f.readlines()
                        if len(lines[0].strip().split("\t")) == 5:
                            table_dict = {"column": [
                                {"field": "sample", "filter": "false", "sort": "false", "title": "Sample name",
                                 "type": "string"},
                                {"field": "ref_id", "filter": "false", "sort": "false", "title": "Reference ID",
                                 "type": "string"},
                                {"field": "s16_similarity", "filter": "false", "sort": "false",
                                 "title": "16SrRNA gene sequence similarity(%)", "type": "float"},
                                {"field": "gtdb_taxon", "filter": "false", "sort": "false", "title": "GTDB taxonomy",
                                 "type": "string"},
                                {"field": "ncbi_taxon", "filter": "false", "sort": "false", "title": "NCBI taxonomy",
                                 "type": "string"}],
                                "condition": {}}
                        elif len(lines[0].strip().split("\t")) == 6:
                            table_dict = {"column": [
                                {"field": "sample", "filter": "false", "sort": "false", "title": "Sample name",
                                 "type": "string"},
                                {"field": "ref_id", "filter": "false", "sort": "false", "title": "Reference ID",
                                 "type": "string"},
                                {"field": "s16_similarity", "filter": "false", "sort": "false",
                                 "title": "16SrRNA gene sequence similarity(%)", "type": "float"},
                                {"field": "ani", "filter": "false", "sort": "false", "title": "ANI(%)",
                                 "type": "float"},
                                {"field": "gtdb_taxon", "filter": "false", "sort": "false", "title": "GTDB taxonomy",
                                 "type": "string"},
                                {"field": "ncbi_taxon", "filter": "false", "sort": "false", "title": "NCBI taxonomy",
                                 "type": "string"}],
                                "condition": {}}
                        for line in lines[1:]:
                            lin = line.strip().split("\t")
                            if len(lin) == 5:
                                insert_data = {
                                    'strain_id': main_id,
                                    "sample": lin[0],
                                    "ref_id": lin[1],
                                    "s16_similarity": lin[2],
                                    "gtdb_taxon": lin[3],
                                    "ncbi_taxon": lin[4],
                                }
                                data_list.append(insert_data)
                            elif len(lin) == 6:
                                insert_data = {
                                    'strain_id': main_id,
                                    "sample": lin[0],
                                    "ref_id": lin[1],
                                    "s16_similarity": lin[2],
                                    "ani": lin[3],
                                    "gtdb_taxon": lin[4],
                                    "ncbi_taxon": lin[5],
                                }
                                data_list.append(insert_data)
            self.create_db_table('sg_strain_detail', data_list)
            table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
            self.update_db_record('sg_bac_id', main_id, main_id=main_id, table_data=table_info)

    def get_sample_from_tree(self, file):
        list = []
        tree = Phylo.read(file, "newick")
        for i in tree.get_terminals():
            list.append(str(i))
        return list