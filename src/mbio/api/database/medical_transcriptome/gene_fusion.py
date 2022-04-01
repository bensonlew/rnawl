# -*- coding: utf-8 -*-
# __author__ = 'fwy'
from __future__ import division
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from bson.objectid import ObjectId
# from cStringIO import StringIO
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
from bson.son import SON
import datetime
from mbio.api.database.medical_transcriptome.api_base import ApiBase
import json
import pandas as pd
import os
import math
import unittest


class GeneFusion(ApiBase):
    def __init__(self, bind_object):
        super(GeneFusion, self).__init__(bind_object)
        self._project_type = 'medical_transcriptome'
        #self._db_name = Config().MONGODB + '_medical_transcriptome'

    def add_fusion_main(self, fusion_result=None,pos_file= None,chr_length=None, circos = None ,task_id=None, main_id=None,params=None, method_type=None, new_output=None,project_sn=None,s3_output = None, group=None):
        """
        导入gene_fusion主表信息
        :return:
        """
        if circos:
            circos = 1
        else:
            circos = 0
        if main_id is None:
            # prepare main table info
            name = "Gene_Fusion" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='snp_indel analysis main table',
                params=params,
                status="start",
                version="v3"
            )
            fusion_id = self.create_db_table('sg_gene_fusion', [main_info])
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            fusion_id = ObjectId(main_id)

        if fusion_result:
            if circos == 1:
                # chr_df = pd.read_table(chr_length,header = None )
                # chr_df.columns = ["chr", "chr_length"]
                # chr_dict = chr_df.set_index("chr").to_dict("index")
                sample_list = self.add_fusion_detail(fusion_result=fusion_result, fusion_id = fusion_id, pos_file =pos_file, circos = circos,chr_length=chr_length, group=group)
                result_dict = OrderedDict()
                for sample in sample_list:
                    result_dict[sample] = os.path.join(s3_output,sample,"star-fusion.fusion_predictions.abridged.tsv")
                    # result_dict[sample] = os.path.join(self.bind_object.sheet.output,'03Gene_structure_analysis', "03GeneFusion","GeneFusion_Star_Fusion","Star_fusion",sample,"star-fusion.fusion_predictions.abridged.tsv")
                self.update_db_record('sg_gene_fusion', fusion_id, status="end",sample =sample_list,has_file = "yes",sample_files = result_dict,
                                      main_id=fusion_id, ciros="true")
            else:
                sample_list = self.add_fusion_detail(fusion_result=fusion_result,fusion_id = fusion_id, circos = circos, group=group)
                result_dict = OrderedDict()
                for sample in sample_list:
                    result_dict[sample] = os.path.join(s3_output, sample, "star-fusion.fusion_predictions.abridged.tsv")
                    # result_dict[sample] = os.path.join(self.bind_object.sheet.output, '03Gene_structure_analysis',
                    #                                    "03GeneFusion", "GeneFusion_Star_Fusion", "Star_fusion", sample,
                    #                                    "star-fusion.fusion_predictions.abridged.tsv")
                self.update_db_record('sg_gene_fusion', fusion_id, status="end", sample=sample_list,has_file = "yes",sample_files = result_dict,
                                      main_id=fusion_id)

        return fusion_id

    def add_fusion_detail(self,fusion_result=None,fusion_id=None,pos_file=None,circos=None,chr_length=None, group=None):
        """
            导入Gene_fusion详情表的函数
            :param fusion_result: 融合分析模块的结果文件夹，即~/output_dir/
            :param fusion_id: Gene_fusion主表的ID
            :return:
            """
        result_dir = os.listdir(fusion_result)
        sample_list=[]
        if circos == 1:
            chr_df = pd.read_table(chr_length, header=None)
            chr_df.columns = ["chr", "chr_length"]
            chr_df["gene_fusion_id"] = fusion_id
            chr_df["attributes"] = "chromosome"
            chr_up_dict = chr_df.to_dict("r")
            self.create_db_table('sg_gene_fusion_circos_detail', chr_up_dict)
        for i in result_dir:
            if i == "fusion_stat.txt":
                fusion_stat = os.path.join(fusion_result,"fusion_stat.txt")
                self.add_fusion_stat(fusion_stat, fusion_id, group=group)
            elif i == "star_fusion":
                sample_result_dirs = os.listdir(os.path.join(fusion_result,i))
                for sample in sample_result_dirs:
                    sample_name = sample
                    sample_list.append(sample_name)
                    sample_detail_result = os.path.join(fusion_result, "star_fusion",sample_name, "star-fusion.fusion_predictions.abridged.tsv")
                    detail_sdf = pd.read_table(sample_detail_result)
                    if detail_sdf.shape[0] > 0:
                        select_columns = ["#FusionName", "JunctionReadCount", "SpanningFragCount", "LeftGene",
                                          "LeftBreakpoint", "RightGene", "RightBreakpoint", "SpliceType", "FFPM"]
                        select_df = detail_sdf[select_columns]
                        select_df.rename({"#FusionName": "FusionName"}, axis=1, inplace=True)
                        select_df.columns = map(str.lower, select_df.columns)
                        select_df["gene_fusion_id"] = fusion_id
                        select_df["sample"] = sample_name
                        select_df["fusion_unique_id"] = select_df["leftbreakpoint"] + select_df["rightbreakpoint"]
                        print(fusion_id)
                        print(select_df.columns)
                        detail_info = select_df.to_dict('records')
                        self.create_db_table('sg_gene_fusion_detail', detail_info)
                        if circos == 1:
                            self.add_circos_detail(sample_name, select_df, pos_file, fusion_id,chr_length=chr_length)
            else:
                sample_name = i
                sample_list.append(i)
                sample_detail_result = os.path.join(fusion_result,i,"star-fusion.fusion_predictions.abridged.tsv")
                detail_sdf = pd.read_table(sample_detail_result)
                if detail_sdf.shape[0] > 0:
                    select_columns = ["#FusionName","JunctionReadCount","SpanningFragCount","LeftGene","LeftBreakpoint","RightGene","RightBreakpoint","SpliceType","FFPM"]
                    select_df = detail_sdf[select_columns]
                    select_df.rename({"#FusionName":"FusionName"},axis=1,inplace=True)
                    select_df.columns = map(str.lower, select_df.columns)
                    select_df["gene_fusion_id"] = fusion_id
                    select_df["sample"] = sample_name
                    select_df["fusion_unique_id"] = select_df["leftbreakpoint"]+select_df["rightbreakpoint"]
                    print(fusion_id)
                    print(select_df.columns)
                    detail_info = select_df.to_dict('records')
                    self.create_db_table('sg_gene_fusion_detail', detail_info)
                    if circos == 1:
                        self.add_circos_detail(sample_name,select_df,pos_file,fusion_id,chr_length=chr_length)
        return sorted(sample_list)


    def add_circos_detail(self,sample_name,select_df,pos_file,fusion_id,chr_length=None):
        # print("sample_name")
        # print(sample_name)
        # print(select_df)
        # print(select_df.columns)
        # select_df["sample"] = sample_name

        print(select_df.columns)
        chr_df = pd.read_table(chr_length,header = None )
        chr_df.columns = ["chr", "chr_length"]
        chr_df["gene_fusion_id"] = fusion_id
        chr_df["attributes"] = "chromosome"
        chr_up_dict = chr_df.to_dict("r")
        # self.create_db_table('sg_gene_fusion_circos_detail', chr_up_dict)
        chr_dict = chr_df.set_index("chr").to_dict("index")
        pos_df = pd.read_table(pos_file,header = None)
        pos_df.columns = ["gene_id","chr","start","end"]
        pos_df = pos_df.set_index("gene_id")
        pos_dict = pos_df.to_dict("index")
        # circos_select_columns = ["sample","FusionName","LeftGene","LeftBreakpoint","RightGene","RightBreakpoint","FFPM"]
        circos_select_columns = ["sample", "fusionname", "leftgene", "leftbreakpoint", "rightgene", "rightbreakpoint",
                                 "ffpm"]
        circos_select_df = select_df[circos_select_columns]
        circos_select_df[["left_chr", "left_local"]] = circos_select_df['leftbreakpoint'].str.split(':',expand=True).iloc[:,:2]
        circos_select_df[["right_chr", "right_local"]] = circos_select_df['rightbreakpoint'].str.split(':',expand=True).iloc[:,:2]
        circos_select_df[["left_start", "left_end"]] = circos_select_df["leftgene"].map(pos_dict).apply(pd.Series).loc[:, ["start","end"]]
        circos_select_df[["right_start", "right_end"]] = circos_select_df["rightgene"].map(pos_dict).apply(pd.Series).loc[:, ["start","end"]]
        circos_select_df["left_local"] = circos_select_df["left_local"].astype("int")
        circos_select_df["right_local"] = circos_select_df["right_local"].astype("int")
        if circos_select_df.shape[0] != 1 :
            if max(circos_select_df["ffpm"]) == min(circos_select_df["ffpm"]):
                circos_select_df["width"] = 3
            else:
                circos_select_df["width"] = circos_select_df.apply(lambda x: 4 * (x["ffpm"] - min(circos_select_df["ffpm"])) / (max(circos_select_df["ffpm"]) - min(circos_select_df["ffpm"])) + 1, axis=1)
        else:
            circos_select_df["width"] = 3
        circos_select_df["attributes"] = "line"
        circos_select_line_df = circos_select_df.drop(["leftbreakpoint", "rightbreakpoint","ffpm"], axis=1)
        circos_select_line_df["gene_fusion_id"] = fusion_id
        circos_select_line_detail = circos_select_line_df.to_dict("r")
        self.create_db_table('sg_gene_fusion_circos_detail', circos_select_line_detail)

        #这一段代码负责通过基因断点位置将融合基因的情况进行分类
        data_list = []
        all_genes  = circos_select_df["leftbreakpoint"].tolist() + circos_select_df["rightbreakpoint"].tolist()
        tree_dict=defaultdict(list)
        for index,item in circos_select_df.iterrows():
            left_gene_name= item["fusionname"].split("-")[0]
            left_gene_chr = item["left_chr"]
            left_gene_loc = item["left_local"]
            left_area_loc =int(math.ceil( left_gene_loc/chr_dict[left_gene_chr]["chr_length"] *100 ))
            left_detail_info=[(left_area_loc-1)*chr_dict[left_gene_chr]["chr_length"]/100,left_gene_loc,left_gene_name]
            tree_dict[left_gene_chr+"_"+str(left_area_loc)].append(left_detail_info)
            right_gene_name = item["fusionname"].split("-")[-1]
            right_gene_chr = item["right_chr"]
            right_gene_loc = item["right_local"]
            right_area_loc = int(math.ceil(right_gene_loc / chr_dict[right_gene_chr]["chr_length"] * 100))
            right_detail_info = [(right_area_loc - 1) * chr_dict[right_gene_chr]["chr_length"] / 100, right_gene_loc, right_gene_name]
            tree_dict[right_gene_chr + "_"+str(right_gene_loc)].append(right_detail_info)
        for cluster in tree_dict:
            sample_name = sample_name
            chr = cluster.split("_")[0]
            root = tree_dict[cluster][0][0]
            branch = [i[1] for i in tree_dict[cluster]]
            gene_name = [i[2] for i in tree_dict[cluster]]
            data = [
                ('attributes', "tree"),
                ('sample', sample_name),
                ('chr', chr),
                ('root', root),
                ('branch', branch),
                ('gene_name',gene_name),
                ('gene_fusion_id',fusion_id)
            ]
            data = SON(data)
            data_list.append(data)
        collection = self.db['sg_gene_fusion_circos_detail']
        collection.insert_many(data_list)
        # except Exception as e:
        #     print("导入circos图信息出错")
        # else:
        #     print("导入circos图信息成功!")


    def add_fusion_stat(self,fusion_stat,fusion_id, group=None):
        if group:
            if os.path.isfile(group):
                sample_list = list()
                with open(group, 'r') as g:
                    for line in g.readlines():
                        if line.startswith('#'):
                            continue
                        sample_list.append(line.strip().split('\t')[0])
            else:
                sample_list = group.split(';')

        stat_df = pd.read_table(fusion_stat)
        stat_df["gene_fusion_id"] = fusion_id
        stat_dict = stat_df.to_dict("r")
        if group:
            stat_dict.sort(key=lambda x: sample_list.index(x['Sample']))
        self.create_db_table('sg_gene_fusion_stat', stat_dict)


    def run(self):
        fusion_result="/mnt/ilustre/users/sanger-dev/workspace/20200820/GeneFusion_tsg_38312_1068_2677/GeneFusion/output"
        pos_file = "/mnt/ilustre/users/sanger-dev/workspace/20200820/GeneFusion_tsg_38312_1068_2677/gene_pos"
        chr_length = "/mnt/ilustre/users/sanger-dev/workspace/20200820/GeneFusion_tsg_38312_1068_2677/chr_length"
        # chr_length = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/data_pre/make_lib/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/chrNameLength.txt",
        self.add_fusion_main(fusion_result=fusion_result, pos_file=pos_file,chr_length=chr_length,
                          task_id="medical_transcriptome", project_sn="medical_transcriptome",circos=1)


    def add_gene_fusion_venn(self, venn_graph, venn_table=None, project_sn='medical_transcriptome', main_id=None,
                      task_id='medical_transcriptome', params=None):
        """
        add venn analysis info
        :param venn_graph: venn_graph.xls resulted from express_venn tool
        :param venn_table: venn_table.xls resulted from express_venn tool
        :param quant_method: exp quant method
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :param exp_level: T or G
        :param params: string of parameters dict, designed for judgement of whether the task is repeated.
        :return: main table id
        """
        # add main table info
        if main_id is None:
            name = "Fusion_Result_Venn" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='venn main table',
                params=params,
                status="start"
            )
            main_id = self.create_db_table('sg_gene_fusion_venn', [main_info])
        else:
            main_id = ObjectId(main_id)
        # add detail table info
        graph_pd = pd.read_table(venn_graph, header=0, sep='\t', keep_default_na=False)
        graph_pd.columns = ["name", "ids"]
        detail_dict_list = graph_pd.to_dict('records')
        if venn_table:
            table_pd = pd.read_table(venn_table, header=None, sep='\t', keep_default_na=False)
            table_pd.columns = ["combination", "num", "only_list"]
            detail_dict_list += table_pd.to_dict('records')
        self.create_db_table('sg_gene_fusion_venn_detail', detail_dict_list, tag_dict={'gene_fusion_venn_id': main_id})
        self.update_db_record('sg_gene_fusion_venn', main_id, status="end", main_id=main_id)
        return main_id

    def run2(self):
        fusion_venn_resut="/mnt/ilustre/users/sanger-dev/workspace/20200818/FusionVenn_tsg_158346_6079_3784/FusionVenn/output/fusion_venn_graph.xls"
        self.add_gene_fusion_venn(fusion_venn_resut)


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
      from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
      from biocluster.wsheet import Sheet
      import random

      data = {
        # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
        "id": "medical_transcriptome",
        "project_sn": "medical_transcriptome",
        "type": "workflow",
        "name": "medical_transcriptome.medical_transcriptome_test_api",
        "options": {
        },
      }
      wsheet = Sheet(data=data)
      wf = MedicalTranscriptomeTestApiWorkflow(wsheet)
      wf.IMPORT_REPORT_DATA = True
      wf.IMPORT_REPORT_AFTER_END = False
      wf.test_api = wf.api.api("medical_transcriptome.gene_fusion")
      wf.test_api.run()
      wf.test_api.run2()
if __name__ == '__main__':
    unittest.main()


