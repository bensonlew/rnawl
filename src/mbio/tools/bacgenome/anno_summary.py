# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# version 1.0
# last_modify: 2019.04.04

import os
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError



class AnnoSummaryAgent(Agent):
    """
    根据各个注释结果的基因id 合成注释汇总表
    """

    def __init__(self, parent):
        super(AnnoSummaryAgent, self).__init__(parent)
        options = [
            {"name": "files_list", "type": "string", "default":""}, # 注释表列表，用分号分割
            {"name": "files_label", "type": "string", "default":""},  #数据库列表，用分号分割
            {"name": "merge_column", "type": "string", "default":""},   # 各表用来合并的列名
            {"name": "result_path", "type": "string"},  #生成的注释总览表的路径
            {"name": "project_name", "type": "string","default":""}  #项目名称，不同的项目保留的列的信息不一样，故设这参数
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("files_list") == "":
            raise OptionError("必须设置参数files_list信息！")
        if self.option("files_label") == "":
            raise OptionError("必须设置参数files_label信息！")
        if self.option('merge_column') == "":
            raise OptionError("必须设置参数merge_column信息！")



    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AnnoSummaryAgent, self).end()


class AnnoSummaryTool(Tool):
    def __init__(self, config):
        super(AnnoSummaryTool, self).__init__(config)
        self.logger.info(self.option('files_list'))
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/Python/bin', LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/Python/lib')
        self.files = self.option("files_list").split(';')
        self.labels = self.option("files_label").split(';')


        if self.option("project_name") == 'bacgenome':
            self.keep_info = {
                "phi":['Gene ID', 'PHI ID','Pathogen Species','Host Species','Gene Function'],
                "vfdb":['Gene ID', 'VFDB ID', 'VFs','Description'],
                "card": ['Gene ID', 'ARO_Accession', 'ARO_description'],## qingchen.zhang 去掉ARO_category @20201109
                "tmhmm" : ['Gene ID', 'Number of predicted TMHs'],
                "tcdb": ['Gene ID','TCDB Description', 'TCDB Class'],
                "cazy":['Gene ID','Family','Class','Class_description'],  # 缺 family desc
                "promoter": ['Gene ID','Upstream Pos','Promoter Len'],
                "island":['Gene Id','GI No.'],
                "prephage":['Gene Id', 'Ph No.'], # 缺 'Possible Phage'
                "sec": ['Gene ID','Sec Type'],  #
                "signal_p":  ['Gene ID','D-score'] ,  #
                "signal_n":  ['Gene ID','D-score'] ,
                "two_sys": ['Gene', 'Type'] , #
                "antismash":['Gene ID', 'Type', 'Most Similar Cluster','Cluster ID'],
                "gene_ori_info":['Gene ID', "gene_ori_name","gene_ori_desc"]

            }

            self.rename_merge_column = {         #列名改成 self.option('merge_column'),统一merge的列名
                "phi":{"Gene Function":"PHI Function"},
                "vfdb": {'Description': "VFDB Description"},
                "island":{'Gene Id':'Gene ID'},
                "prephage":{'Gene Id':'Gene ID'},
                "two_sys" : {"Gene":"Gene ID", 'Type': "Two_sys_Type"},
                "antismash" : {"Type": "antismash_type", "Most Similar Cluster":"antismash_Most_similar_cluster","Cluster ID":"antismash_cluster"},
                "signal_p":{"D-score": "signal_p"},
                "signal_n":{"D-score": "signal_n"},

            }

        else:
            self.keep_info = {}
            self.rename_merge_column = {}


    def run_summary(self):
        data1 = pd.read_table(self.files[0])
        lab = self.labels[0]
        if lab in self.keep_info.keys():  # 如果某数据库名没在self.keep_info ，则保留所有的列
            data1 = data1[self.keep_info[lab]]
        if lab in self.rename_merge_column.keys():
            data1.rename(columns=self.rename_merge_column[lab], inplace = True)
        data1.drop_duplicates(self.option('merge_column'),'first',inplace=True)

        for id in range(1,len(self.labels)):
            lab = self.labels[id]
            file = self.files[id]
            if not os.path.exists(file):
                continue
            if lab == 'island':  #island 的detail 基因有重复，所有把1个基因的信息合并
                new_file = self.work_dir + '/new_island.xls'
                data = pd.read_table(file, sep='\t', header=0)
                new_data = data[['GI No.','Gene Id']]
                n2 = new_data.groupby(by='Gene Id').apply(lambda x : ','.join(x['GI No.']))
                n3=n2.reset_index()
                n3.columns=['Gene Id','GI No.']
                n3.to_csv(new_file,index=False, sep='\t')
                file = new_file
            data = pd.read_table(file, sep='\t', header=0)
            if lab in self.keep_info.keys():
                data = data[self.keep_info[lab]]
            if lab in self.rename_merge_column.keys():
                data.rename(columns=self.rename_merge_column[lab], inplace = True)

            data.drop_duplicates(self.option('merge_column'),'first',inplace=True)
            data1 = data1.merge(data, on=self.option('merge_column'),how='left')

        data1.to_csv(self.output_dir + '/summary.xls', sep='\t', index=False)
        self.option('result_path', self.output_dir + '/summary.xls')


    def run(self):
        super(AnnoSummaryTool, self).run()
        self.run_summary()
        self.end()
