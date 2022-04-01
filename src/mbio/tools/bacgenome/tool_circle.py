# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __last_modify__ = '2021.07.27'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd


class ToolCircleAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(ToolCircleAgent, self).__init__(parent)
        options = [
            {"name": "file_dir", "type": "string"},
            {"name": "data_type", "type": "string"},
        ]
        self.add_option(options)

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "5G"


class ToolCircleTool(Tool):
    def __init__(self, config):
        super(ToolCircleTool, self).__init__(config)
        self.python = '/program/Python/bin/python'

    def run_cog_total(self):
        """
        description
        :return:
        """
        list2 = []
        funs = []
        samples = []
        for file in os.listdir(self.option("file_dir")):
            list1 =[]
            sample = file.split(".cog.xls")[0]
            samples.append(sample)
            with open(self.option("file_dir") + "/" + file, "r") as f:
                lines = f.readlines()
                for line in lines[1:]:
                    lin = line.strip().split("\t")
                    des =lin[1].split(']')[1].strip()
                    list1.append([des, lin[2]])
                    funs.append(des)
            df = pd.DataFrame(list1, columns=['fun', sample])
            df = df.set_index('fun')
            list2.append(df)
        df2 = pd.DataFrame(list(set(funs)), columns=['fun'])
        df2 = df2.set_index('fun')
        list3 =[df2]+list2
        all_data = pd.concat(list3, axis=1, join='outer')
        all_data = all_data.fillna(0)
        all_data['fun'] = all_data.index
        names = ['fun'] + samples
        all_data.to_csv(self.output_dir+ "/all.table.xls", sep='\t', index=False, columns=names)


    def run_go_total(self):
        """
        description
        :return:
        """
        list1 =[]
        list2 = []
        samples = []
        for file in os.listdir(self.option("file_dir")):
            sample = file.split(".go.xls")[0]
            samples.append(sample)
            data =pd.read_table(self.option("file_dir") + "/" +file, sep='\t', header=0)
            data = data[data['GO (Lev1)'] != '-']
            data= data.groupby(data['GO (Lev1)'], as_index=False).agg({"Seq Number": 'sum'})
            data.columns = ['fun', sample]
            des = list(data['fun'])
            for i in des:
                list2.append(i)
            data = data.set_index('fun')
            list1.append(data)
        df2 = pd.DataFrame(list(set(list2)), columns=['fun'])
        df2 = df2.set_index('fun')
        list3 = [df2] + list1
        all_data = pd.concat(list3, axis=1, join='outer')
        all_data = all_data.fillna(0)
        all_data['fun'] = all_data.index
        names = ['fun'] + samples
        all_data.to_csv(self.output_dir+ "/all.table.xls", sep='\t', index=False, columns=names)

    def run_kegg_total(self):
        """
        description
        :return:
        """
        list1 = []
        list2 = []
        samples = []
        for file in os.listdir(self.option("file_dir")):
            sample = file.split(".kegg.xls")[0]
            samples.append(sample)
            data = pd.read_table(self.option("file_dir") + "/" + file, sep='\t', header=0)
            data = data[data['Level1'] != '-']
            data = data.groupby(data['Level1'], as_index=False).agg({"Gene nu": 'sum'})
            data.columns = ['fun', sample]
            des = list(data['fun'])
            for i in des:
                list2.append(i)
            data = data.set_index('fun')
            list1.append(data)
        df2 = pd.DataFrame(list(set(list2)), columns=['fun'])
        df2 = df2.set_index('fun')
        list3 = [df2] + list1
        all_data = pd.concat(list3, axis=1, join='outer')
        all_data = all_data.fillna(0)
        all_data['fun'] = all_data.index
        names = ['fun'] + samples
        all_data.to_csv(self.output_dir+ "/all.table.xls", sep='\t', index=False, columns=names)

    def run(self):
        super(ToolCircleTool, self).run()
        if self.option("data_type") in ['COG']:
            self.run_cog_total()
        elif  self.option("data_type") in ['GO']:
            self.run_go_total()
        elif  self.option("data_type") in ['KEGG']:
            self.run_kegg_total()
        self.end()