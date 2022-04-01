# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
import pandas as pd
import numpy as np
import math
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class BarChartAgent(Agent):
    def __init__(self, parent):
        super(BarChartAgent, self).__init__(parent)
        options = [
            {"name": "data_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "group_file", "type": "infile", "format": "ref_rna_v2.common"},
            {'name': 'ishape', 'type': 'string', 'default': ''},
            {'name': 'out_file', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)
        self.step.add_steps("bar_chart")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.bar_chart.start()
        self.step.update()

    def step_finish(self):
        self.step.bar_chart.finish()
        self.step.update()

    def check_options(self):
        if not self.option("data_file").is_set:
            raise OptionError("必须设置输入数据文件")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        super(BarChartAgent, self).end()


class BarChartTool(Tool):
    def __init__(self, config):
        super(BarChartTool, self).__init__(config)
        self._version = "v1.0"
        self.data_file = ''
        self.group_file = ''

    def run(self):
        super(BarChartTool, self).run()
        self.format_check()
        self.run_data_process()
        self.set_output()
        self.end()

    def format_check(self):
        if os.path.exists(self.option('data_file').prop['path']):
            data_format = os.path.basename(self.option('data_file').prop['path']).split('.')[-1]
            self.data_file = self.format_convert(path=self.option('data_file').prop['path'],
                                                 f_format=data_format, f_type='data')
        if self.option('group_file').is_set and os.path.exists(self.option('group_file').prop['path']):
            group_format = os.path.basename(self.option('group_file').prop['path']).split('.')[-1]
            self.group_file = self.format_convert(path=self.option('group_file').prop['path'],
                                                  f_format=group_format, f_type='group')

    def format_convert(self, path, f_format, f_type):
        if f_format.lower() == 'txt':
            file_path = path
        elif f_format.lower() in ['xls', 'xlsx']:
            try:
                data_pd = pd.read_excel(path, header=None)
            except:
                file_path = path
            else:
                file_path = os.path.join(self.work_dir, '{}_file.txt'.format(f_type))
                data_pd.to_csv(file_path, sep='\t', header=False, index=False)
        else:
            file_path = ''
            self.logger.info('{}文件格式不为txt/xls，请检查'.format(f_type))
        return file_path

    def run_data_process(self):
        out_file = os.path.join(self.work_dir, 'data4bar.txt')
        if self.group_file:
            self.data_process(outfile=out_file)
        else:
            data_pd = pd.read_table(self.data_file, header=0)
            data_pd.rename(columns={data_pd.columns[0]: 'name', data_pd.columns[1]: 'value'}, inplace=True)
            data_pd.to_csv(out_file, sep='\t', header=True, index=False)

    def data_process(self, outfile):
        data_dict = dict()
        group_dict = dict()
        with open(self.data_file, 'r') as d:
            d.readline()
            for l in d:
                each = l.strip().split('\t')
                if each[0] not in data_dict.keys():
                    data_dict[each[0]] = each[1]
                else:
                    self.logger.info('第一列{}重复出现，请检查。仅保留首次出现的值。'.format(each[0]))
        with open(self.group_file, 'r') as g:
            for line in g:
                items = line.strip().split('\t')
                if items[1] not in group_dict.keys():
                    group_dict[items[1]] = list()
                group_dict[items[1]].append(items[0])
        with open(outfile, 'w') as o:
            if self.option('ishape') == 'sd':
                o.write('name\tmean\tstd\n')
                for g_name in group_dict.keys():
                    g_list = [float(data_dict[x]) for x in group_dict[g_name]]
                    mean = round(np.mean(g_list), 2)
                    std = round(np.std(g_list, ddof=1), 2)
                    o.write(g_name + "\t" + str(mean) + '\t' + str(std) + '\n')
            elif self.option('ishape') == 'sem':
                o.write('name\tmean\tstd\n')
                for g_name in group_dict.keys():
                    g_list = [float(data_dict[x]) for x in group_dict[g_name]]
                    mean = round(np.mean(g_list), 2)
                    std = round(np.std(g_list, ddof=1)/math.sqrt(len(g_list)), 2)
                    o.write(g_name + "\t" + str(mean) + '\t' + str(std) + '\n')
            else:
                o.write('name\tvalue\n')
                for g_name in group_dict.keys():
                    g_list = [float(data_dict[x]) for x in group_dict[g_name]]
                    mean = round(np.mean(g_list), 2)
                    o.write(g_name + "\t" + str(mean) + '\n')

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        data4bar = os.path.join(self.work_dir, 'data4bar.txt')
        self.option('out_file', data4bar)
