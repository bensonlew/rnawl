# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'

""""""

import datetime
from biocluster.workflow import Workflow
import re
import os
import json
import shutil
from mbio.packages.meta.save_params import save_params


class EnvironmentalRegressionWorkflow(Workflow):
    """
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EnvironmentalRegressionWorkflow, self).__init__(wsheet_object)
        options = [
        
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "group_table", "type": "infile","format": "meta.otu.group_table"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "env_labs", "type": "string", "default": ""},
            {"name": "PCAlabs", "type": "string", "default": ""},
            #{"name": "method", "type": "string", "default": "sum"},
            #{"name": "problem_type", "type": "int", "default": 2},
            #{"name": "top_n", "type": "int", "default": 100},
            {"name": "environmental_regression_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": 'string'},
            {"name": "otu_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            {"name": "env_id","type":"string"},
            #{"name": "env_labs","type":"string"}
        ]
        self.add_option(options)
        #newtable = os.path.join(self.work_dir, 'otutable1.xls')
        #f2 = open(newtable, 'w+')
        #with open(tablepath, 'r') as f:

        self.set_options(self._sheet.options())
        self.environmental_regression = self.add_tool("meta.beta_diversity.environmental_regression")
        #self.samples = re.split(',', self.option("samples"))
        self.output_dir = self.environmental_regression.output_dir

    def change_otuname(self, tablepath):
        newtable = os.path.join(self.work_dir, 'otutable1.xls')
        f2 = open(newtable, 'w+')
        with open(tablepath, 'r') as f:
            i = 0
            for line in f:
                if i == 0:
                    i = 1
                    f2.write(line)
                else:
                    line = line.strip().split('\t')
                    line_data = line[0].strip().split(' ')
                    line_he = "".join(line_data)
                    line[0] = line_he
                    #line[0] = line_data[-1]
                    for i in range(0, len(line)):
                        if i == len(line)-1:
                            f2.write("%s\n"%(line[i]))
                        else:
                            f2.write("%s\t"%(line[i]))
        f2.close()
        return newtable



    def run_environmental_regression(self):
        newtable = self.change_otuname(self.option('otu_table').prop['path'])
        options = {
            'otu_table': newtable,
            #'otutable':self.option('otutable'),
            'level': self.option('level'),
            'env_labs':self.option('env_labs'),
            'group_table':self.option('group_table'),
            'env_table':self.option('envtable'),
            'PCAlabs':self.option('PCAlabs'),
            #'problem_type':self.option('problem_type'),
            #'top_n':self.option('top_n')
        }
        self.environmental_regression.set_options(options)
        self.environmental_regression.on('end',self.set_db)
        self.output_dir = self.environmental_regression.output_dir
        self.environmental_regression.run()

    def end(self):
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "environmentalregression分析结果目录"],
            ["./sites.xls", "xls", "metagenomeseq坐标数据"],
            ["./rotation.xls", "xls", "pca的全部输出"],
            ["./rotationone.xls", "xls", "一个pca的输出"],
            ["./messages.xls", "xls", "主要信息"],
            #["./linesdata.xls", "xls", "直线信息"],
            ["./importance.xls", "xls", "百分率返回"]
        ])
        super(EnvironmentalRegressionWorkflow, self).end()

    def set_db(self):
        api_environmental_regression = self.api.environmental_regression
        datasite = self.output_dir + '/sites.xls'
        datamessages = self.output_dir + '/messages.xls'
        if not os.path.isfile(datasite):
            self.logger.error("找不到报告文件:{}".format(datasite))
            self.set_error("找不到报告文件", code="12701301")
        if not os.path.isfile(datamessages):
            self.logger.error("找不到报告文件:{}".format(datamessages))
            self.set_error("找不到报告文件", code="12701301")
        api_environmental_regression.add_environmental_regression_site(file_path=datasite, table_id=self.option("environmental_regression_id"))
        api_environmental_regression.add_environmental_regression_messages(file_path=datamessages, table_id=self.option("environmental_regression_id"))
        self.end()

    def run(self):
        self.run_environmental_regression()
        super(EnvironmentalRegressionWorkflow, self).run()        

