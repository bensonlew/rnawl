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


class MetagenomeseqWorkflow(Workflow):
    """
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetagenomeseqWorkflow, self).__init__(wsheet_object)
        options = [
        
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "group_table", "type": "infile","format": "meta.otu.group_table"},
            
            #{"name": "method", "type": "string", "default": "sum"},
            #{"name": "problem_type", "type": "int", "default": 2},
            #{"name": "top_n", "type": "int", "default": 100},
            {"name": "metagenomeseq_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": 'string'},
            {"name": "otu_id", "type": 'string'},
            {"name": "group_detail", "type": "string"}
        ]
        self.add_option(options)
        #newtable = os.path.join(self.work_dir, 'otutable1.xls')
        #f2 = open(newtable, 'w+')
        #with open(tablepath, 'r') as f:

        self.set_options(self._sheet.options())
        self.metagenomeseq = self.add_tool("meta.beta_diversity.metagenomeseq")
        #self.samples = re.split(',', self.option("samples"))
        self.output_dir = self.metagenomeseq.output_dir

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



    def run_metagenomeseq(self):
        newtable = self.change_otuname(self.option('otu_table').prop['path'])
        options = {
            'otu_table': newtable,
            #'otutable':self.option('otutable'),
            'level': self.option('level'),
            #'envlabs':self.option('envlabs'),
            'group_table':self.option('group_table'),
            #'method':self.option('method'),
            #'problem_type':self.option('problem_type'),
            #'top_n':self.option('top_n')
        }
        self.metagenomeseq.set_options(options)
        self.metagenomeseq.on('end',self.set_db)
        self.output_dir = self.metagenomeseq.output_dir
        self.metagenomeseq.run()

    def end(self):
        save_params(self.output_dir,self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出文件目录", 0, "110148"],
            ["./diff.xls", "xls", "坐标数据", 0, "110149"],
            #["./roc_auc.xls", "xls", "面积"]
        ])
        super(MetagenomeseqWorkflow, self).end()

    def set_db(self):
        api_metagenomeseq = self.api.metagenomeseq
        datadiff = self.output_dir + '/diff.xls'
        datalist = self.output_dir + '/list.xls'
        if not os.path.isfile(datadiff):
            self.logger.error("找不到报告文件:{}".format(datadiff))
            self.set_error("找不到报告文件", code="12702501")
        if not os.path.isfile(datalist):
            self.logger.error("找不到报告文件:{}".format(datalist))
            self.set_error("找不到报告文件", code="12702501")
        api_metagenomeseq.add_metagenomeseq_diff(file_path=datadiff, table_id=self.option("metagenomeseq_id"))
        api_metagenomeseq.add_metagenomeseq_list(file_path=datalist, table_id=self.option("metagenomeseq_id"))
        os.remove(datalist)  # 删除不上传的文件 GHD @20180918
        self.end()

    def run(self):
        self.run_metagenomeseq()
        super(MetagenomeseqWorkflow, self).run()        

