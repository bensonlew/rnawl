# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
import pandas as pd
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class DbrdaWorkflow(Workflow):
    """
    dbRDA小工具
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DbrdaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "method", "type": "string", "default": "bray_curtis"},
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table, meta.otu.group_table"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.ellipse = self.add_tool("graph.ellipse")
        self.matrix = self.add_tool('meta.beta_diversity.distance_calc')
        self.dbrda = self.add_tool('meta.beta_diversity.dbrda')
    
    def check_option(self):
        """
        参数检查
        """
        if not self.option("method"):
            raise OptionError("必须设置距离算法")
        if not self.option("otutable"):
            raise OptionError("必须输入OTU表")
        if not self.option("envtable"):
            raise OptionError("必须输入环境因子表")
        if not self.option("group_table"):
            raise OptionError("必须输入分组表")
        
    def run_tools(self):
        self.otu_table = self.option('otutable').prop['path']
        self.matrix.set_options({'method' :self.option('method'),
                                         'otutable': self.otu_table})
        self.matrix.on('end', self.dbrda_run)
        self.matrix.run()
    
    def dbrda_run(self):
        self.dbrda.set_options ({"otutable": self.otu_table,
                         'envtable': self.env_remove(), 
                         'method':self.option('method'),
                         'dis_matrix': self.matrix.option('dis_matrix'),
                         'ellipse': 'T',
                         'group_table':self.option("group_table")})
        with open(self.option("group_table").prop['path'],'r') as group:
            gl = []
            group.readline()
            while 1:
                line = group.readline()
                if not line:
                    break
                gl.append(line.rstrip().split('\t')[1])
            if len(set(gl)) >1 or len(set(gl)) < len(gl):
                self.skip_ellipse = False
        if not self.skip_ellipse:
            self.dbrda.on("end",self.run_ellipse)
        else:
            self.dbrda.on("end",self.set_output)
        self.dbrda.run()
        
    def run_ellipse(self):
        options = {
            "group_table" : self.option("group_table").prop['path'],
            "analysis" : "dbrda",
            "meta" : "dbrda",
            'pc_table' : self.dbrda.output_dir + "/db_rda_sites.xls"
            }
        self.ellipse.set_options(options)
        self.ellipse.on('end', self.set_output)
        self.ellipse.run()
    
    def env_remove(self):
        env_remove = os.path.join(self.work_dir,"env_new.txt")
        with open(self.option("envtable").prop["path"],"r") as old, open(env_remove,"w") as new:
            lines = []
            while 1:
                line = old.readline().rstrip()
                if not line:
                    break
                lines.append(line)
            new.write("\n".join(lines))
                
        return env_remove




    def set_output(self):
        self.logger.info('strat set_output as {}'.format(
            self.__class__.__name__))
        plot_species_path = self.dbrda.output_dir.rstrip('/') + '/db_rda_plot_species_data.xls'
        plot_feature_path = self.dbrda.output_dir.rstrip('/') + '/db_rda_plot_feature_data.xls'
        os.rename(plot_species_path,plot_feature_path)
        try:
            shutil.copytree(self.matrix.output_dir,
                            os.path.join(self.output_dir, "distance"))
            shutil.copytree(self.dbrda.output_dir,
                            os.path.join(self.output_dir, "dbrda"))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))
        self.set_db()

    # def replace_name(self, input):
    #     """
    #     对将修改的物种名称替换为原来的名称
    #     qingchen.zhang @20200611
    #     :return:
    #     """
    #     out_table = os.path.join(self.work_dir, "out_table.xls")
    #     with open(input, "r") as f, open(out_table, "w") as w:
    #         lines = f.readlines()
    #         w.write(lines[0])
    #         for line in lines[1:]:
    #             line = line.strip().split("\t")
    #             sp_name = line[0].split("; ")[-1].strip()
    #             if sp_name in self.rename:
    #                 line[0] = self.rename[sp_name]
    #             w.write("\t".join(line) + "\n")
    #     os.remove(input)
    #     os.rename(out_table, input)

    def set_db(self):
        # self.replace_name(os.path.join(self.output_dir, "Dbrda" ,"db_rda_plot_species_data.xls"))
        # self.replace_name(os.path.join(self.output_dir, "Dbrda" ,"db_rda_species.xls"))
        sn_group = {}
        feature_list = self.get_species_name()
        feature_index = {}
        n = 0
        for i in feature_list:
            feature_index[i] = n
            n += 1 
        with open(self.option('group_table').prop['path'],'r') as gt:
            gt.readline()
            while 1:
                line = gt.readline()
                if not line:
                    break
                fd = line.rstrip('\r\n').split('\t')
                sn_group[fd[0]] = fd[1]
        api = self.api.api("tool_lab.dbrda")
        dir_path = self.dbrda.output_dir
        api.add_detail(dir_path,self.option("main_id"),sn_group,feature_index)
        ellispe_table = os.path.join(self.ellipse.work_dir, 'ellipse_out.xls')
        api.add_ellipse(ellispe_table,self.option("main_id"))
        self.end()

    def get_species_name(self):  # 20170122 add by zhouxuan , last_modify by zhujuan 1017.10.09
        """
        判断丰度表中的物种数量是否大于30 ，如果大于30，筛选出丰度在前30的物种
        :return: 丰度为前30的物种或者 空的列表
        """
        old_abund_file_path = self.otu_table
        df = pd.DataFrame(pd.read_table(old_abund_file_path, sep='\t', index_col=0))
        df['Col_sum'] = df.apply(lambda x: x.sum(), axis=1)
        new_otu_file = df.sort_values(by=['Col_sum'], ascending=0).head(30)
        species_list = list(new_otu_file.index)
        return species_list

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(DbrdaWorkflow, self).end()

    def run(self):
        """
        运行
        """
        self.run_tools()
        super(DbrdaWorkflow, self).run()    

if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    # data = {
    #     'name': 'test_enrichbubble',
    #     'id': 'enrichbubble_' +  str(random.randint(1, 10000)),
    #     'type': 'workflow',
    #     'options': {
    #         "table_file1": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/go_pvalue.xls",
    #         "table_file2": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/go_pvalue1.xls",
    #         "group_num":2,
    #         "output_form":"intersection",
    #         "p_value" : 0.05, 
    #         "result_num":20,
    #         "show_order":"table_file1",
    #         "main_id" : "5e9e6a6017b2bf2049a81be3"
    #     }
    # }
    data = {
        'name': 'test_dbrda',
        'id': 'dbrda_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        'options': {
        "otutable": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/dbrda/otu_taxon.txt",
        "envtable": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/dbrda/env_table.txt",
        "group_table": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/dbrda/group.txt",
        "method":'canberra',
        "main_id" : "5e9e6a6017b2bf2049a81be1"
        }
    }
    wsheet = Sheet(data=data)
    wf = DbrdaWorkflow(wsheet)
    wf.run()

            