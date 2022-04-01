# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class AnosimWorkflow(Workflow):
    """
    Anosim分析小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnosimWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"otu_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name":"group_table","type": "infile", "format": 'toolapps.group_table,meta.otu.group_table'},
            {"name":"method","type":"string", "default": "bray_curtis"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.matrix = self.add_tool('meta.beta_diversity.distance_calc')
        self.anosim = self.add_tool('meta.beta_diversity.anosim')
        self.anosim_box = self.add_tool('meta.beta_diversity.anosim_box')

    def check_option(self):
        """
        参数检查
        """
        if not self.option("otu_table"):
            raise OptionError("必须输入数据表")
        if not self.option("group_table"):
            raise OptionError("必须输入分组表")
        if not self.option("method"):
            raise OptionError("必须输入距离算法")
        elif self.option("method") not in ["abund_jaccard","binary_chisq",
                    "binary_chord","binary_euclidean","binary_hamming","binary_jaccard","binary_lennon","binary_pearson","binary_sorensen_dice",
                    "bray_curtis","bray_curtis_faith","bray_curtis_magurran","canberra","chisq","chord",
                    "euclidean","gower","hellinger","kulczynski","manhattan","morisita_horn","pearson","soergel","spearman_approx","specprof"]:
            raise OptionError("{} 不在可选的距离算法内".format(self.option("method")))

    def run_tools(self):
        self.group=self.check_group()
        self.otu_table = self.option('otu_table').prop['path']
        self.matrix.set_options({'method' :self.option('method'),
                                         'otutable': self.otu_table})
        self.matrix.on('end', self.anosim_run)
        self.matrix.run()

    def anosim_run(self):
        self.anosim.set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            "group" :self.group,
            "permutations": 999
        })
        self.anosim.on('end',self.anosim_box_run)
        self.anosim.run()

    def anosim_box_run(self):
        self.anosim_box.set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            "group":self.group,
            "permutations": 999
        })
        self.anosim_box.on('end',self.set_output)
        self.anosim_box.run()

    def set_output(self):
        self.logger.info('strat set_output as {}'.format(
            self.__class__.__name__))
        try:
            with open(os.path.join(self.anosim.output_dir,"format_results.xls"),"r") as old_re,open(os.path.join(self.output_dir,"anosim_results.txt"),"w") as new_re:
                line = old_re.readline()
                title = line.rstrip().split('\t')
                title_out = []
                for i in title:
                    if "p-value" in i:
                        title_out.append("Pvalue")
                    elif "permutation" in i :
                        title_out.append("Permutation_number")
                    else:
                        title_out.append(i.capitalize())
                new_re.write("\t".join(title_out))
                new_re.write('\n')
                line1 = old_re.readline()
                fd = line1.rstrip().split('\t')
                fd[0] = fd[0].upper()
                new_re.write("\t".join(fd))

            # os.link(os.path.join(self.anosim.output_dir,"format_results.xls")
            #         ,os.path.join(self.output_dir,"anosim_results.txt"))

        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))
        self.set_db()
    
    def set_db(self):
        self.logger.info("开始导表")
        api_anosim = self.api.api("tool_lab.anosim")
        api_anosim.add_box_detail(self.option("main_id")
                    ,os.path.join(self.anosim_box.output_dir,"box_data.xls"),self.option("method"))
        api_anosim.add_table_detail(self.option("main_id")
                    ,os.path.join(self.output_dir,"anosim_results.txt"))
        self.logger.info("导表结束")
        self.end()

    def run(self):
        """
        运行
        """
        self.run_tools()
        super(AnosimWorkflow, self).run()

    def check_group(self):
        group_path = os.path.join(self.work_dir,"group_tem.xls")
        with open(group_path,"w") as new_group:
            with open(self.option("group_table").prop['path'],"r") as group:
                line = group.readline()
                if line[1] == "#":
                    new_group.write(line)
                else:
                    new_group.write("#")
                    new_group.write(line)
                title = line.rstrip().split('\t')[1:]
                group_num = {}
                group_list = {}
                for i in title:
                    group_num[i] = []
                    group_list[i] = []
                while 1:
                    line = group.readline()
                    if not line:
                        break
                    new_group.write(line)
                    fd = line.rstrip().split('\t')
                    for i in range(len(title)):
                        group_list[title[i]].append(fd[i+1])
                        if fd[i+1] not in group_num[title[i]]:
                            group_num[title[i]].append(fd[i+1])
                for i in title:
                    if len(group_num[i]) < 2:
                        self.set_error("分组列表的组{}的分组数小于2".format(i))
                    for j in group_num[i]:
                        if group_list[i].count(j) < 3:
                            self.set_error("分组列表的组{}的{}组小于3".format(i,j))
        return group_path

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AnosimWorkflow, self).end()

if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    data = {
        'name': 'test_anosim',
        'id': 'anosim_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        'options': {
        "otu_table": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/anosim/otu.txt",
        "group_table":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/anosim/group.txt",
        "method":"bray_curtis",
        "main_id" : "5e9e6a6017b2bf2049a81be6"
        }
    }
    wsheet = Sheet(data=data)
    wf = AnosimWorkflow(wsheet)
    wf.run()