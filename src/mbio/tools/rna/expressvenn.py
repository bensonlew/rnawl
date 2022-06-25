#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
# import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.graph.venn_table import venn_graph
import unittest

class ExpressvennAgent(Agent):
    """
    version 1.0
    author: khl
    last_modify: 2017.07.10
    需要R软件
    """
    def __init__(self, parent):
        super(ExpressvennAgent, self).__init__(parent)
        options = [
            {"name": "express_matrix", "type": "infile", "format": "rna.express_matrix"}, #表达量矩阵表
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table, meta.otu.group_table, sample.group_table"},  # 输入的group表格
            {"name": "threshold","type":"float", "default":1}  #对表达量的值进行过滤
        ]
        self.add_option(options)
        self.step.add_steps('express_venn')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.express_venn.start()
        self.step.update()

    def step_end(self):
        self.step.express_venn.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        :return:
        """
        if not self.option("express_matrix"):
            raise OptionError("参数express_matrix不能为空")
        if not self.option("group_table").is_set:
            raise OptionError("参数group_table不能为空")
        # if self.option("group_table").format == 'meta.otu.otu_table':    # add by wzy 2017.6.23
        #     group_file = self.option("group_table").prop['path']

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 10
        self._memory = '20G'

    def end(self):  # add by wzy 20170608
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "venn图结果目录"],
        ])
        super(ExpressvennAgent, self).end()


class ExpressvennTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(ExpressvennTool, self).__init__(config)
        self.R_path = '/program/R-3.3.1/bin/'
        self.R_path_total = self.config.SOFTWARE_DIR + self.R_path + "R"
        self.venn_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/scripts/'
        self.python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self._version = 1.0

    def _create_venn_table(self):
        """
        调用脚本venn_table.py,输出venn表格
        """
        group_file = self.option("group_table").prop['path']
        express_matrix = self.option("express_matrix").prop['path']
        num_group_file = sum(1 for line in open(group_file))
        if num_group_file == 1:
            tmp_group_path = os.path.join(self.work_dir,"new_group.txt")
            with open(express_matrix,'r+') as f1,open(tmp_group_path,'w+') as f2:
                header = f1.readline().strip().split("\t")
                f2.write("#control\tgroup\n")
                for line in header:
                    f2.write(line+"\t"+line+"\n")

        num_lines = sum(1 for line in open(express_matrix))
        if num_lines < 11:
            self.set_error("输入文件的行数小于10个！请更换输入文件！")
            raise Exception("输入文件的行数小于10个！请更换输入文件！")
        if num_group_file >1:
            venn_cmd = '%spython %sexpressvenntable.py -i %s -g %s -f %s -o %s -out %s' % (self.python_path, self.venn_path, express_matrix, group_file,self.option("threshold"),self.R_path_total,                                                                                 self.work_dir)
        # add by qiuping, for denovo_rna venn, 20160728
        else:
            # self.option('group_table').sub_group(self.work_dir + '/venn_group', self.option("group_table").prop['group_scheme'][0])
            venn_cmd = '%spython %sexpressvenntable.py -i %s -g %s -f %s -o %s -out %s' % (self.python_path, self.venn_path, express_matrix, self.work_dir + '/new_group.txt',self.option("threshold"),
                                                                                   self.R_path_total,self.work_dir)
        # add end
        self.logger.info(venn_cmd)
        os.system(venn_cmd)
        cmd = self.R_path + 'Rscript cmd.r'
        # print cmd
        self.logger.info("开始运行venn_table")
        command = self.add_command("venn_table", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行venn_table完成")
        else:
            self.set_error("运行venn_table运行出错!")
            raise Exception("运行venn_table运行出错，请检查输入的otu表和group表是否正确")
        # 统计各组所有otu/物种名 add by qindanhua
        venn_graph(express_matrix, group_file, "venn_graph.xls")
        self.set_output()

    def set_output(self):
        """
        将结果文件链接至output
        """
        self.logger.info("set out put")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        os.link(self.work_dir + '/venn_table.xls', self.output_dir + '/venn_table.xls')
        if os.path.exists(self.work_dir + "/venn_graph.xls"):
            os.link(self.work_dir + '/venn_graph.xls', self.output_dir + '/venn_graph.xls')
        # self.option('venn_table.xls').set_path(self.output_dir+'/venn_table.xls')
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(ExpressvennTool, self).run()
        self._create_venn_table()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/rna/??'
        data = {
            "id": "ExpVenn" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "rna.expressvenn",
            "instant": True,
            "options": {
                "express_matrix": "/mnt/ilustre/users/sanger-dev/workspace/20190701/Refrna_tsg_34668/Quant" + '/' + 'transcript.tpm.matrix',
                "group_table": "/mnt/ilustre/users/sanger-dev/workspace/20190701/Refrna_tsg_34668/remote_input/group_table" + '/' + 'default_group.txt',
                "threshold": 1,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()