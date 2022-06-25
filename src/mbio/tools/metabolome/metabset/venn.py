# -*- coding: utf-8 -*-
# __author__ = 'liulinmeng'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import unittest
import glob

class VennAgent(Agent):
    """
    代谢组学代谢集Venn分析
    version 1.0
    author: liulinmeng
    last_modify: 2018.6.25
    需要R软件
    """
    def __init__(self, parent):
       super(VennAgent, self).__init__(parent)
       options = [
           {'name': 'list_file', 'type': 'infile', 'format': 'metabolome.mul_metabset', 'required': True}
       ]
       self.add_option(options)
       self.step.add_steps("venn_list")
       self.on('start', self.step_start)
       self.on('end', self.step_end)


    def step_start(self):
        self.step.venn_list.start()
        self.step.update()


    def step_end(self):
        self.step.venn_list.finish()
        self.step.update()


    def check_options(self):
        if not self.option("list_file"):
            raise OptionError("输入文件不能为空")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = '8G'


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "venn分析结果目录"],
        ])
        super(VennAgent, self).end()





class VennTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(VennTool, self).__init__(config)
        self.R_path = '/program/R-3.3.1/bin/'
        self.venn_path = self.config.PACKAGE_DIR + '/graph/scripts/venn_list.py'
        self.python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        #self.gcc = software_dir + '/gcc/5.1.0/bin'
        #self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        #self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self._version = 1.0



    def _create_venn_table(self):
        """
        调用脚本venn_list.py,输出venn分析结果表
        """
        list_file = self.option("list_file").prop['path']
        venn_cmd = '%spython %s -i %s  -o %s/cmd.r -n 7' % (
                self.python_path, self.venn_path, list_file, self.work_dir)
        self.logger.info(venn_cmd)
        os.system(venn_cmd)
        self.logger.info('运行venn_cmd,生成 R script')
        self.logger.info("开始运行venn_list")
        cmd = self.R_path + 'Rscript cmd.r'
        command = self.add_command("venn_list", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行venn_list完成")
        else:
            self.set_error("运行venn_list运行出错!")
            raise Exception("运行venn_list运行出错，请检查输入的文件是否正确")

    def set_output(self):
        self.logger.info("set output files")
        self.link_file("venn_table.xls","venn_table.xls")
        self.link_file("new_venn_graph.xls","metabset_detail.xls")

    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.work_dir, oldfile)
        newfile = os.path.join(self.output_dir,newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(oldfile, newfile)

    def run(self):
        super(VennTool, self).run()
        self._create_venn_table()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Venn" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "metabolome.metabset.venn",
            "instant": True,
            "options": dict(
                list_file="/mnt/ilustre/users/sanger-dev/sg-users/liulinmeng/metabolome/package/venn/test/metabset_file3",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

