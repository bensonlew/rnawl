# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest


class GtdbTkAgent(Agent):
    '''
    GTDB-TK的树种注释分类，数据库是GTDB数据库
    '''
    def __init__(self, parent):
        super(GtdbTkAgent, self).__init__(parent)
        options = [
            {'name': 'genome_dir', 'type': 'infile', 'format': 'sequence.fasta_dir'},  ## fasta格式的基因组文件，一个基因组一个文件

        ]
        self.add_option(options)

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option('genome_dir').is_set:
            raise OptionError("必须输入genome_dir的文件夹！")
        return True

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = '2'
        self._memory = '300G'

    def end(self):
        """
        运行结束
        :return:
        """
        super(GtdbTkAgent, self).end()

class GtdbTkTool(Tool):
    """
    GTDB-TK运行的命令
    """
    def __init__(self, config):
        super(GtdbTkTool, self).__init__(config)
        self.conda = self.config.SOFTWARE_DIR + "/program/miniconda3/etc/profile.d/conda.sh"
        self.database = self.config.SOFTWARE_DIR + "/database/GTDB/release95"
        self.shell = "program/sh"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "toolapps/gtdb_tk.sh")
        self.genomes = self.option("genome_dir").prop['path']


    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("开始运行tool!")
        super(GtdbTkTool, self).run()
        self.run_gtdbtk()
        self.set_output()
        self.end()

    def run_gtdbtk(self):
        """
        运行gtdbtk的sh，里面包含分析的三步
        :return:
        """
        self.logger.info("分析开始!")
        self.identify = self.work_dir + "/identify"
        self.align = self.work_dir + "/align"
        self.classify = self.work_dir + "/classify"
        tmp = self.work_dir + "/tmp"
        cmd = '{} {} {} {} {} {} {} {} {}'.format(self.shell, self.shell_path, self.conda, self.database, self.genomes, self.identify, self.align, self.classify, tmp)
        self.logger.info(cmd)
        command = self.add_command("run_gtdbtk", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('GTDB-tk软件运行成功！')
        else:
            self.set_error('GTDB-tk软件运行失败！')

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('开始设置结果文件目录')
        for i in ['all.bac120.tree','all.taxon.xls']:
            if os.path.exists(self.output_dir + "/" + i):
                os.remove(self.output_dir + "/" + i)
        if os.path.exists(self.work_dir + "/classify/classify/gtdbtk.bac120.summary.tsv"):
            os.link(self.work_dir + "/classify/classify/gtdbtk.bac120.summary.tsv", self.output_dir + "/all.taxon.xls")
        if os.path.exists(self.work_dir + "/classify/classify/gtdbtk.bac120.classify.tree"):
            os.link(self.work_dir + "/classify/classify/gtdbtk.bac120.classify.tree", self.output_dir + "/gtdbtk.bac120.classify.tree")


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "GTDB_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "toolapps.gtdb_tk",
            "options": dict(
                genome_dir="/mnt/ilustre/users/sanger-dev/home/gaohao/GTDB/genome3"

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()