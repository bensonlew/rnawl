# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.ref_rna.filter_gtf import FilterGtf
import re
import unittest


class NcbiGffAgent(Agent):
    """
    根据ncbi基因组提取信息  enterz, gene name, gene description
    version v1.0.1
    author: liubinxu
    last_modify: 2019.01.09
    """
    def __init__(self, parent):
        super(NcbiGffAgent, self).__init__(parent)
        options = [
            {"name": "gff", "type": "infile", "format": "ref_genome_db.gtf"},  # 参考基因的注释文件
            {'name': 'source', 'type': 'string', 'default': 'unknown'},
            {'name': 'output', 'type': 'string', 'default': 'ncbi'},
            {"name": "cpu", "type": "int", "default": 10},  #cufflinks软件所分配的cpu数量
        ]
        self.add_option(options)
        self.step.add_steps("cufflinks")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.cufflinks.start()
        self.step.update()

    def stepfinish(self):
        self.step.cufflinks.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('gff').is_set:
            raise OptionError('必须输入参考序列gtf 或 gff')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "30G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(NcbiGffAgent, self).end()


class NcbiGffTool(Tool):
    def __init__(self, config):
        super(NcbiGffTool, self).__init__(config)
        self._version = "v1.0.1"
        self.ncbi_gff = self.config.PACKAGE_DIR + "/ref_genome_db/ncbi_gff.py"
        self.python = "miniconda2/bin/python"


    def run(self):
        """
        运行
        :return:
        """
        super(NcbiGffTool, self).run()
        self.run_ncbi_gff()
        self.set_output()
        self.end()

    def run_ncbi_gff(self):
        '''
        提取NCBI gff相关信息
        '''
        cmd = "{} {} {} {}".format(self.python, self.ncbi_gff, self.option("gff").prop['path'], self.option("output"))
        self.logger.info('运行ncbi_gff')
        command = self.add_command("ncbi_gff", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("ncbi运行完成")
        else:
            self.set_error("ncbi运行出错!")


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            all_files = [self.option("output") + ".gene2enterz", self.option("output") +  '.tran2name', self.option("output") +  '.tran2des']
            for each in all_files:
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Trichoderma_virens/NCBI"
        data = {
            "id": "ncbi_gff" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db.ncbi_gff",
            "instant": True,
            "options": dict(
                gff=test_dir + "/" + "GCF_000170995.1_TRIVI_v2.0_genomic.gff"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
