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


class BedToBiomartAgent(Agent):
    """
    根据ncbi基因组提取信息  enterz, gene name, gene description
    version v1.0.1
    author: liubinxu
    last_modify: 2019.01.09
    bed: 基因结构注释bed格式
    gene2tran: 基因转录本对应关系文件
    tran2des: 转录本功能描述文件
    tran2name: 转录本和基因名对应关系 提供输出为type2 否则为type3
    source: 基因名来源
    """
    def __init__(self, parent):
        super(BedToBiomartAgent, self).__init__(parent)
        options = [
            {'name': 'bed', 'type': 'string', 'default': None},
            {'name': 'tran2des', 'type': 'string', 'default': None},
            {'name': 'tran2name', 'type': 'string', 'default': None},
            {'name': 'type', 'type': 'string', 'default': None},
            {'name': 'gene2tran', 'type': 'string', 'default': None},
            {'name': 'gene2tran2pep', 'type': 'string', 'default': None},
            {'name': 'source', 'type': 'string', 'default': "unknown"},
            {'name': 'out_pre', 'type': 'string', 'default': "biomart"},
            {"name": "cpu", "type": "int", "default": 1},
        ]
        self.add_option(options)
        self.step.add_steps("bed_to_biomart")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.bed_to_biomart.start()

    def stepfinish(self):
        self.step.bed_to_biomart.finish()


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('bed'):
            raise OptionError('必须输入bed')
        if not self.option('gene2tran') and not self.option('gene2tran2pep'):
            raise OptionError('必须输入gene2tran 或 gene2tran2pep')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(BedToBiomartAgent, self).end()


class BedToBiomartTool(Tool):
    def __init__(self, config):
        super(BedToBiomartTool, self).__init__(config)
        self._version = "v1.0.1"
        self.bed_to_biomart = self.config.PACKAGE_DIR + "/ref_genome_db/bed2biomart2.pl"
        self.bed_to_biomart3 = self.config.PACKAGE_DIR + "/ref_genome_db/bed2biomart3.pl"
        self.perl = "miniconda2/bin/perl"


    def run(self):
        """
        运行
        :return:
        """
        super(BedToBiomartTool, self).run()
        self.run_bed_to_biomart()
        self.set_output()
        self.end()

    def run_bed_to_biomart(self):
        '''
        提取NCBI gff相关信息
        '''
        if self.option("gene2tran"):
            cmd = "{} {} ".format(self.perl, self.bed_to_biomart)
            cmd += " {} {}".format("--m", self.option("gene2tran"))
        elif self.option("gene2tran2pep"):
            cmd = "{} {} ".format(self.perl, self.bed_to_biomart3)
            cmd += " {} {}".format("--m", self.option("gene2tran2pep"))

        cmd += " {} {}".format("--i", self.option("bed"))
        cmd += " {} {}".format("--o", self.option("out_pre"))


        if self.option("tran2des"):
            cmd += " {} {}".format("--r", self.option("tran2des"))
        else:
            pass

        if self.option("tran2name"):
            cmd += " {} {}".format("--n", self.option("tran2name"))
            cmd += " {} {}".format("--nf", self.option("source"))
        else:
            pass

        self.logger.info('运行bed2biomart')
        command = self.add_command("bed2biomart", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行bed2biomart成功")
        else:
            self.set_error("运行bed2biomart出错")


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            all_files = [self.option("out_pre")]
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
            "id": "bed_to_biomart" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db.bed_to_biomart",
            "instant": True,
            "options": dict(
                bed="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Glarea_lozoyensis/NCBI/gtf/GCF_000409485.1_GLAREA.bed",
                gene2tran2pep="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Glarea_lozoyensis/NCBI/biomart/gene2trans2pep",
                tran2des = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Glarea_lozoyensis/NCBI/GCF_000409485.1_GLAREA_genomic.gff.tran2des",

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
