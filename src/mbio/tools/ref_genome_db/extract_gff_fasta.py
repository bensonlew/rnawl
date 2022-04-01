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


class ExtractGffFastaAgent(Agent):
    """
    有参转录组拼接
    version v1.0.1
    author: liubinxu
    last_modify: 2019.01.09
    """
    def __init__(self, parent):
        super(ExtractGffFastaAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "infile", "format": "ref_genome_db.fasta"},  # 参考基因文件
            {"name": "gtf", "type": "infile", "format": "ref_genome_db.gtf"},  # 参考基因的注释文件
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
        if not self.option('ref_fa'):
            raise OptionError('必须输入参考序列ref.fa')
        if not self.option('gtf'):
            raise OptionError('必须输入参考序列gtf')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(ExtractGffFastaAgent, self).end()


class ExtractGffFastaTool(Tool):
    def __init__(self, config):
        super(ExtractGffFastaTool, self).__init__(config)
        self._version = "v1.0.1"
        self.cufflinks_path = '/bioinfo/rna/cufflinks-2.2.1/'
        self.bioawk_path =  Config().SOFTWARE_DIR + '/bioinfo/seq/bioawk/'
        self.bedtools_path = Config().SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/'
        self.script_path = '/bioinfo/rna/scripts/'

    def run(self):
        """
        运行
        :return:
        """
        super(ExtractGffFastaTool, self).run()
        self.run_gtf_to_fa()
        self.set_output()
        self.end()


    def run_gtf_to_fa(self):
        """
        运行gtf_to_fasta，转录本gtf文件转fa文件
        """
        cmd = self.cufflinks_path + "gffread %s -g %s -w %s" % (
            self.option("gtf").prop['path'],
            self.option('ref_fa').prop['path'],
            "transcript.fa")
        self.logger.info('运行gtf_to_transcripts，形成fasta文件')
        command = self.add_command("gtf_to_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!")

        cmd = self.cufflinks_path + "gffread %s -g %s -x %s -y %s" % (
            self.option("gtf").prop['path'],
            self.option('ref_fa').prop['path'],
            "cds.fa",
            "pep.fa"
        )
        self.logger.info('运行gtf_to_cds，形成cds pep文件')
        command = self.add_command("gtf_to_fa2cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_cds运行完成")
        else:
            self.set_error("gtf_to_cds运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            all_files = ['transcript.fa', 'cds.fa', 'pep.fa']
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
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "extract_gff_fasta" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db.extract_gff_fasta",
            "instant": True,
            "options": dict(
                ref_fa=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                gtf=test_dir + "/" + "gtf/Saccharomyces_cerevisiae.R64-1-1.39.gtf",

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
