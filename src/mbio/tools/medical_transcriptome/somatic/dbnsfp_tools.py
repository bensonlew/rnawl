# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# modified 2020.10.26

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import os
import subprocess
from biocluster.config import Config

class DbnsfpToolsAgent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(DbnsfpToolsAgent, self).__init__(parent)
        options = [
            {"name": "combine_vcf", "type": "infile", "format": "ref_rna_v2.common"},  # 输入的bam
            {"name": "line_num", "type": "int", "default": 10000},  # 输入格式  bam/cram 20191231
            {"name": "fa_file", "type": "infile", "format": "ref_rna_v2.common"},  # 参考基因组文件
            {"name": "file_format", "type": "string", "default": "bam"},  # 输入格式  bam/cram 20191231
            {"name": "name", "type": "string"},  # 生成文件名字，测试文件中名字为DE1_10.g.vcf，其中.g.vcf为固定。

        ]
        self.add_option(options)
        self._memory_increase_step = 200

    def check_options(self):
        # if not self.option("bam_file"):
        #     raise OptionError("请设置bam路径")
        # if not self.option("fa_file"):
        #     raise OptionError("请设置ref.fa路径")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(DbnsfpToolsAgent, self).end()


class DbnsfpToolsTool(Tool):
    def __init__(self, config):
        super(DbnsfpToolsTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.gatk_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v2/gatk-4.1.4.1/gatk-4.1.4.1/"
        # self.gatk_path = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/gatk_4.1.4test/software/gatk-4.1.4.1/"
        self.java_path = self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/bin/java"
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.file = {
            'identity_file': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/stem/identity_file')
        }
        self.script = {
            "search_dbnsfp": "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/somatic_mutation/software/dbNSFP4/search_dbNSFP41a.jar",
            'bash': os.path.join(self.work_dir, 'dbfbsfp.sh')
        }
        self.conf = Config()

    def search_vcfnew(self):
        hostname = {'sanger-dev': '192.168.12.101', 'sanger': '10.2.0.110', 'isanger': '10.2.0.115'}[
            self.conf.wpm_user]
        cmd = '{} -mx4096M -jar {} -v hg38 -i {} -o {} -p '.format(
            self.java_path, self.script['search_dbnsfp'], self.option("combine_vcf").prop["path"],
            os.path.join(self.output_dir, "result")
        )
        flag = True
        n = {'sanger-dev': 70, 'sanger': 60, 'isanger': 50}[self.conf.wpm_user]
        while flag and n:
            n -= 1
            s = '''ssh -i {} {}@{} "cd {}; export DISPLAY=localhost:{}.0; {}; exit"\n'''.format(
                self.file['identity_file'], self.conf.wpm_user, hostname, self.work_dir, n, cmd
            )
            open(self.script['bash'], 'w').write(s)
            spc = subprocess.Popen(
                'sh {}'.format(self.script['bash']), shell=True,
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            ret = spc.wait()
            self.logger.debug('stdout of stem:\n{}'.format(spc.stdout.read()))
            self.logger.debug('stderr of stem:\n{}'.format(spc.stderr.read()))
            if ret:
                self.logger.debug('fail to run stem at localhost:{}'.format(n))
            else:
                self.logger.info('succeed in running at localhost:{}'.format(n))
                flag = False
                self.logger.info('stop iteration after succeed in running dbnsfp')
                self.logger.info('#' * 64)


    def search_vcf(self):
        hostname = {'sanger-dev': '192.168.12.101', 'sanger': '10.2.0.110', 'isanger': '10.2.0.115'}[
            self.conf.wpm_user]
        cmd = '{} -mx10240M -jar {} -i {} -o {}'.format(
            self.java_path, self.script['search_dbnsfp'], self.option("combine_vcf").prop["path"],
            os.path.join(self.output_dir, "result")
        )
        flag = True
        n = {'sanger-dev': 70, 'sanger': 60, 'isanger': 50}[self.conf.wpm_user]
        while flag and n:
            n -= 1
            s = '''ssh -i {} {}@{} "cd {}; export DISPLAY=localhost:{}.0; {}; exit"\n'''.format(
                self.file['identity_file'], self.conf.wpm_user, hostname, self.work_dir, n, cmd
            )
            open(self.script['bash'], 'w').write(s)
            spc = subprocess.Popen(
                'sh {}'.format(self.script['bash']), shell=True,
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            ret = spc.wait()
            self.logger.debug('stdout of stem:\n{}'.format(spc.stdout.read()))
            self.logger.debug('stderr of stem:\n{}'.format(spc.stderr.read()))
            if ret:
                self.logger.debug('fail to run stem at localhost:{}'.format(n))
            else:
                self.logger.info('succeed in running at localhost:{}'.format(n))
                flag = False
                self.logger.info('stop iteration after succeed in running dbnsfp')
                self.logger.info('#' * 64)


    def run(self):
        super(DbnsfpToolsTool, self).run()
        self.search_vcfnew()
        self.end()



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "dbnsfp" + str(random.randint(1, 10000))+"yyyy",
            "type": "tool",
            "name": "medical_transcriptome.somatic.dbnsfp_tools",
            "instant": False,
            "options": dict(
                combine_vcf="/mnt/ilustre/users/sanger-dev/workspace/20201105/MedicalTranscriptome_v6dpivfr84k4ooq2lmvh47dpfs/CallSnpIndelSentieon/VcfFilterGatk/output/final.vcf",
                # fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                # name="add_sort",
                # file_format="bam",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()