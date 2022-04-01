# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class BwaAgent(Agent):
    def __init__(self, parent):
        super(BwaAgent, self).__init__(parent)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {'name': 'fastq_s', 'type': 'infile', 'format': 'sequence.fastq'},  # 输入文件SE序列
            {'name': 'fastq_r', 'type': 'infile', 'format': 'sequence.fastq'},  # 输入文件PE的右端序列
            {'name': 'fastq_l', 'type': 'infile', 'format': 'sequence.fastq'},  # PE的左端序列
            {'name': 'sample_name', 'type': 'string', 'default': None},
            {'name': 'fq_type', 'type': 'string', 'default': None},
            {"name": "seq_file", "type": "infile", "format": "ref_rna_v2.fasta"},   # target sequence in fasta format
            {'name': 'out', 'type': 'string'},
        ]
        self.add_option(options)
        self.step.add_steps("bwa")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.bwa.start()
        self.step.update()

    def step_finish(self):
        self.step.bwa.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fastq_l").is_set and not self.option("fastq_r").is_set and not self.option("fastq_s").is_set:
            raise OptionError("必须设置输入原始reads文件")
        if self.option("fastq_l").is_set and not self.option("fastq_r").is_set:
            raise OptionError("必须设置输入右端原始reads文件")
        if self.option("fastq_r").is_set and not self.option("fastq_l").is_set:
            raise OptionError("必须设置输入左端原始reads文件")
        if not self.option("seq_file").is_set:
            raise OptionError("必须设置输入目的基因序列")
        return True

    def set_resource(self):
        self._cpu = 8
        self._memory = "20G"

    def end(self):
        super(BwaAgent, self).end()


class BwaTool(Tool):
    def __init__(self, config):
        super(BwaTool, self).__init__(config)
        self._version = "v1.0"
        self.program = {
            'samtools': 'bioinfo/align/samtools-1.8/samtools',
            'minimap2': 'bioinfo/align/minimap2/minimap2-2.17_x64-linux/minimap2',
            'bwa': 'bioinfo/align/bwa-0.7.17/bwa',
        }
        self.file = {
            'outfile':  os.path.join(self.option('out'), self.option('sample_name') + '_map.sam')
        }

    def run(self):
        super(BwaTool, self).run()
        self.bwa_index()
        self.bwa_aln()
        self.set_output()
        self.end()

    def bwa_index(self):
        self.logger.info(self.option('sample_name'))
        cmd = "{} index {}".format(self.program['bwa'], self.option('seq_file').prop['path'])
        self.logger.info("使用bwa对fasta文件建立索引")
        command = self.add_command("bwa_index", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用bwa对fasta文件建立索引完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用bwa对fasta文件建立索引出错！")

    def bwa_aln(self):
        self.logger.info(self.option('sample_name'))
        cmd = "{} mem -t 8 {} ".format(os.path.join(self.config.SOFTWARE_DIR, self.program['bwa']),
                                       self.option("seq_file").prop["path"])
        if self.option('fastq_s').is_set:
            cmd += '{} '.format(self.option('fastq_s').prop['path'])
        else:
            cmd += '{} {} '.format(self.option('fastq_l').prop['path'],
                                   self.option('fastq_r').prop['path'])
        cmd += '> {}'.format(self.file['outfile'])
        self.logger.info("使用bwa进行序列比对")
        command = self.add_command("bwa_align", cmd, ignore_error=True, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用bwa进行序列比对完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用bwa进行序列比对引出错！")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        pass


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "target_depth_bwa_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.target_depth.bwa",
            "instant": False,
            "options": dict(
                sample_name='HY_1',
                fastq_l='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/308_TagGeneInRawSeq/script/bwa/HY_1.1.fq',
                fastq_r='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/308_TagGeneInRawSeq/script/bwa/HY_1.2.fq',
                seq_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/308_TagGeneInRawSeq/script/bwa/NDM_5.fasta',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)