# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class Minimap2Agent(Agent):
    def __init__(self, parent):
        super(Minimap2Agent, self).__init__(parent)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {'name': 'fastq_s', 'type': 'infile', 'format': 'sequence.fastq'},  # 输入文件SE序列
            {'name': 'sample_name', 'type': 'string', 'default': None},
            {'name': 'fq_type', 'type': 'string', 'default': None},
            {"name": "seq_file", "type": "infile", "format": "ref_rna_v2.fasta"},  # target sequence in fasta format
            {'name': 'out', 'type': 'string'},
        ]
        self.add_option(options)
        self.step.add_steps("minimap2")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.minimap2.start()
        self.step.update()

    def step_finish(self):
        self.step.minimap2.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fastq_s").is_set:
            raise OptionError("必须设置输入原始reads文件夹")
        if not self.option("seq_file").is_set:
            raise OptionError("必须设置输入目的基因序列")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(Minimap2Agent, self).end()


class Minimap2Tool(Tool):
    def __init__(self, config):
        super(Minimap2Tool, self).__init__(config)
        self._version = "v1.0"
        self.program = {
            'samtools': 'miniconda2/bin/samtools',
            'minimap2_tsg': 'bioinfo/align/minimap2/minimap2-2.17_x64-linux/minimap2',
            'minimap2_others': 'bioinfo/align/minimap2/minimap2',
            'bwa': 'bioinfo/align/bwa-0.7.17/bwa',
        }
        self.file = {
            'outfile':  os.path.join(self.option('out'), self.option('sample_name') + '_map.sam')
        }

    def run(self):
        super(Minimap2Tool, self).run()
        if self.config.SOFTWARE_DIR == '/mnt/ilustre/users/sanger-dev/app':
            self.minimap_path = self.program['minimap2_tsg']
        else:
            self.minimap_path = self.program['minimap2_others']
        self.minimap_index()
        self.minimap_aln()
        self.set_output()
        self.end()

    def minimap_index(self):
        self.logger.info(self.option('sample_name'))
        ref_mmi = os.path.splitext(os.path.basename(self.option('seq_file').prop['path']))[0] + '.mmi'
        if os.path.isfile(os.path.join(self.work_dir, ref_mmi)):
            return
        cmd = "{} -d {} {}".format(self.minimap_path, os.path.join(self.work_dir, ref_mmi),
                                   self.option('seq_file').prop['path'])
        self.logger.info("使用minimap2对fasta文件建立索引")
        command = self.add_command("minimap2_index", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用minimap2对fasta文件建立索引完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用minimap2对fasta文件建立索引出错！")

    def minimap_aln(self):
        self.logger.info(self.option('sample_name'))
        cmd = "{} -ax map-pb ".format(os.path.join(self.config.SOFTWARE_DIR, self.minimap_path))
        cmd += "{} {} ".format(self.option("seq_file").prop["path"], self.option('fastq_s').prop['path'])
        cmd += '> {}'.format(self.file['outfile'])
        self.logger.info("使用minimap2进行序列比对")
        command = self.add_command("minimap2_align", cmd, ignore_error=True, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用minimap2进行序列比对完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用minimap2进行序列比对引出错！")

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
            "id": "target_depth_minimap2_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.target_depth.minimap2",
            "instant": False,
            "options": dict(
                sample_name='HY_1',
                fastq_s='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/308_TagGeneInRawSeq/script/minimap/HY_1.fastq',
                seq_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/308_TagGeneInRawSeq/script/minimap/NDM_5.fasta',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)