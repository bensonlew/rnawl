# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os,glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import json
from collections import OrderedDict

class ExtractDetailSeqAgent(Agent):
    """
    This script is used to extracts mapped/unmapped reads from input bam/sam file.
    """
    def __init__(self, parent):
        super(ExtractDetailSeqAgent, self).__init__(parent)
        options = [
            {"name": "seq_detail_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {'name': 'seq_db_mirna_novol', 'type': 'string', 'default': None},
            {'name': 'seq_db_mirna_known', 'type': 'string', 'default': None},
            {"name": "geneset_extract", "type": "infile", "format": 'ref_rna_v2.common'},
            {"name": "level", "type": "string", "default": None},
            {"name": "extract_info", "type": "string", "default": None},
            {"name": 't2g_file', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.step.add_steps("extract")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.extract.start()
        self.step.update()

    def stepfinish(self):
        self.step.extract.finish()
        self.step.update()

    def check_options(self):
        # if not self.option('input_file').is_set:
        #     raise OptionError('SAM/BAM文件必须输入')
        # if self.option('seq_type') not in ["PE", "SE"]:
        #     raise OptionError('测序类型参数输入有误')
        if self.option('level').lower() not in ["g", "m"]:
            raise OptionError('暂不支持提取该类型: {}'.format(self.option('level')))
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        super(ExtractDetailSeqAgent, self).end()


class ExtractDetailSeqTool(Tool):
    def __init__(self, config):
        super(ExtractDetailSeqTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'seq_extract1': os.path.join(self.config.PACKAGE_DIR, 'small_rna_v2/seq_extract.py'),
            'seq_extract2': os.path.join(self.config.PACKAGE_DIR, 'small_rna_v2/seq_extract2.py')
        }
        self.extract_info = ""



    def run(self):
        """
        运行
        :return:
        """
        super(ExtractDetailSeqTool, self).run()
        self.run_tool()
        self.set_output()
        self.end()

    def run_tool(self):
        self.extract_info = json.loads(self.option("extract_info"), object_pairs_hook=OrderedDict)
        extract_types = ",".join(self.extract_info["seq_type"])
        if self.option('seq_db_mirna_novol') and self.option('seq_db_mirna_known'):
            seq_db = os.path.join(self.work_dir, 'seq.fasta')
            os.system('cat {} {} > {}'.format(self.option('seq_db_mirna_novol'), self.option('seq_db_mirna_known'), seq_db))
            cmd = '{} {} '.format(self.program['python'], self.script['seq_extract1'])
            cmd += ' -s {}'.format(seq_db)
            cmd += ' -i {}'.format(self.option('geneset_extract').prop['path'])
            cmd += ' -l {}'.format(self.option('level'))
            cmd += ' -t {}'.format(extract_types)
            cmd += ' -o {}'.format(self.output_dir)
            cmd += ' -g {}'.format(self.option('t2g_file'))
            cmd_name = 'extract_{}'.format(self.option('level').lower())
            runcmd(self, cmd_name, cmd)
        else:
            cmd = '{} {} '.format(self.program['python'], self.script['seq_extract2'])
            cmd += ' -r {}'.format(self.option('seq_detail_dir').prop['path'])
            cmd += ' -i {}'.format(self.option('geneset_extract').prop['path'])
            cmd += ' -l {}'.format(self.option('level'))
            cmd += ' -t {}'.format(extract_types)
            cmd += ' -o {}'.format(self.output_dir)
            cmd_name = 'extract_{}'.format(self.option('level').lower())
            runcmd(self, cmd_name, cmd)
    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        # try:
        #     output_files = glob.glob(self.work_dir + "/*.fq")
        #     for file_path in output_files:
        #         os.link(file_path, os.path.join(self.output_dir, os.path.basename(file_path)))
        # except Exception as e:
        #     self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/medical_transcriptome/test_files'
        data = {
            "id": "Extract_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna_v2.extract_detail_seq",
            "instant": False,
            "options":
                dict(
                seq_db ="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/test_data/refrna_seqs.db",
                geneset_extract = "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/test_data/geneset_list",
                level = "G",
                extract_info = json.dumps({"seq_type":["transcript", "cds", "pep"]}),
                t2g_file = '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/test_data/tran2gene.txt'
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
