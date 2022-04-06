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
            {"name": "geneset_extract", "type": "infile", "format": 'ref_rna_v2.common'},
            {"name": "level", "type": "string", "default": None},
            {"name": "extract_info", "type": "string", "default": None},
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
        if self.option('level').lower() not in ["g", "t"]:
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
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'seq_extract': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/seq_extract_detail_new.py')
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
        if self.extract_info["seq_type"]:
            extract_types = ",".join(self.extract_info["seq_type"])
        else:
            extract_types = None
        cmd = '{} {} '.format(self.program['python'], self.script['seq_extract'])
        cmd += ' -r {}'.format(self.option('seq_detail_dir').prop['path'])
        cmd += ' -i {}'.format(self.option('geneset_extract').prop['path'])
        cmd += ' -l {}'.format(self.option('level'))
        if extract_types:
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
            "id": "ExpPca" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "medical_transcriptome.extract_detail_seq",
            "instant": False,
            "options": dict(
                seq_detail_dir ="/mnt/ilustre/users/sanger-dev/workspace/20201104/DownloadDetailSeq_gbnijnfcuslii6rieegpbdgpt5_9045_87/remote_input/seq_detail_dir/SequenceDetail/",
                geneset_extract = "/mnt/ilustre/users/sanger-dev/workspace/20201104/DownloadDetailSeq_gbnijnfcuslii6rieegpbdgpt5_9045_87/geneset_extract_gene.list",
                level = "G",
                extract_info = json.dumps({"seq_type":["transcript", "cds", "pep"]}).replace('"', '\\"'),
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
