# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
from biocluster.core.exceptions import OptionError

class SeqStructureWorkflow(Workflow):
    """
    hypergeometric test function for gene set enrichment analysis that are designed to accept user defined annotation
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SeqStructureWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "seq", "type": "string"},  # consists of gene id/gene name/transcript id, separating by comma
            {"name": "annotation", "type": "string"},   # an annotation file in GTF/GFF format
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.seq_structure")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(SeqStructureWorkflow, self).run()

    def check_options(self):
        if not self.option('seq'):
            raise OptionError('基因集文件必须输入')
        if not self.option('annotation'):
            raise OptionError("必须设置输入gtf/gff格式的注释文件。")
        return True

    def run_tool(self):
        file_format = os.path.basename(self.option('annotation')).split(".")[1]
        opts = {
            'seq': self.option('seq'),
        }
        if file_format.lower() == 'gff':
            opts.update({'annotation_gff': self.option('annotation')})
        if file_format.lower() == 'gtf':
            opts.update({'annotation_gtf': self.option('annotation')})
        self.tool.set_options(opts)
        self.tool.on('end', self.set_output)
        self.tool.run()

    def set_output(self):
        for each in os.listdir(self.tool.output_dir):
            if os.path.exists(os.path.join(self.output_dir, each)):
                os.remove(os.path.join(self.output_dir, each))
            os.link(os.path.join(self.tool.output_dir, each), os.path.join(self.output_dir, each))
        self.end()

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "基因结构结果文件",0],
        #     [r'*.pdf', 'pdf', '基因结构结果文件', 0],
        # ])
        super(SeqStructureWorkflow, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the workflow. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitoollabtest.py '
        cmd += 'post toollabpipeline '
        cmd += '-c {} '.format("client03")
        cmd += "-b http://bcl.tsg.com "
        cmd += "-n \"params;basis\" -d \"{"
        args = dict(
            seq="ARV1,AT1G01020.1,tryab",
            annotation="/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/Practice1/TAIR10.gtf",
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.seq_structure",
            main_table_name="seq_structure",
            task_id="seq_structure",
            project_sn="seq_structure",
            submit_location="seq_structure"
        )
        for arg in args:
            cmd += "\\\""
            cmd += arg
            cmd += "\\\":\\\""
            cmd += args[arg]
            cmd += "\\\","
        cmd = cmd.rstrip(",")
        cmd += "};{"
        for arg in config:
            cmd += "\\\""
            cmd += arg
            cmd += "\\\":\\\""
            cmd += config[arg]
            cmd += "\\\","
        cmd = cmd.rstrip(",")
        cmd += "}\""

        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
