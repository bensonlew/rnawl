# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from biocluster.core.exceptions import OptionError


class MirandaPicWorkflow(Workflow):
    """
    miRanda result visualisation
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MirandaPicWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "mirna", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "ref", "type": "infile", "format": "ref_rna_v2.fasta"},
            {'name': 'score', 'type': 'string', 'default': '160'},
            {'name': 'energy', 'type': 'string', 'default': '-20'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.miranda_pic")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(MirandaPicWorkflow, self).run()

    def check_options(self):
        if not self.option("mirna").is_set:
            raise OptionError("必须设置输入miRNA序列文件")
        if not self.option("ref").is_set:
            raise OptionError("必须设置输入靶基因序列文件")
        return True

    def run_tool(self):
        opts = {
            "mirna": self.option('mirna'),
            "ref": self.option('ref'),
            'score': self.option('score'),
            'energy': self.option('energy'),
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        miranda_pic_api = self.api.api('tool_lab.miranda_pic')
        pics = glob.glob(os.path.join(self.tool.output_dir, '*.png'))
        pic_s3 = list()
        for i in pics:
            f_name = os.path.basename(i)
            f_path = os.path.join(self._sheet.output, f_name)
            pic_s3.append(f_path)
        out_table = glob.glob(os.path.join(self.tool.output_dir, 'result_table.xls'))
        miranda_pic_api.add_miranda_pic(pics=pic_s3, out_table=out_table, main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "miRanda绘图结果文件",0],
            [r'*.png', 'png', 'miRanda结合位点图片', 0],
            ['miRanda_result.xls', 'xls', 'miRanda预测结果文件', 0],
            ['result_table.xls', 'xls', 'miRanda结合位点结果表', 0],
        ])
        super(MirandaPicWorkflow, self).end()


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
            mirna='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_042021/mirna_test.fa',
            ref='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_042021/ref_no_test.fa',
            score='160',
            energy='-20',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.miranda_pic",
            main_table_name="sg_miranda_pic",
            task_id="miranda_pic",
            project_sn="miranda_pic",
            submit_location="miranda_pic"
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