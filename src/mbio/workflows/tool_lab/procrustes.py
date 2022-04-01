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


class ProcrustesWorkflow(Workflow):
    """
    Procrustes Analysis
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProcrustesWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "df1", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "df2", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {'name': 'group_names', 'type': 'string', 'default': ''},
            {'name': 'distance', 'type': 'string', 'default': 'euclidean'},
            {'name': 'method', 'type': 'string', 'default': 'pcoa'},
            {'name': 'df1_log', 'type': 'string', 'default': 'log2'},
            {'name': 'df2_log', 'type': 'string', 'default': 'log2'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.procrustes = self.add_module("tool_lab.procrustes")
        self.ref_select = self.add_tool("tool_lab.procrustes.table_select")
        self.query_select = self.add_tool("tool_lab.procrustes.table_select")
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("df1").is_set:
            raise OptionError("必须设置输入数据表一文件")
        if not self.option("df2").is_set:
            raise OptionError("必须设置输入数据表二文件")
        if not self.option("group_file").is_set:
            raise OptionError("必须设置输入样本分组文件")
        return True

    def run(self):
        self.on_rely([self.ref_select, self.query_select], self.run_procrustes)
        self.procrustes.on('end', self.set_db)
        self.select_ref()
        self.select_query()
        super(ProcrustesWorkflow, self).run()

    def select_ref(self):
        options = {
            "origin_table": self.option('df1').prop['path'],
            "group": self.option("group_file").prop['path'],
            'group_names': self.option('group_names'),
            "log": self.option('df1_log'),
            "method": self.option('method'),
        }
        self.ref_select.set_options(options)
        self.ref_select.run()

    def select_query(self):
        options = {
            "origin_table": self.option('df2').prop['path'],
            "group": self.option("group_file").prop['path'],
            'group_names': self.option('group_names'),
            "log": self.option('df2_log'),
            "method": self.option('method'),
        }
        self.query_select.set_options(options)
        self.query_select.run()

    def run_procrustes(self):
        options = dict(
            ref_table=self.ref_select.option("select_table").prop['path'],
            dist=self.option('distance'),
            query_table=self.query_select.option("select_table").prop['path'],
            method=self.option('method'),
            group_file=self.option('group_file').prop['path'],
        )
        self.procrustes.set_options(options)
        self.procrustes.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        procrustes_api = self.api.api('tool_lab.procrustes')
        summary = os.path.join(self.procrustes.output_dir, 'procrustes_results.txt')
        ref = glob.glob(os.path.join(self.procrustes.output_dir, '*transformed_reference.txt'))[0]
        query = glob.glob(os.path.join(self.procrustes.output_dir, '*transformed_q*.txt'))[0]
        procrustes_api.add_procrustes_main(summary=summary, ref=ref, query=query,
                                           group=self.option('group_file').prop['path'], main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.procrustes.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "普氏分析结果文件",0],
            [r'*.txt', 'txt', '普氏分析结果文件', 0],
        ])
        super(ProcrustesWorkflow, self).end()


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
            df1='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/procrustes/IAH_vs_Control.exp.nodiff.txt',
            df2='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/procrustes/IAH_vs_Control.Genus.otu.txt',
            group_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/procrustes/IAH_vs_Control.group.txt',
            method='pcoa',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.procrustes",
            main_table_name="sg_procrustes",
            task_id="procrustes",
            project_sn="procrustes",
            submit_location="procrustes"
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