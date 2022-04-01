# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError
import pandas as pd


class SelectiveSweepWorkflow(Workflow):
    """
    carry out selective sweep analysis  using VCFtools
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SelectiveSweepWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "group_file", "type": "infile", "format": 'ref_rna_v2.common'},  # 分组文件,输出文件以分组文件名称.分割命名
            {"name": "max_missing", "type": "string", 'default': '30%'},  # --max-missing
            {"name": "maf", "type": "string", 'default': '0.05'},  # --maf 0.05
            {"name": "win_size", "type": "string", "default": '10'},  # --window-pi,窗口大小
            {"name": "win_step", "type": "string", "default": '2'},  # --window-pi-step,窗口步长
            {'name': 'threshold', 'type': 'string', 'default': '0.05'},
            {'name': 'advanced_params', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.filter = self.add_tool("tool_lab.selective_sweep.vcftools_filter")
        self.stats = self.add_tool("tool_lab.selective_sweep.vcftools_stat")
        self.graph = self.add_tool("tool_lab.selective_sweep.sweep_graph")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_filter()
        super(SelectiveSweepWorkflow, self).run()

    def check_options(self):
        if not self.option("vcf_path").is_set:
            raise OptionError("必须设置输入VCF文件")
        if not self.option("group_file").is_set:
            raise OptionError("必须设置输入分组方案文件")
        if self.option("vcf_path").is_set and self.option("group_file").is_set:
            group_df = pd.read_table(self.option("group_file").prop['path'], header=None)
            sample = list(group_df.iloc[:, 0])
            with open(self.option("vcf_path").prop['path'], 'r') as vcf:
                for line in vcf:
                    if line.startswith('#CHROM'):
                        header = line.strip().split('\t')[8:]
                        break
                    else:
                        continue
            total_len = len(sample) + len(header)
            diff = set(sample).symmetric_difference(header)
            if total_len == len(diff):
                raise OptionError("分组group文件与vcf文件不对应，请检查数据！")
        return True

    def run_filter(self):
        missing = self.option('max_missing').split('%')[0]
        opts = {
            'vcf_path': self.option('vcf_path').prop['path'],
            'max_missing': float(missing)/100,
            'maf': float(self.option('maf')),
        }
        self.filter.set_options(opts)
        self.filter.on('end', self.run_stats)
        self.filter.run()

    def run_stats(self):
        opts = {
            'vcf_file': self.filter.option('filtered_vcf').prop['path'],
            'group_file': self.option('group_file').prop['path'],
            'window_size': int(self.option('win_size'))*10**6,
            'window_step': int(self.option('win_step'))*10**6,
        }
        self.stats.set_options(opts)
        self.stats.on('end', self.run_graph)
        self.stats.run()

    def run_graph(self):
        opts = {
            'vcftools_dir': self.stats.output_dir,
            'threshold': float(self.option('threshold')),
        }
        self.graph.set_options(opts)
        self.graph.on('end', self.set_db)
        self.graph.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        sweep_api = self.api.api('tool_lab.selective_sweep')
        png = glob.glob(os.path.join(self.graph.output_dir, '*.png'))[0]
        png = os.path.join(self._sheet.output, os.path.basename(png))
        pdf = glob.glob(os.path.join(self.graph.output_dir, '*.pdf'))[0]
        pdf = os.path.join(self._sheet.output, os.path.basename(pdf))
        result = os.path.join(self.graph.output_dir, 'result_table.xls')
        sweep_api.add_selective_sweep(result=result, png=png, pdf=pdf, main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.graph.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "基于VCF进行选择性消除分析结果",0],
            ['./result_table.xls', 'xls', '选择性消除差异表', 0],
            ['./*png', 'png', '选择性消除分布图', 0],
            ['./*pdf', 'pdf', '选择性消除分布图', 0],
            ['./*select', '', '选择性消除结果文件', 0],
            ['./*detail', '', '选择性消除结果文件', 0],
        ])
        super(SelectiveSweepWorkflow, self).end()


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
            vcf_path='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/selective_sweep/final.vcf',
            group_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/selective_sweep/group.txt',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.selective_sweep",
            main_table_name="sg_selective_sweep",
            task_id="selective_sweep",
            project_sn="selective_sweep",
            submit_location="selective_sweep"
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