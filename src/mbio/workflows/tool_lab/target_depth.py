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


class TargetDepthWorkflow(Workflow):
    """
    Target genes in raw reads
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TargetDepthWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'ref_rna_v2.common_dir'},
            {"name": "seq_file", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "seq_str", "type": "string"},
            {'name': 'reads_type', "type": 'string'},
            {'name': 'seq_type', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.target_depth = self.add_module("tool_lab.target_depth")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_module()
        super(TargetDepthWorkflow, self).run()

    def check_options(self):
        if not self.option("fastq_dir").is_set:
            raise OptionError("必须设置输入fq序列文件夹")
        if not self.option('seq_file').is_set and not self.option('seq_str'):
            raise OptionError("必须设置输入目的基因序列")
        if not self.option('reads_type'):
            raise OptionError("必须设置输入原始reads类型")
        return True

    def generate_seq_file(self):
        seq_path = ''
        if self.option('seq_file').is_set:
            with open(self.option('seq_file').prop['path'], 'r') as f:
                seq = f.read()
            seq_str = ''.join(seq.split('\n')[1:])
            seq_num = seq.count('>')
            if seq_num == 1:
                seq_path = self.option('seq_file')
            else:
                self.set_error('检测到多条目的基因序列，请核查。')
        elif self.option('seq_str'):
            instr = self.option('seq_str').replace('&gt;', '>')
            seq_num = instr.count('>')
            seq_path = os.path.join(self.work_dir, 'target_seq.fa')
            if seq_num == 0:
                seq_str = instr
                with open(seq_path, 'w') as fa:
                    fa.write('>ref_seq\n')
                    seq = instr
                    while len(seq) > 60:
                        fa.write(seq[0:60] + '\n')
                        seq = seq[60:len(seq)]
                    fa.write(seq + '\n')
            elif seq_num == 1:
                seq_str = ''.join(instr.split('\n')[1:])
                with open(seq_path, 'w') as fa:
                    fa.write(instr)
            else:
                self.set_error('检测到多条目的基因序列，请核查。')
        self.seq_len = len(seq_str)
        return seq_path

    def run_module(self):
        seq_file = self.generate_seq_file()
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'seq_file': seq_file,
        }
        if self.option('reads_type') == 'second':
            opts.update({'method': 'bwa'})
        elif self.option('reads_type') == 'third':
            opts.update({'method': 'minimap2'})
        self.target_depth.set_options(opts)
        self.target_depth.on('end', self.set_db)
        self.target_depth.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        target_depth_api = self.api.api('tool_lab.target_depth')
        data = glob.glob(os.path.join(self.target_depth.output_dir, 'uploads', '*.xls'))
        target_depth_api.add_target_depth(out_table=data, seq_len=self.seq_len, main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(os.path.join(self.target_depth.output_dir, 'uploads'))
        result_dir.add_relpath_rules([
            [".", "", "原始reads目的基因查询结果", 0],
            [r'*.xls', 'xls', '目的基因测序深度详情表', 0],
        ])
        super(TargetDepthWorkflow, self).end()


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
            fastq_dir='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/target_depth/minimap2',
            reads_type='third',
            seq_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/308_TagGeneInRawSeq/script/minimap/NDM_5.fasta',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.target_depth",
            main_table_name="sg_target_depth",
            task_id="target_depth",
            project_sn="target_depth",
            submit_location="target_depth"
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