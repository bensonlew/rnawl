# -*- coding: utf-8 -*-
# __author__ = 'xuxi_20210729'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError


class BiclusterWorkflow(Workflow):
    """
    bicluster base R package bicluster
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BiclusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "infile", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "preprocess", "type": "string", 'default': 'no'},
            {"name": "format", "type": "string", 'default': "read_counts"},
            {"name": "missing", "type": "string", 'default': "geneMedian"},
            {"name": "ngenes", "type": "int", 'default': 2000},
            {"name": "method", "type": "string", 'default': 'BCCC'},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.bicluster = self.add_tool("tool_lab.bicluster")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_bicluster()
        super(BiclusterWorkflow, self).run()

    def check_options(self):
        if not self.option("infile").is_set:
            raise OptionError("必须设置输入表达量文件")
        return True

    def run_bicluster(self):
        opts = {
            "infile":self.option('infile').prop['path'],
            "preprocess":self.option('preprocess'),
            "format":self.option('format'),
            "missing":self.option('missing'),
            "ngenes":self.option('ngenes'),
            "method":self.option('method')
        }
        self.bicluster.set_options(opts)
        self.bicluster.on('end', self.set_db)
        self.bicluster.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        try:
            result_pdf_num = int(open(os.path.join(self.bicluster.output_dir, 'bicluster.num'),'r').read().strip())
        except:
            result_pdf_num = 0
        bicluster_api = self.api.api('tool_lab.bicluster')
        bicluster_api.add_bicluster(
            pdf_num = result_pdf_num, 
            tool_output_dir = self.bicluster.output_dir,
            s3_dir = self._sheet.output,
            main_id = self.option('main_id'),
        )
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.bicluster.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "双聚类",0],
            ['./bicluster.num', '', '子聚类数', 0],
            ['./sub_cluster_*.tsv', '', '子聚类矩阵', 0],
            ['./sub_cluster_*.pdf', '', '子聚类热图', 0],
        ])
        super(BiclusterWorkflow, self).end()


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
            vcf_path='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/bsa/pop.final.vcf',
            # wp='',
            # mp='HQS1',
            mb='XS11_1',
            wb='F44_mix',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.bsa",
            main_table_name="sg_bsa",
            task_id="bsa",
            project_sn="bsa",
            submit_location="bsa"
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