# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class EnrichmentWorkflow(Workflow):
    """
    hypergeometric test function for gene set enrichment analysis that are designed to accept user defined annotation
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EnrichmentWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene", "type": "infile", "format": "ref_rna_v2.common"},  # gene list file
            {"name": "universe", "type": "infile", "format": "ref_rna_v2.common"},
            # background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background
            {"name": "term2gene", "type": "infile", "format": "ref_rna_v2.common"},
            # user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
            {"name": "term2name", "type": "infile", "format": "ref_rna_v2.common"},
            # user input of TERM TO NAME mapping, a data.frame of 2 column with term and name
            {"name": "pvalue_cutoff", "type": "float", "default": 1.0},
            {"name": "qvalue_cutoff", "type": "float", "default": 1.0},
            {"name": "padjust_method", "type": "string", "default": 'BH'},
            # one of holm, hochberg, hommel, bonferroni, BH, BY, fdr, none
            {"name": "min", "type": "int", "default": 10},  # minimal size of genes annotated for testing
            {"name": "max", "type": "int", "default": 500},  # maximal size of genes annotated for testing
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.enrichment")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(EnrichmentWorkflow, self).run()

    def check_options(self):
        if not self.option('gene').is_set:
            raise OptionError('基因集文件必须输入')
        if not self.option('term2gene').is_set:
            raise OptionError('基因注释信息文件必须输入')
        return True

    def run_tool(self):
        if self.option('padjust_method') == "BH/fdr":
            self.option('padjust_method', 'BH')
        opts = {
            'gene': self.option('gene'),
            'term2gene': self.option('term2gene'),
            'pvalue_cutoff': self.option('pvalue_cutoff'),
            'qvalue_cutoff': self.option('qvalue_cutoff'),
            'padjust_method': self.option('padjust_method'),
            'min': self.option('min'),
            'max': self.option('max'),
        }
        if self.option("universe").is_set:
            opts.update({'universe': self.option('universe')})
        if self.option("term2name").is_set:
            opts.update({'term2name': self.option('term2name')})
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            if os.path.exists(os.path.join(self.output_dir, file)):
                os.remove(os.path.join(self.output_dir, file))
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def set_db(self):
        api = self.api.api('tool_lab.enrichment')
        api.add_enrich_detail(self.tool.output_dir + "/enrichment.txt", self.option('main_id'))
        self.set_output()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "富集分析结果文件",0],
            [r'enrichment.txt', 'txt', '富集分析结果文件', 0],
        ])
        super(EnrichmentWorkflow, self).end()


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
            gene="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/query.list",
            term2gene="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/kegg_query.list",
            term2name='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/TERM2NAME',
            universe='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/backgroud',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.enrichment",
            main_table_name="enrichment",
            task_id="enrichment",
            project_sn="enrichment",
            submit_location="enrichment"
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
