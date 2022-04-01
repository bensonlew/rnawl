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


class FeatureSelectionWorkflow(Workflow):
    """
    hypergeometric test function for gene set enrichment analysis that are designed to accept user defined annotation
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FeatureSelectionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "clinical_file", "type": "infile", 'format': 'medical_transcriptome.common'},  # clinical data must have columns of sample, survival_status, survival_time
            {"name": "exp_file", "type": "infile", 'format': 'medical_transcriptome.common'},   # file of expression profile
            {'name': 'type', 'type': 'string'},
            {"name": "geneset_str", "type": "string", 'default': ''},  # string of gene names, separate by comma
            {"name": "geneset_file", "type": "infile", 'format': 'medical_transcriptome.common'},   # file of genelist
            {"name": "method", "type": "string", 'default': 'lasso'},   # one of Lasso, Stepwise, Backward, Forward
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.tool = self.add_tool("tool_lab.feature_selection")
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.run_tool()
        super(FeatureSelectionWorkflow, self).run()

    def check_options(self):
        if not self.option("clinical_file").is_set:
            raise OptionError("必须设置输入临床信息文件。")
        else:
            clin_form = os.path.basename(self.option("clinical_file").prop['path']).split('.')[-1]
            if clin_form not in ['txt', 'TXT']:
                raise OptionError("必须设置输入txt格式的临床信息文件。")
        if not self.option("exp_file").is_set:
            raise OptionError("必须设置输入表达谱")
        else:
            exp_form = os.path.basename(self.option("exp_file").prop['path']).split('.')[-1]
            if exp_form not in ['txt', 'TXT']:
                raise OptionError("必须设置输入txt格式的表达谱。")
        if not self.option("geneset_str") and not self.option('geneset_file').is_set:
            raise OptionError("必须设置输入基因集")
        if self.option('geneset_file').is_set:
            geneset_form = os.path.basename(self.option("geneset_file").prop['path']).split('.')[-1]
            if geneset_form not in ['txt', 'TXT']:
                raise OptionError("必须设置输入txt格式的基因集文件。")
        return True

    def run_tool(self):
        opts = {
            'clinical_file': self.option('clinical_file'),
            'exp_file': self.option('exp_file'),
            'method': self.option('method').lower()
        }
        if self.option('geneset_file').is_set:
            opts.update({'geneset_file': self.option('geneset_file')})
        elif self.option('geneset_str'):
            opts.update({'geneset_str': self.option('geneset_str')})
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        feature_selection = self.api.api("tool_lab.feature_selection")
        # add result info
        if self.option('method').lower() == 'lasso':
            r1 = os.path.join(self._sheet.output, 'Lasso_result1.pdf')
            r2 = os.path.join(self._sheet.output, 'Lasso_result2.pdf')
            feature_selection.add_selection(self.option('main_id'), result1=r1, result2=r2)
        else:
            feature_selection.add_selection(self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "生存特征基因筛选结果",0],
            [r'*.xls', 'XLS', '生存特征基因筛选结果文件', 0],
            [r'*.pdf', 'PDF', '生存特征基因筛选结果图', 0],
        ])
        super(FeatureSelectionWorkflow, self).end()


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
            clinical_file="/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/lung_test.txt",
            exp_file="/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/count_test.txt",
            geneset_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/geneset_test.txt',
            method='Lasso'
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.feature_selection",
            main_table_name="feature_selection",
            task_id="feature_selection",
            project_sn="feature_selection",
            submit_location="feature_selection"
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
