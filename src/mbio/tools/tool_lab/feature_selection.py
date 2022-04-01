# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import unittest
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class FeatureSelectionAgent(Agent):
    def __init__(self, parent):
        super(FeatureSelectionAgent, self).__init__(parent)
        options = [
            {"name": "clinical_file", "type": "infile", 'format': 'medical_transcriptome.common'},
            {"name": "exp_file", "type": "infile", "format": "medical_transcriptome.common"},
            {"name": "geneset_file", "type": "infile", 'format': 'medical_transcriptome.common'},
            {"name": "geneset_str", "type": "string", "default": ""},
            {"name": "method", "type": "string", "default": "lasso"},
            {'name': 'genelist_file', 'type': 'outfile', 'format': 'medical_transcriptome.common'},
        ]
        self.add_option(options)
        self.step.add_steps("feature_selection")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.feature_selection.start()
        self.step.update()

    def step_finish(self):
        self.step.feature_selection.finish()
        self.step.update()

    def check_options(self):
        if not self.option("clinical_file").is_set:
            raise OptionError("必须设置输入临床信息文件。")
        if not self.option("exp_file").is_set:
            raise OptionError("必须设置输入表达谱")
        if not self.option("geneset_str") and not self.option('geneset_file').is_set:
            raise OptionError("必须设置输入基因集")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"]
        # ])
        # result_dir.add_regexp_rules([
        #     [r"disgenet_enrichment.xls$", "xls", "DisGeNET富集分析结果"]
        # ])
        super(FeatureSelectionAgent, self).end()


class FeatureSelectionTool(Tool):
    def __init__(self, config):
        super(FeatureSelectionTool, self).__init__(config)
        self._version = "v1.0"
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.feature_selection_path = self.config.PACKAGE_DIR + "/tool_lab/cox_feature_selection.r"
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin/Rscript"

    def run(self):
        super(FeatureSelectionTool, self).run()
        self.get_geneset()
        self.run_feature_selection()
        self.check_result()
        self.set_output()
        self.end()

    def get_geneset(self):
        if self.option('geneset_file').is_set:
            geneset_l = list()
            with open(self.option('geneset_file').prop['path'], 'r') as f:
                for line in f:
                    gene = line.strip()
                    geneset_l.append(gene)
            self.geneset = ','.join(geneset_l)
        elif self.option('geneset_str'):
            self.geneset = self.option('geneset_str')
            geneset_l = self.option('geneset_str').split(',')
            if len(geneset_l) < 2:
                self.set_error("输入基因集过少,或正确未使用逗号分隔基因。请检查。")
        if len(geneset_l) > 1000:
            self.set_error("输入基因集超过最大基因个数限制：1000。请检查。")
        elif len(geneset_l) < 2:
            self.set_error("输入基因集过少。请检查。")

    def run_feature_selection(self):
        clinical = self.option("clinical_file").prop["path"]
        exp = self.option("exp_file").prop["path"]
        cmd = '{} {}'.format(self.r_path, self.feature_selection_path)
        cmd += ' -c {}'.format(clinical)
        cmd += ' -e {}'.format(exp)
        if self.geneset:
            cmd += ' -g {}'.format(self.geneset)
        cmd += ' -m {}'.format(self.option('method'))
        cmd_name = 'cox_feature_selection'
        command = self.add_command(cmd_name, cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))

    def check_result(self):
        genelist = glob.glob(self.work_dir + '/*_filtered_genelist.xls')
        with open(genelist[0], 'r') as f:
            count = 0
            for line in f:
                count += 1
                if count > 1:
                    break
        if count < 2:
            self.set_error("未筛选到生存特征基因，请重新上传数据进行分析！")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        files = glob.glob(self.work_dir + '/*.pdf')
        files += glob.glob(self.work_dir + '/*_filtered_genelist.xls')
        files += glob.glob(self.work_dir + '/genelist_exp.xls')
        for each in files:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
        genelist_file = os.path.join(self.output_dir, '{}_filtered_genelist.xls'.format(self.option('method')))
        self.option("genelist_file", genelist_file)
