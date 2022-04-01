# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
import re
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class PathviewAgent(Agent):
    def __init__(self, parent):
        super(PathviewAgent, self).__init__(parent)
        options = [
            {"name": "gene_file", "type": "infile", 'format': 'ref_rna_v2.common'},
            {"name": "compound_file", "type": "infile", 'format': 'ref_rna_v2.common'},
            {"name": "type", "type": "string", 'default': "gene"},
            {"name": "data_type", "type": "string", 'default': "continuous"},
            {"name": "gene_id_type", "type": "string", 'default': "entrez"},
            {"name": "compound_id_type", "type": "string", 'default': "kegg"},
            {'name': 'species', 'type': 'string', 'default': 'hsa'},
            {'name': 'pathway', 'type': 'string', 'default': None},
        ]
        self.add_option(options)
        self.step.add_steps("pathview")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.pathview.start()
        self.step.update()

    def step_finish(self):
        self.step.pathview.finish()
        self.step.update()

    def check_options(self):
        if self.option('type') == 'gene' and not self.option("gene_file").is_set:
            raise OptionError("必须设置输入基因文件")
        if self.option('type') == 'compound' and not self.option("compound_file").is_set:
            raise OptionError("必须设置输入化合物文件")
        if self.option('type') == 'both' and not self.option("gene_file").is_set:
            raise OptionError("数据类型与输入文件不匹配，请检查。")
        if self.option('type') == 'both' and not self.option("compound_file").is_set:
            raise OptionError("数据类型与输入文件不匹配，请检查。")
        if self.option('gene_id_type').lower() not in ['entrezid', 'entrez'] and self.option('data_type') == 'continuous' and not self.option('pathway'):
            raise OptionError("暂不支持该基因ID类型自动检测KEGG pathway。请输入指定pathwayID，以逗号分割，或转换成EntrezID。")
        if self.option('compound_id_type').lower() != 'kegg' and self.option('data_type') == 'continuous' and not self.option('pathway'):
            raise OptionError("暂不支持该化合物ID类型自动检测KEGG pathway。请输入指定pathwayID，以逗号分割，或转换成KEGG Compound ID。")
        if self.option('data_type') != 'continuous' and not self.option('pathway'):
            raise OptionError("请输入指定pathwayID，以逗号分割。")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        super(PathviewAgent, self).end()


class PathviewTool(Tool):
    def __init__(self, config):
        super(PathviewTool, self).__init__(config)
        self._version = "v1.0"
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/bin/'))
        self.program = {
            'rscript': 'bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/bin/Rscript'
        }
        self.script = {
            'pathview': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/pathview.r')
        }
        self.species = self.option('species')

    def run(self):
        super(PathviewTool, self).run()
        self.file_check()
        self.run_pathview()
        self.set_output()
        self.end()

    def file_check(self):
        if self.option('type') == 'both' and self.option('data_type') == 'continuous':
            g_df = pd.read_table(self.option('gene_file').prop['path'])
            c_df = pd.read_table(self.option('compound_file').prop['path'])
            if g_df.shape[1] != c_df.shape[1]:
                self.logger.info('样本数量在基因文件与化合物文件中不一致，')
        if self.option('pathway'):
            pattern = re.compile(r'[a-z]+')
            species = pattern.findall(self.option('pathway'))
            if len(set(species)) > 1:
                self.set_error('输入pathwayID包含多个物种，请检查')
            elif len(set(species)) == 1 and species[0] not in [self.option('species'), 'ko']:
                self.logger.info('输入pathwayID可能与自定义物种不一致，后续分析以pathwayID物种进行。')
                self.species = species[0]

    def run_pathview(self):
        ko_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Annotation/all/KEGG/version_202007/data', self.option('species'))
        kegg_path = os.path.join(self.config.SOFTWARE_DIR, 'database', 'pathview')
        if not os.path.exists(kegg_path):
            os.makedirs(kegg_path)
        cmd = '{} {} '.format(self.program['rscript'], self.script['pathview'])
        if self.option('gene_file').is_set:
            cmd += '-g {} '.format(self.option('gene_file').prop['path'])
            cmd += '-i {} '.format(self.option('gene_id_type'))
        if self.option('compound_file').is_set:
            cmd += '-c {} '.format(self.option('compound_file').prop['path'])
            cmd += '-d {} '.format(self.option('compound_id_type'))
        if self.option('pathway'):
            pattern = re.compile(r'\d+')
            pathways = pattern.findall(self.option('pathway'))
            cmd += '-p {} '.format(','.join(pathways))
        if not self.option('pathway') or self.option('pathway').startswith('ko'):
            cmd += '-k {} '.format(ko_path)
        if self.option('pathway') and self.option('pathway').startswith('ko'):
            cmd += '-y '
        cmd += '-e {} '.format(kegg_path)
        cmd += '-s {} -t {}'.format(self.species, self.option('data_type'))
        cmd_name = 'run_pathview'
        command = self.add_command(cmd_name, cmd)
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

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        files = glob.glob(os.path.join(self.work_dir, '*pathview.multi.png'))
        files += glob.glob(os.path.join(self.work_dir, '*_GAGE_result.txt'))
        for each in files:
            new_path = os.path.join(self.output_dir, os.path.basename(each))
            os.link(each, new_path)

