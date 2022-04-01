# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modifies 2017.11.16

"""anosim/adonis以及箱线图和有参检验无参检验"""

import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
from biocluster.core.exceptions import OptionError
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name
import json

class AnosimWorkflow(CommTableWorkflow):
    """
    报告中使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(AnosimWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "distance_method", "type": "string", "default": 'bray_curtis'},
            {"name": "permutation", "type": "int", "default": 999},
            {"name": "anno_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "profile_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "env_labs", "type": "string", "default": ""},
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "env_id", "type": "string", "default": ""},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        #self.abundance = self.add_tool("meta.create_abund_table")
        self.anosim = self.add_module("meta.beta_diversity.beta_diversity")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        if not self.option("profile_table").is_set:
            #self.abundance.on('end', self.run_anosim)
            self.run_abundance(self.run_anosim)
        if self.option("profile_table").is_set:
            self.run_anosim()
        else:
            #self.run_abundance()
            self.abundance.run()
        self.output_dir = self.anosim.output_dir
        super(AnosimWorkflow, self).run()

    def run_anosim(self):
        if self.option("profile_table").is_set:
            otutable = self.option("profile_table")
        else:
            otutable = self.abundance.option('out_table')
        options = {
            'dis_method': self.option('distance_method'),
            'otutable': otutable
        }
        options['otutable'] = self.filter_otu_sample(options['otutable'].path,
                                                     self._get_samplenames(self.option('group_table').path),
                                                     self.work_dir + '/temp.xls')
        options['analysis'] = 'anosim'
        options['permutations'] = self.option('permutation')
        options['group'] = self.option('group_table')
        self.anosim.set_options(options)
        self.anosim.on('end', self.set_db)
        self.anosim.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        if os.path.exists(self.output_dir + '/Anosim/adonis_results.txt'):
            os.remove(self.output_dir + '/Anosim/adonis_results.txt')
        if os.path.exists(self.output_dir + '/Anosim/new_adonis_results.xls'):
            os.remove(self.output_dir + '/Anosim/new_adonis_results.xls')
        api_anosim = self.api.api('metagenomic.anosim')
        if not (os.path.isdir(self.output_dir + '/Anosim') and os.path.isdir(
                    self.output_dir + '/AnosimBox')):  # change by wzy
            self.set_error("找不到报告文件夹:%s", variables=(self.output_dir), code="12801001")
        api_anosim.add_anosim(self.output_dir + '/Anosim', self.output_dir + '/AnosimBox', main=False,
                              main_id=self.option('main_id'))
        self.logger.info('运行self.end')
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "anosim")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = json.loads(self.option('params'))["submit_location"]
            # self.logger.info(submit_loc)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def filter_otu_sample(self, otu_path, filter_samples, newfile):
        if not isinstance(filter_samples, types.ListType):
            self.set_error('过滤otu表样本的样本名称应为列表', code="12801002")
        try:
            with open(otu_path, 'rb') as f, open(newfile, 'wb') as w:
                one_line = f.readline()
                all_samples = one_line.rstrip().split('\t')[1:]
                if not ((set(all_samples) & set(filter_samples)) == set(filter_samples)):
                    raise OptionError('提供的过滤样本存在otu表中不存在的样本all:%s,filter_samples:%s', variables=(all_samples, filter_samples), code="12801001")
                if len(all_samples) == len(filter_samples):
                    return otu_path
                samples_index = [all_samples.index(i) + 1 for i in filter_samples]
                w.write('#OTU\t' + '\t'.join(filter_samples) + '\n')
                for line in f:
                    all_values = line.rstrip().split('\t')
                    new_values = [all_values[0]] + [all_values[i] for i in samples_index]
                    w.write('\t'.join(new_values) + '\n')
                return newfile
        except IOError:
            raise OptionError('无法打开OTU相关文件或者文件不存在', code="12801002")

    def _get_samplenames(self, groupfile):
        try:
            with open(groupfile, 'rb') as f:
                alllines = f.readlines()
                all_names = [i.split('\t')[0] for i in alllines]
            return all_names[1:]
        except IOError:
            raise OptionError('无法打开分组文件或者文件不存在', code="12801003")

    def end(self):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                temp_dir = self.output_dir + "/AnosimBox"
                os.system("cp -r {}/* {}/".format(pdf_outs, temp_dir))
        repaths = [
            [".", "", "ANOSIM分析结果目录", 0, "120149"],
            ["Anosim", "", "ANOSIM结果输出目录", 0, "120150"],
            ["Anosim/anosim_results.txt", "txt", "ANOSIM分析结果", 0, "120151"],
            # ["Anosim/adonis_results.txt", "txt", "Adonis分析结果"],
            ["Anosim/format_results.xls", "xls", "ANOSIM结果表，包含统计R值和p值", 0, "120152"],
            ["AnosimBox", "", "ANOSIM箱式图数据表", 0, "120153"],
            ["AnosimBox/box_data.xls", "xls", "箱式图数据", 0, "120154"],
            ["AnosimBox/anosim_box.pdf", "pdf", "组间距离箱线图"],
            # ["Anosim/new_adonis_results.xls", "xls", "Adonis分析结果"],
            # ["Box/Distances.xls", "xls", "组内组间距离值统计结果"],
            ["Distance", "", "距离矩阵计算结果输出目录", 0, "120155"]
        ]
        regexps = [
            [r'Distance/%s.*\.xls$' % self.option('distance_method'), 'xls', '样本距离矩阵文件', 0, "120156"]
        ]

        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(AnosimWorkflow, self).end()
