# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'

"""PERMANOVA分析"""

import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
from mbio.packages.beta_diversity.filter_newick import get_level_newicktree
import datetime
from mbio.packages.meta.common_function import envname_restore
from mbio.packages.meta.save_params import save_params


class PermanovaWorkflow(Workflow):
    """
    报告中使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(PermanovaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "method", "type": "string", "default": 'bray_curtis'},
            {"name": "update_info", "type": "string"},
            {"name": "otu_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "group_detail", "type": "string"},
            {"name": "permutation", "type": "int", "default": 999},
            {"name": "main_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "group_detail_list", "type": "string"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "env_labs", "type": "string", "default": ""},
            {"name": "env_id", "type": "string", "default": ""}
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.permanova = self.add_tool("meta.beta_diversity.permanova")

    def run(self):
        [binary, distance_method] = self.convert_dis_method(self.option('method'))
        options = {
            'dis_method': distance_method,
            'permutations': self.option('permutation'),
            'envtable': self.option('envtable'),
            'binary': binary
        }
        options['abu_table'] = self.filter_otu_sample(self.option('otu_file').prop['path'],
                                                      self._get_samplenames(self.option('group_table').path),
                                                      self.work_dir + '/temp.xls')
        self.permanova.set_options(options)
        self.permanova.on('end', self.set_db)
        self.permanova.run()
        self.output_dir = self.permanova.output_dir
        super(PermanovaWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        api_anosim = self.api.anosim
        if not (os.path.isdir(self.output_dir + '/Permanova') and os.path.isdir(self.output_dir + '/PermanovaBox')): # change by wzy
            self.logger.error("找不到报告文件夹:{}".format(self.output_dir))
            self.set_error("找不到报告文件夹", code="12703401")
        api_anosim.add_beta_anosim_result(self.output_dir, main=False, main_id=self.option('main_id'))
        self.logger.info('运行self.end')
        self.end()

    def convert_dis_method(self, dis_method):
        distance_name = ['euclidean', 'binary_euclidean', 'manhattan', 'binary_manhattan', 'gowerM', 'binary_gowerM', 
                              'altGower', 'binary_altGower', 'canberraNZ', 'binary_canberraNZ', 'bray_curtis', 'binary_bray_curtis', 
                              'kulczynski', 'binary_kulczynski', 'morisita_horn', 'binary_morisita_horn', 'morisita', 
                              'binomial', 'binary_binomial', 'cao', 'binary_cao', 'chao', 'jaccard', 'binary_jaccard',
                              'raup_crick', 'mountford', 'mahalanobis']
        if dis_method not in distance_name:
            self.set_error('不支持所选的距离算法', code="12703402")
        splite = dis_method.split("_")
        if splite[0] == "binary":
            binary = "true"
            distance_method = dis_method[7:]
        else:
            binary = "false"
            distance_method = dis_method
        if distance_method == 'gowerM':
            distance_method = 'gower'
        elif distance_method == 'bray_curtis':
            distance_method = 'bray'
        elif distance_method == 'morisita_horn':
            distance_method = 'horn'
        elif distance_method == 'raup_crick':
            distance_method = 'raup'
        elif distance_method == 'canberraNZ':
            distance_method = 'canberra'
        return binary, distance_method

    def filter_otu_sample(self, otu_path, filter_samples, newfile):
        if not isinstance(filter_samples, types.ListType):
            self.logger.error('过滤otu表样本的样本名称应为列表')
            self.set_error("filter_otu_sample错误", code="12703403")
        try:
            with open(otu_path, 'rb') as f, open(newfile, 'wb') as w:
                one_line = f.readline()
                all_samples = one_line.rstrip().split('\t')[1:]
                if not ((set(all_samples) & set(filter_samples)) == set(filter_samples)):
                    self.logger.error('提供的过滤样本存在otu表中不存在的样本all:%s,filter_samples:%s' % (all_samples, filter_samples))
                    self.set_error("otu表中不存在过滤样本", code="12703404")
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
            self.set_error('无法打开OTU相关文件或者文件不存在', code="12703405")

    def _get_samplenames(self, groupfile):
        try:
            with open(groupfile, 'rb') as f:
                alllines = f.readlines()
                all_names = [i.split('\t')[0] for i in alllines]
            return all_names[1:]
        except IOError:
            self.set_error('无法打开分组文件或者文件不存在', code="12703406")

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        api_permanova = self.api.permanova
        api_permanova.add_permanova(self.output_dir + '/permanova.xls', main=False,
                                    main_id=self.option('main_id'))
        self.logger.info('运行self.end')
        self.end()

    @envname_restore
    def end(self):
        repaths = [
            [".", "", "PERMANOVA分析结果目录", 0, "110127"],
            ["./permanova.xls", "xls", "PERMANOVA分析结果", 0, "110128"]
        ]
        save_params(self.output_dir, self.id)
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(PermanovaWorkflow, self).end()
