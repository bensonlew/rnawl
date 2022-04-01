# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
from mbio.packages.metaasv.filter_newick import get_level_newicktree
from mbio.packages.metaasv.common_function import link_dir
import datetime


class PermanovaWorkflow(Workflow):
    """
    metaasv PERMANOVA分析
    报告中使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(PermanovaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", "format": "meta.otu.otu_table"},#ASV表
            {"name": "method", "type": "string", "default": 'bray_curtis'},#距离方法
            {"name": "update_info", "type": "string"},
            {"name": "asv_id", "type": "string"},#asv主表的id
            {"name": "env_id", "type": "string"},#asv主表的id
            {"name": "main_id", "type": "string"},#插入主表的ID
            {"name": "level", "type": "int"},#分类水平
            {"name": "env_labs", "type": "string"},#分类水平
            {"name": "group_detail", "type": "string", "default": ""},#group_detail
            {"name": "group_detail_list", "type": "string"},
            {"name": "permutation", "type": "int", "default": 999},#置换次数
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},#环境因子表
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},##分组表
            {"name": "group_id", "type": "string"},#插入主表的group_id
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.permanova = self.add_tool("metaasv.permanova")

    def check_option(self):
        """
        参数二次检查
        :return:
        """
        if self.option("level") in range(10):
            self.set_error("分类水平不存在，请检查！")
        if self.option("otu_file").is_set:
            self.set_error("ASV表不存在，请输入正确的丰度文件！")

    def run_permanova(self):
        """
        运行
        :return:
        """
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

    def run(self):
        """
        运行
        :return:
        """
        self.run_permanova()
        super(PermanovaWorkflow, self).run()

    def convert_dis_method(self, dis_method):
        """
        将距离算法进行转换
        :param dis_method:
        :return:
        """
        distance_name = ['euclidean', 'binary_euclidean', 'manhattan', 'binary_manhattan', 'gowerM', 'binary_gowerM',
                        'altGower', 'binary_altGower', 'canberraNZ', 'binary_canberraNZ', 'bray_curtis', 'binary_bray_curtis',
                        'kulczynski', 'binary_kulczynski', 'morisita_horn', 'binary_morisita_horn', 'morisita',
                        'binomial', 'binary_binomial', 'cao', 'binary_cao', 'chao', 'jaccard', 'binary_jaccard',
                        'raup_crick', 'mountford', 'mahalanobis']
        if dis_method not in distance_name:
            self.set_error('不支持所选的距离算法')
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
        """
        根据otu表和group表过滤掉样本
        :param otu_path: otu表
        :param filter_samples: 要过滤的样本list
        :param newfile: 生成的结果表
        :return:
        """
        if not isinstance(filter_samples, types.ListType):
            self.logger.error('过滤otu表样本的样本名称应为列表')
            self.set_error("filter_otu_sample错误")
        try:
            with open(otu_path, 'rb') as f, open(newfile, 'wb') as w:
                one_line = f.readline()
                all_samples = one_line.strip().split('\t')[1:]
                if not ((set(all_samples) & set(filter_samples)) == set(filter_samples)):
                    self.logger.error('提供的过滤样本存在otu表中不存在的样本all:%s,filter_samples:%s' % (all_samples, filter_samples))
                    self.set_error("otu表中不存在过滤样本")
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
            self.set_error('无法打开OTU相关文件或者文件不存在')

    def _get_samplenames(self, groupfile):
        try:
            with open(groupfile, 'rb') as f:
                alllines = f.readlines()
                all_names = [i.split('\t')[0] for i in alllines]
            return all_names[1:]
        except IOError:
            self.set_error('无法打开分组文件或者文件不存在')

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        link_dir(self.permanova.output_dir, self.output_dir)
        api_permanova = self.api.api("metaasv.permanova")
        api_permanova.add_permanova(self.output_dir, main=False,main_id=self.option('main_id'))
        self.logger.info('运行self.end')
        self.end()

    def end(self):
        repaths = [
            ["../", "", "PERMANOVA分析结果文件", 0, ""],
            [".", "", "PERMANOVA分析结果目录", 0, ""],
            ["./One-factor_PERMANOVA.xls", "xls", "单因素PERMANOVA分析结果表", 0, ""],
            ["./Multi-factor_PERMANOVA.xls", "xls", "多因素PERMANOVA分析结果表", 0, ""],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(PermanovaWorkflow, self).end()
