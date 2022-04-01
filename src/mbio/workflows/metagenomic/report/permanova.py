# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modifies 2017.11.16

import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
from comm_table import CommTableWorkflow
from mbio.packages.meta.common_function import envname_restore


class PermanovaWorkflow(CommTableWorkflow):
    """
    报告中使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(PermanovaWorkflow, self).__init__(wsheet_object)
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
            {"name": "group_detail_list", "type": "string"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "env_labs", "type": "string", "default": ""},
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "env_id", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        #self.abundance = self.add_tool("meta.create_abund_table")
        self.permanova = self.add_tool("meta.beta_diversity.permanova")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        if not self.option("profile_table").is_set:
            #self.abundance.on('end', self.run_permanova)
            self.run_abundance(self.run_permanova)
        if self.option("profile_table").is_set:
            self.run_permanova()
        else:
            #self.run_abundance()
            self.abundance.run()
        self.output_dir = self.permanova.output_dir
        super(PermanovaWorkflow, self).run()

    def run_permanova(self):
        if self.option("profile_table").is_set:
            otutable = self.option("profile_table")
        else:
            otutable = self.abundance.option('out_table')
        [binary, distance_method] = self.convert_dis_method(self.option('distance_method'))
        options = {
            'dis_method': distance_method,
            'abu_table': otutable,
            'permutations': self.option('permutation'),
            'envtable': self.option('envtable'),
            'binary': binary,
        }
        options['abu_table'] = self.filter_otu_sample(otutable.path,
                                                      self._get_samplenames(self.option('group_table').path),
                                                      self.work_dir + '/temp.xls')
        self.permanova.set_options(options)
        self.permanova.on('end', self.set_db)
        self.permanova.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        api_permanova = self.api.api('metagenomic.permanova')
        api_permanova.add_permanova(self.output_dir + '/permanova.xls', main=False,
                                    main_id=self.option('main_id'))
        self.logger.info('运行self.end')
        self.end()

    @envname_restore
    def end(self):
        repaths = [
            [".", "", "PERMANOVA分析结果目录", 0, "120166"],
            ["./permanova.xls", "xls", "PERMANOVA分析结果", 0, "120167"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(PermanovaWorkflow, self).end()

    def convert_dis_method(self, dis_method):
        distance_name = ['euclidean', 'binary_euclidean', 'manhattan', 'binary_manhattan', 'gowerM', 'binary_gowerM',
                              'altGower', 'binary_altGower', 'canberraNZ', 'binary_canberraNZ', 'bray_curtis', 'binary_bray_curtis',
                              'kulczynski', 'binary_kulczynski', 'morisita_horn', 'binary_morisita_horn', 'morisita',
                              'binomial', 'binary_binomial', 'cao', 'binary_cao', 'chao', 'jaccard', 'binary_jaccard',
                              'raup_crick', 'mountford', 'mahalanobis']
        if dis_method not in distance_name:
            raise OptionError('不支持所选的距离算法', code="12802301")
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
            self.set_error('过滤otu表样本的样本名称应为列表', code="12802301")
        try:
            with open(otu_path, 'rb') as f, open(newfile, 'wb') as w:
                one_line = f.readline()
                all_samples = one_line.rstrip().split('\t')[1:]
                if not ((set(all_samples) & set(filter_samples)) == set(filter_samples)):
                    raise OptionError('提供的过滤样本存在otu表中不存在的样本all:%s,filter_samples:%s',variables=(all_samples, filter_samples), code="12802302")
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
            raise OptionError('无法打开OTU相关文件或者文件不存在', code="12802303")

    def _get_samplenames(self, groupfile):
        try:
            with open(groupfile, 'rb') as f:
                alllines = f.readlines()
                all_names = [i.split('\t')[0] for i in alllines]
            return all_names[1:]
        except IOError:
            raise OptionError('无法打开分组文件或者文件不存在', code="12802304")
