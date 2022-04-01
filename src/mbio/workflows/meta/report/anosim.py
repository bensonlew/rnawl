# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

"""anosim/adonis以及箱线图和有参检验无参检验"""

import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
from mbio.packages.beta_diversity.filter_newick import get_level_newicktree
import datetime
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class AnosimWorkflow(Workflow):
    """
    报告中使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(AnosimWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "method", "type": "string", "default": 'bray_curtis'},
            {"name": "update_info", "type": "string"},
            {"name": "otu_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "group_file", "type": "infile", 'format': 'meta.otu.group_table'},
            {"name": "group_detail", "type": "string"},
            {"name": "permutations", "type": "int", "default": 999}
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        task = self.add_module("meta.beta_diversity.beta_diversity")
        if 'unifrac' in self.option('method'):  # sanger_bioinfo/src/mbio/workflows/meta/report/distance_calc.py中的解释
            if self.option('level') != 9:
                newicktree = get_level_newicktree(self.option('otu_id'), level=self.option('level'),
                                                  tempdir=self.work_dir, return_file=False, bind_obj=self)
                all_find = re.findall(r'\'.+?\'', newicktree)
                for n, m in enumerate(all_find):
                    all_find[n] = m.strip('\'')
                all_find = dict((i[1], i[0]) for i in enumerate(all_find))

                def match_newname(matchname):
                    if hasattr(match_newname, 'count'):
                        match_newname.count = match_newname.count + 1
                    else:
                        match_newname.count = 1
                    return 'OTU' + str(match_newname.count)
                newline = re.sub(r'\'.+?\'', match_newname, newicktree)
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newline)
                tempfile.close()
                otu_table = self.option('otu_file').path
                temp_otu_file = self.option('otu_file').path + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                if len(all_lines) < 3:
                    self.bind_object.logger.error('分类水平：%s,otu表数据少于2行：%s' % (self.option('level'), len(all_lines)))
                    self.set_error("数据表少于两行", code="12700501")
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:
                    name = line.split('\t')
                    if name[0] not in all_find:
                        self.set_error('OTU表中存在不是直接通过组合原始表分类名称的OTU名：%s', variables=(name[0]), code="12700502")
                    name[0] = 'OTU' + str(all_find[name[0]] + 1)
                    new_all.append('\t'.join(name))
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                options = {
                    'dis_method': self.option('method'),
                    'otutable': temp_otu_file,
                    'phy_newick': temp_tree_file
                }
                options['otutable'] = self.filter_otu_sample(options['otutable'],
                                                             self._get_samplenames(self.option('group_file').path),
                                                             options['otutable'] + '.temp')
            else:
                newicktree = get_level_newicktree(self.option('otu_id'), level=self.option('level'),
                                                  tempdir=self.work_dir, return_file=False, bind_obj=self)
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newicktree)
                tempfile.close()
                otu_table = self.option('otu_file').path
                temp_otu_file = self.option('otu_file').path + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:  # OTU表中有复杂的名称OTU名称，包含进化物种类型，进化树种只有OTU名称
                    name = line.split('\t')
                    name[0] = name[0].split(';')[-1].strip()
                    new_all.append('\t'.join(name))
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                options = {
                    'dis_method': self.option('method'),
                    'otutable': temp_otu_file,
                    'phy_newick': temp_tree_file
                    }
                options['otutable'] = self.filter_otu_sample(options['otutable'],
                                                             self._get_samplenames(self.option('group_file').path),
                                                             options['otutable'] + '.temp')
        else:
            options = {
                'dis_method': self.option('method'),
                'otutable': self.option('otu_file')
                }
            options['otutable'] = self.filter_otu_sample(options['otutable'].path,
                                                         self._get_samplenames(self.option('group_file').path),
                                                         options['otutable'].path + '.temp')
        options['analysis'] = 'anosim'
        options['permutations'] = self.option('permutations')
        options['group'] = self.option('group_file')
        task.set_options(options)
        task.on('end', self.set_db)
        task.run()
        self.output_dir = task.output_dir
        super(AnosimWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        api_anosim = self.api.anosim
        if not (os.path.isdir(self.output_dir + '/Anosim') and os.path.isdir(self.output_dir + '/AnosimBox')): # change by wzy
            self.logger.error("找不到报告文件夹:{}".format(self.output_dir))
            self.set_error("找不到报告文件夹", code="12700503")
        api_anosim.add_beta_anosim_result(self.output_dir, main=False, main_id=self.option('main_id'))
        self.logger.info('运行self.end')
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.logger.info("get_save_pdf_status：{}".format(get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))))
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_beta_multi_anosim")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "beta_anosim",
                "interaction": 1,
                "main_table": "sg_beta_multi_anosim",
            })
            self.figsave.run()
        else:
            self.end()


    def filter_otu_sample(self, otu_path, filter_samples, newfile):
        if not isinstance(filter_samples, types.ListType):
            self.logger.error('过滤otu表样本的样本名称应为列表')
            self.set_error("filter_otu_sample错误", code="12700504")
        try:
            with open(otu_path, 'rb') as f, open(newfile, 'wb') as w:
                one_line = f.readline()
                all_samples = one_line.rstrip().split('\t')[1:]
                if not ((set(all_samples) & set(filter_samples)) == set(filter_samples)):
                    self.logger.error('提供的过滤样本存在otu表中不存在的样本all:%s,filter_samples:%s' % (all_samples, filter_samples))
                    self.set_error("样本名称不存在otu表中", code="12700505")
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
            self.set_error('无法打开OTU相关文件或者文件不存在', code="12700506")

    def _get_samplenames(self, groupfile):
        try:
            with open(groupfile, 'rb') as f:
                alllines = f.readlines()
                all_names = [i.split('\t')[0] for i in alllines]
            return all_names[1:]
        except IOError:
            self.set_error('无法打开分组文件或者文件不存在', code="12700507")

    def end(self):
        save_params(self.output_dir, self.id)
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir + "/AnosimBox"))
        repaths = [
            [".", "", "ANOSIM&Adonis分析结果目录", 0, "110112"],
            ["Anosim", "", "ANOSIM&Adonis结果输出目录", 0, "110113"],
            ["Anosim/anosim_results.txt", "txt", "ANOSIM分析结果", 0, "110115"],
            ["Anosim/adonis_results.txt", "txt", "Adonis分析结果", 0, "110114"],
            ["Anosim/format_results.xls", "xls", "ANOSIM&Adonis整理结果表", 0, "110116"],
            ["AnosimBox", "", "ANOSIM箱式图数据表", 0, "110117"],
            ["AnosimBox/box_data.xls", "xls", "箱式图数据", 0, "110118"],
            ["Anosim/new_adonis_results.xls", "xls", "Adonis分析结果", 0, "110114"],
            #["Box/Distances.xls", "xls", "组内组间距离值统计结果"],
            ["Distance", "", "距离矩阵计算结果输出目录", 0, "110089"],
            ["AnosimBox/ANOSIM分析箱式图.pdf", "pdf", "ANOSIM分析箱式图", 0, ""]
            ]
        regexps = [
            [r'Distance/%s.*\.xls$' % self.option('method'), 'xls', '样本距离矩阵文件', 0, "110090"]
            ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(AnosimWorkflow, self).end()
