# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang '


import os
import shutil
import re
import pandas as pd
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
from mbio.packages.metaasv.common_function import link_dir
from mbio.packages.metaasv.filter_newick import get_level_newicktree


class CompositionHeatmapWorkflow(Workflow):
    """
    metaasv 群落组成分析heatmap图模块
    做的改变是：调用module--heatmap去做，输入otu表带有上级水平的物种名称，便于进行着色的筛选
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CompositionHeatmapWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "input_otu_id", "type": "string"},  # 输入的OTU id
            {"name": "level", "type": "string", "default": "9"},  # 输入的ASV level
            {"name": "group_detail", "type": "string"},  # 输入的group_detail 示例如下
            {"name": "species_number", "type": "string", "default": ""},  # 物种数目，默认全部物种
            {"name": "method", "type": "string", "default": ""},  # 物种层次聚类方式，默认不聚类
            {"name": "sample_method", "type": "string", "default": ""},  # 样本层次聚类方式，默认不聚类
            {"name": "add_Algorithm", "type": "string", "default": ""},  # 分组样本求和算法，默认不求和
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},##group表
            {"name": "level_color", "type": "string", "default": "3"} ,##颜色水平
            {"name": "sample_distance", "type": "string", "default": "bray_curtis"},
            {"name": "species_distance", "type": "string", "default": "bray_curtis"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.heatmap = self.add_module("metaasv.heatmap")
        self.rename = {}
        self.level_dict = {
            "1": "d__",
            "2": "k__",
            "3": "p__",
            "4": "c__",
            "5": "o__",
            "6": "f__",
            "7": "g__",
            "8": "s__",
            "9": "ASV",
        }

    def check_options(self):
        if self.option('method') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的物种层次聚类方式：%s', variables=(self.option('method')))
        if self.option('sample_method') not in ['average', 'single', 'complete', ""]:
            raise OptionError('错误的样本层次聚类方式：%s', variables=(self.option('sample_method')))
        if self.option('add_Algorithm') not in ['sum', 'average', 'middle', ""]:
            raise OptionError('错误的样本求和方式：%s', variables=(self.option('add_Algorithm')))
        if self.option("method") != "" :
            if (self.option("species_number") != "" and self.option("species_number") != "0"):
                if int(self.option("species_number")) == 1:
                    raise OptionError('物种聚类的个数不能为：%s', variables=(self.option('species_number')))

    def run_heatmap(self):
        """
        运行 heatmap module
        """
        input_table = self.change_otuname(self.option("in_otu_table").prop['path'], "asv")
        # input_table = self.option("in_otu_table").prop['path']
        opts = ({
            "abundtable": input_table,
            "group": self.option("group"),
            "species_number": self.option("species_number"),
            "method": self.option("method"),
            "analysis": "heatmap",
            "fill_zero": "true",
            "sample_method": self.option("sample_method"),
            "add_Algorithm": self.option("add_Algorithm"),
            "sample_distance": self.option('sample_distance'),
            "species_distance": self.option('species_distance'),
        })
        if 'unifrac'in self.option("species_distance"):
            self.logger.info("开始从原始asv_id中挑选出进化树的距离文件！")
            input_table, tree_file= self.get_newicktree(input_table)
            opts["phy_newick"] = tree_file
            opts["abundtable"] = input_table
            opts["species_number"] = "" ## unifrac计算方法不能选择top物种，因为这个地方需要根据筛选的进化树距离做，不能根据top物种做
        self.heatmap.set_options(opts)
        self.heatmap.on("end", self.set_db)
        self.heatmap.run()

    def get_newicktree(self, input_table):
        """
        从进化树中抽取距离的方法，然后后面进行计算
        :return:
        """
        if 'unifrac' in self.option('species_distance'):  # sanger_bioinfo/src/mbio/workflows/meta/report/distance_calc.py中的解释
            if self.option('level') != 9:
                newicktree = get_level_newicktree(self.option('input_otu_id'), level=int(self.option('level')),tempdir=self.work_dir, return_file=False,topN=self.option("species_number"), bind_obj=self)
                all_find = re.findall(r'\'.+?\'', newicktree)
                for n, m in enumerate(all_find):
                    all_find[n] = m.strip('\'')
                all_find = dict((i[1], i[0]) for i in enumerate(all_find))

                def match_newname(matchname):
                    if hasattr(match_newname, 'count'):
                        match_newname.count = match_newname.count + 1
                    else:
                        match_newname.count = 1
                    return 'ASV' + str(match_newname.count)
                newline = re.sub(r'\'.+?\'', match_newname, newicktree)
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newline)
                tempfile.close()
                self.logger.info('get_newick:' + temp_tree_file)
                otu_table = input_table
                temp_otu_file = input_table + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                if len(all_lines) < 3:
                    self.logger.error('分类水平：%s,ASV表数据少于2行：%s' % (self.option('level'), len(all_lines)))
                    self.set_error("ASV表数据少于2行")
                self.logger.info(len(all_lines))
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:
                    name = line.split('\t')
                    origin_name = name[0].split("; ")[-1].strip()
                    if name[0] in all_find:
                        name[0] = 'ASV' + str(all_find[name[0]] + 1)
                    new_all.append('\t'.join(name))
                    if name[0] not in self.rename:
                        self.rename[name[0]] = origin_name
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                new_input_table = os.path.join(self.work_dir, "otu_file.xls")
                os.rename(otu_table, os.path.join(self.work_dir, "otu_table2.xls"))
                os.rename(temp_otu_file, new_input_table)
            else:
                newicktree = get_level_newicktree(self.option('input_otu_id'), level=int(self.option('level')),
                                                  tempdir=self.work_dir, return_file=False, bind_obj=self)
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newicktree)
                tempfile.close()
                otu_table = input_table
                temp_otu_file = input_table + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:  # OTU表中有复杂的名称OTU名称，包含进化物种类型，进化树种只有OTU名称
                    name = line.split('\t')
                    origin_name = name[0].split("; ")[-1].strip()
                    name[0] = name[0].split(';')[-1].strip()
                    new_all.append('\t'.join(name))
                    if name[0] not in self.rename:
                        self.rename[name[0]] = origin_name
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                new_input_table = os.path.join(self.work_dir, "otu_file.xls")
                os.rename(otu_table, os.path.join(self.work_dir, "otu_table2.xls"))
                os.rename(temp_otu_file, new_input_table)
            return (new_input_table, temp_tree_file)

    def change_otuname(self, tablepath, file_name):
        """
        改换asv名称
        :param tablepath:
        :param file_name:
        :return:
        """
        self.rename_dict = {}
        newtable = self.work_dir + "/" + file_name + "_input_abund.xls"
        with open(tablepath, "r") as f, open(newtable+"_tmp", "w") as g:
            head = f.readline()
            g.write(head)
            for line in f:
                lines = line.split("\t", 1)
                origin_name = lines[0]
                specimen = re.subn("^.*; ", "", lines[0])[0]
                if specimen not in self.rename_dict:
                    self.rename_dict[specimen] = origin_name
                g.write(specimen + "\t" + lines[1])
        with open(newtable+"_tmp", "r") as y, open(newtable, "w") as t:
            all_name = []
            data_dict = {}
            all_data = y.readlines()
            for i in all_data:
                all_name.append(i.strip().split("\t")[0])
            raw_num = len(all_name)
            if len(set(all_name)) == raw_num:
                for v in all_data:
                    t.write(v)
            else:
                t.write(all_data[0])
                for m in all_data[1:]:
                    if m.strip().split("\t")[0] not in data_dict.keys():
                        data_dict[m.strip().split("\t")[0]] = m.strip().split("\t")[1:]
                    else:
                        for x in range(len(data_dict[m.strip().split("\t")[0]])):
                            data_dict[m.strip().split("\t")[0]][x] = str(int(data_dict[m.strip().split("\t")[0]][x]) + int(m.strip().split("\t")[1:][x]))
                for n in data_dict.keys():
                    t.write(n + "\t" + "\t".join(data_dict[n]) + "\n")
        return newtable

    def replace_name(self, input):
        """
        对文件替换为原来的名称
        :return:
        """
        out_table = os.path.join(self.work_dir, "out_table.xls")
        with open(input, "r") as f, open(out_table, "w") as w:
            lines = f.readlines()
            w.write(lines[0])
            for line in lines[1:]:
                line = line.strip().split("\t")
                sp_name = line[0].split("; ")[-1].strip()
                if sp_name in self.rename:
                    line[0] = self.rename[sp_name]
                w.write("\t".join(line) + "\n")
        os.remove(input)
        os.rename(out_table, input)

    def convert_format(self, input_file, outout_file):
        """
        metaasv 筛选leve_color
        :return:
        """
        level_color = self.option('level_color')
        self.logger.info("+++++++++{}".format(self.rename_dict))
        level_color_number = int(level_color) - 1
        with open(input_file, 'r') as f, open(outout_file, "w") as w:
            lines = f.readlines()
            w.write(lines[0].strip() + "\tlevel_color\n")
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[0] in self.rename_dict.keys():
                    asv_name = self.rename_dict[line[0]]
                    asv_name_list = asv_name.strip().split("; ")
                    level_name = asv_name_list[level_color_number]
                    w.write("\t".join(line) + "\t{}\n".format(level_name))

    def set_db(self):
        """
        数据导入MongoDB和链接结果文件
        :return:
        """
        link_dir(self.heatmap.output_dir, self.output_dir)
        os.remove(os.path.join(self.output_dir,"taxa.percents.mongo.xls"))
        insert_asv_table = os.path.join(self.output_dir, "taxa.table.xls")
        output_asv_table = os.path.join(self.work_dir, "taxa.table.mongo.xls")
        insert_asv_percents_table = os.path.join(self.output_dir, "taxa.percents.table.xls")
        output_asv_percents_table = os.path.join(self.work_dir, "taxa.percents.mongo.xls")
        self.convert_format(insert_asv_table, output_asv_table)
        self.convert_format(insert_asv_percents_table, output_asv_percents_table)
        if 'unifrac' in self.option('species_distance'):
            self.replace_name(output_asv_table)
            self.replace_name(output_asv_percents_table)
        sample_list = []
        species_list = []
        self.logger.info("正在写入mongo数据库")
        species_tree_path = os.path.join(self.output_dir, "species_hcluster.tre")
        if self.option("method") != "":
            species_tree_path2 = os.path.join(self.output_dir, "species_hcluster2.tre")
            if os.path.exists(species_tree_path):
                if 'unifrac' in self.option('species_distance'):
                    with open(species_tree_path, "r") as f, open(species_tree_path2, "w") as w:
                        species_tree = f.readline().strip()
                        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', species_tree)
                        species_list = [i[1].split("; ")[-1].strip() for i in raw_samp]
                        self.logger.info("species_tree: {} \n {}".format(raw_samp, species_list))
                        for i in raw_samp:
                            spe_index = raw_samp.index(i)
                            new_name = str(species_list[spe_index])
                            if new_name in self.rename:
                                new_name = self.rename[new_name]
                            species_tree = species_tree.replace(str(i), new_name)
                        w.write(species_tree+"\n")
                else:
                    with open(species_tree_path, "r") as f, open(species_tree_path2, "w") as w:
                        species_tree = f.readline().strip()
                        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', species_tree)
                        species_list = [i[1].split("; ")[-1].strip() for i in raw_samp]
                        self.logger.info("species_tree: {} \n {}".format(raw_samp, species_list))
                        for i in raw_samp:
                            spe_index = raw_samp.index(i)
                            new_name = str(species_list[spe_index])
                            species_tree = species_tree.replace(str(i), new_name)
                        w.write(species_tree+"\n")
            os.remove(species_tree_path)
            os.rename(species_tree_path2, species_tree_path)
        sample_tree_path = os.path.join(self.output_dir, "specimen_hcluster.tre")
        if self.option("sample_method") != "":
            if os.path.exists(sample_tree_path):
                with open(sample_tree_path, "r") as f:
                    sample_tree = f.readline().strip()
                    raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', sample_tree)
                    sample_list = [i[1].split("; ")[-1].strip() for i in raw_samp]

        if (self.option("add_Algorithm") != "" and self.option("sample_method") == ""):
            sample_name = open(self.option("group").prop["path"], "r")
            content = sample_name.readlines()
            for f in content:
                f = f.strip("\n")
                arr = f.strip().split("\t")
                if arr[0] != "#sample":
                    if arr[0] not in sample_list:
                        sample_list.append(arr[0].split("; ")[-1].strip())

        api_heatmap = self.api.api("metaasv.composition_heatmap")
        api_heatmap.add_heatmap_detail(output_asv_table, self.option("main_id"), "absolute",specimen_sorts=sample_list,species_sorts=species_list)
        api_heatmap.add_heatmap_detail(output_asv_percents_table, self.option("main_id"), "relative",specimen_sorts=sample_list,species_sorts=species_list)
        if os.path.exists(species_tree_path):
            api_heatmap.insert_tree_table(species_tree_path,self.option("main_id"), "species")
        if os.path.exists(sample_tree_path):
            api_heatmap.insert_tree_table(sample_tree_path,self.option("main_id"), "specimen")
        self.end()

    def end(self):
        """
        结束和上传文件
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "群落Heatmap分析结果输出目录", 0, ""],
            ["./heatmap.taxa.table.xls", "xls", "群落Heatmap分析可视化结果数据表", 0, ""],
            ["./sample_hcluster.tre", "tre", "样本聚类树", 0, ""],
            ["./species_hcluster.tre", "tre", "物种聚类树", 0, ""]
        ])
        super(CompositionHeatmapWorkflow, self).end()

    def run(self):
        """
        运行
        :return:
        """
        self.run_heatmap()
        super(CompositionHeatmapWorkflow, self).run()
