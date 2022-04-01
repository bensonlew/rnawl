# -*- coding: utf-8 -*-
# __author__ = 'Binbin Zhao'
# modified 20190304

import os
import json
import re
import gevent
import time
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class VariantCompareWorkflow(Workflow):
    """
    交互分析：样本比较分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VariantCompareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "sample", "type": "string"},
            {"name": "genotype", "type": "string"},
            {"name": "group", "type": "string"},  # group用于存放所有的group信息。
            # 这里所有需要default值的都需要下面在判断的时候一个个列出，而且需要考虑list当中有的地方是空值，有的地方不空的情况。
            {"name": "marktype", "type": "string"}, # 是diff还是same
            {"name": "vcf_file", "type": "string"},
            {"name": "dep", "type": "string"},
            {"name": "maf", "type": "string"},
            {"name": "ad", "type": "string"},     # 为什么ad在下面的时候变成int格式了
            {"name": "max_miss", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_sn", "type": "string"},
            {"name": "analysis_model", "type": "string"},  # 用于判断分析的是但条件还是多条件, 为single或者multiple
            {"name": "alle_number", "type": "string"},  # 用于显示等位基因数目
            {"name": "eff_type", "type": "string"},  # 用于功效过滤
            # {"name": "fun_type", "type": "string"},  # 用于功能过滤
            {"name": "region", "type": "string"},  # 用于基因组位置过滤
            {"name": "variant_type", "type": "string"},  # 用于显示突变类型
            {"name": "region_type", "type": "string"},  # 用于显示突变类型

          ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.group_path = ""

    def check_options(self):
        if not self.option("vcf_file"):
            raise OptionError("缺少vcf_file参数", code="15500109")

    def get_config(self):
        """
        将tool所有使用的参数包成一个config文件传给Tool。
        这里还需要补充染色体的信息。
        :return:
        """
        os.mkdir(self.work_dir + "/config_file")
        if os.path.exists(os.path.join(self.work_dir + "/config_file", "diff.config")):
            os.remove(os.path.join(self.work_dir + "/config_file", "config.txt"))
        if self.option("analysis_model") == "single":
            sample = json.loads(self.option("sample"))  # 接口中
            # sample = self.option("sample").strip().split(",")
            for s in sample:
                try:
                    sample1 = s.strip().split("|")[0]
                    sample2 = s.strip().split("|")[1]
                except:
                    self.set_error("样本传入格式不正确！")
                with open(os.path.join(self.work_dir + "/config_file", sample1 + "_vs_" + sample2 + "_compare"), "w") \
                        as w:
                    w.write("Variant Type=" + self.option("variant_type") + "\n")
                    w.write("Variant Eff=" + self.option("eff_type") + "\n" )
                    w.write("allele_num=" + self.option("alle_number") + "\n")
                    if self.option("region_type") == "allregion":    # 这里还需要补充为loci的情况。
                        w.write("Region=" + self.option("region") + "\n")
                    elif self.option("region_type") == "custom":
                        region_list = json.loads(self.option("region"))
                        self.logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
                        self.logger.info(type(self.option("region")))
                        self.logger.info(self.option("region"))
                        for region in region_list:
                            w.write("Region=" + region + "\n")
                    elif self.option("region_type") == "location":
                        region_list1 = json.loads(self.option("region"))
                        self.logger.info("-----------------------------------------------")
                        # self.logger.info(region_list1)
                        self.logger.info(type(region_list1))
                        for region1 in region_list1:
                            w.write("Region=" + region1 + "\n")
                    min_dep = self.option("dep").strip().split(",")[0]
                    max_dep = self.option("dep").strip().split(",")[1]
                    if max_dep == "":
                        max_dep = 100000
                    w.write("Sample Diff=" + sample1 + "," + str(min_dep) + "," + str(max_dep) + "," +
                            self.option("genotype") + "," +
                            sample2 + "," + str(min_dep) + "," + str(max_dep) + "," + self.option("genotype") + "," +
                            self.option("marktype"))

        else:
            with open(os.path.join(self.work_dir, "diff.config"), "w") as w:
                w.write("Variant Type=" + self.option("variant_type") + "\n")
                w.write("Variant Eff=" + self.option("eff_type") + "\n")
                w.write("allele_num=" + self.option("alle_number") + "\n")
                # w.write("Region=" + self.option("region") + "\n")
                if self.option("region_type") == "allregion":  # 这里还需要补充为loci的情况。
                    w.write("Region=" + self.option("region") + "\n")
                elif self.option("region_type") == "custom":
                    region_list = json.loads(self.option("region"))
                    self.logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
                    self.logger.info(type(self.option("region")))
                    self.logger.info(self.option("region"))
                    for region in region_list:
                        w.write("Region=" + region + "\n")
                elif self.option("region_type") == "location":
                    region_list1 = json.loads(self.option("region"))
                    self.logger.info("-----------------------------------------------")
                    self.logger.info(region_list1)
                    self.logger.info(type(region_list1))
                    for region1 in region_list1:
                        w.write("Region=" + region1 + "\n")
                sample = json.loads(self.option("sample"))   # 接口中是这个
                # sample = self.option("sample").strip().split(",")
                if not sample == [""]:
                    dep = json.loads(self.option("dep"))
                    # dep = self.option("dep").strip().split(",")
                    genotype = json.loads(self.option("genotype"))
                    # genotype = self.option("genotype").strip().split(",")
                    is_same = json.loads(self.option("marktype"))
                    # is_same = self.option("marktype").strip().split(",")
                    for (x, y, z, m)in zip(sample, dep, genotype, is_same):
                        try:
                            sample1 = x.strip().split("|")[0]
                            sample2 = x.strip().split("|")[1]
                        except BaseException:
                            self.set_error("样本传入格式不正确！")
                        try:
                            dep1_min = y.strip().split(
                                "|")[0].strip().split("-")[0]
                            dep1_max = y.strip().split(
                                "|")[0].strip().split("-")[1]
                            dep2_min = y.strip().split(
                                "|")[1].strip().split("-")[0]
                            dep2_max = y.strip().split(
                                "|")[1].strip().split("-")[1]
                            self.logger.info(dep1_min)
                            self.logger.info(dep1_max)
                            if float(dep1_min) < 0:
                                self.set_error("深度信息必须为正整数！")
                            if dep1_max == "":
                                dep1_max = 100000
                            elif float(dep1_max) < 0:
                                self.set_error("深度信息必须为正整数！")
                            if float(dep2_min) < 0:
                                self.set_error("深度信息必须为正整数！")
                            if dep2_max == "":
                                dep2_max = 100000
                            elif float(dep2_max) < 0:
                                self.set_error("深度信息必须为正整数！")
                        except BaseException:
                            self.set_error("测序深度信息传入格式不正确！")
                        try:
                            genetype_1 = z.strip().split("|")[0]
                            genetype_2 = z.strip().split("|")[1]
                        except:
                            self.set_error("基因型信息传入格式不正确！")
                        w.write("Sample Diff=" + sample1 + "," + str(dep1_min) + "," + str(dep1_max) + "," + genetype_1 + "," +
                                sample2 + "," + str(dep2_min) + "," + str(dep2_max) + "," +
                            genetype_2 + "," + m + "\n")

                group = json.loads(self.option("group"))
                # group = self.option("group").strip().split("|")
                ad = json.loads(self.option("ad"))
                # ad = self.option("ad").strip().split(",")
                miss = json.loads(self.option("max_miss"))
                # miss = self.option("max_miss").strip().split(",")
                maf = json.loads(self.option("maf"))
                # maf = self.option("maf").strip().split(",")
                for (x, y, z, m) in zip(group, ad, miss, maf):
                    try:
                        group_split = x.strip().split(":")[0]
                        group_sample = x.strip().split(":")[1]
                    except BaseException:
                        self.set_error("样本组传入格式不正确！")
                    try:
                        ad_min = y.strip().split("-")[0]
                        ad_max = y.strip().split("-")[1]
                        if ad_min < 0 and not isinstance(ad_min, int):
                            self.set_error("平均深度信息必须为正整数！")
                        if ad_max < 0 and not isinstance(ad_max, int):
                            self.set_error("平均深度信息必须为正整数！")
                        elif ad_max == "":
                            ad_max = 100000
                    except BaseException:
                        self.set_error("平均深度信息传入格式不正确！")
                    if float(z) > 1 or float(z) < 0:
                        self.set_error("最大缺失率必须在0和1之间！")
                    try:
                        maf_min = m.strip().split("-")[0]
                        maf_max = m.strip().split("-")[1]
                        if float(maf_min) > 1 and float(maf_max) < 0:
                            self.set_error("平均频率必须在0和1之间！")
                        if float(maf_max) > 1 or float(maf_max) < 0:
                            self.set_error("平均频率必须在0和1之间！")
                    except BaseException:
                        self.set_error("平均频率传入格式不正确！")
                    w.write("Group Info=" + group_split + "," + str(ad_min) + "," + str(ad_max) + "," + z + "," + str(maf_min) + "," +
                            str(maf_max) + "\n")

    def get_group_path(self):
        """
        group参数是tool的一个参数，这里需要group参数的路径
        :return:
        """
        sample_list = []  # 用于存放用户选择的group文件。
        # if self.option("analysis_type") == "1":
        if not self.option("group") or json.loads(self.option("group")) == []:
            # with open(os.path.join(self.work_dir, "config.txt")) as f:
            #     lines = f.readlines()
            #     for line in lines:
            #         if re.match("Sample Diff=", line):
            #             item = re.split("[=,]", line)
            #             if not item[1] in sample_list:
            #                 sample_list.append(item[1])
            #             if not item[5] in sample_list:
            #                 sample_list.append(item[5])
            if os.path.exists(os.path.join(self.work_dir, "group_new.txt")):
                os.remove(os.path.join(self.work_dir, "group_new"))
            with open(os.path.join(self.work_dir, "group_new"), "w") as w:
                # for i in sample_list:
                #     w.write(i + " " + "default" + "\n")
                w.write("None:None")
            self.group_path = os.path.join(self.work_dir, "group_new")
        else:
            if os.path.exists(os.path.join(self.work_dir, "group_new.txt")):
                os.remove(os.path.join(self.work_dir, "group_new"))
            group1 = {}  # 用于存储group name和 sample的对应值。
            with open(os.path.join(self.work_dir, "group_new"), "w") as w:
                for i in json.loads(self.option("group")):   # 接口中使用
                # for i in self.option("group").strip().split("|"):
                    if i.strip().split(":")[0] not in group1.keys():
                        group1[i.strip().split(":")[0]] = []
                        group1[i.strip().split(":")[0]] = i.strip().split(":")[1]
                for j in group1.keys():
                    sample_list1 = group1[j]
                    # samples = sample_list1.strip().split(",")
                    # for m in samples:
                    #     w.write(m + " " + j + "\n")
                    w.write(j + ":" + sample_list1 + "\n")
            self.group_path = os.path.join(self.work_dir, "group_new")
            return self.group_path

    def run_variant_compare(self):
        if self.option("analysis_model") == "single":
            self.variant_compare_tools = []
            config_file = os.listdir(self.work_dir + "/config_file")
            for config in config_file:
                variant_compare = self.add_tool("wgs_v2.variant_compare")
                variant_compare.set_options(
                    {
                        "filter_recode_vcf": self.option("vcf_file"),
                        "variant_compare_config": os.path.join(self.work_dir + "/config_file", config),
                        "group_table": self.group_path,
                        "name": config,
                    }
                )
                self.variant_compare_tools.append(variant_compare)
            for j in range(len(self.variant_compare_tools)):
                self.variant_compare_tools[j].on("end", self.set_output, "variant_compare_s")
            if self.variant_compare_tools:
                if len(self.variant_compare_tools) > 1:
                    self.on_rely(self.variant_compare_tools, self.end)
                elif len(self.variant_compare_tools) == 1:
                    self.variant_compare_tools[0].on('end', self.end)
            else:
                self.set_error("variant_compare_tools为空！")
            for tool in self.variant_compare_tools:
                gevent.sleep(1)
                tool.run()

        else:
            self.variant_compare = self.add_tool("wgs_v2.variant_compare")
            options = {
                "filter_recode_vcf": self.option("vcf_file"),
                "group_table": self.group_path,
                "variant_compare_config": os.path.join(
                    self.work_dir,
                    "diff.config"),
                "name": "variant_compare",
            }
            self.variant_compare.set_options(options)
            self.variant_compare.on("end", self.set_output, "variant_compare_m")
            self.variant_compare.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'variant_compare_s':
            self.linkdir(obj.output_dir, self.output_dir + "/variant_compare")
        elif event['data'] == 'variant_compare_m':
            self.linkdir(obj.output_dir, self.output_dir + "/variant_compare")
            self.end()

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        self.logger.info("保存结果到mongo")
        api_path = self.api.api("wgs_v2.variant_compare")
        if self.option("analysis_model") == "single":
            list_vcf = []
            vcf_path = self._sheet.output + "/variant_compare/"
            path = os.listdir(self.output_dir+"/variant_compare/")
            for i in path:
                name = i.strip().split(".")[0]
                if re.match(".*detail$", i):
                    api_path.sg_varian_compare_detail(self.option("main_id"), self.output_dir + '/variant_compare/' +
                                                        i, name)
                elif re.match(".*filter\.vcf$", i):
                    # api_path.add_sg_variant_compare_stat(self.option("main_id"), self.output_dir + '/variant_compare/' +
                    #                                      i,  name)
                    list_vcf.append(vcf_path + i)
                elif re.match(".*snp_indel_stat\.txt$", i):
                    api_path.add_sg_variant_compare_stat_v2(self.option("main_id"),
                                                            self.output_dir + '/variant_compare/' + i, name)
                elif re.match(".*eff$", i):
                    api_path.add_sg_variant_compare_effect(self.option("main_id"), self.output_dir + '/variant_compare/'
                                                           + i,  name)
                    if self.option("variant_type") == "SNP,INDEL":
                        api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
                                                               self.output_dir + '/variant_compare/' + i, name,
                                                               "all")
                        api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
                                                                self.output_dir + '/variant_compare/' + i, name,
                                                                "snp")
                        api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
                                                                self.output_dir + '/variant_compare/' + i, name,
                                                                "indel")
                    else:
                        api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
                                                                self.output_dir + '/variant_compare/' + i, name,
                                                                self.option("variant_type").lower())

                elif re.match(".*func$", i):
                    api_path.add_sg_variant_compare_func(self.option("main_id"), self.output_dir + '/variant_compare/'
                                                           + i, name)
                    if self.option("variant_type") == "SNP,INDEL":
                        api_path.sg_variant_compare_impact_bar(self.option("main_id"), self.option("task_id"),
                                                               self.output_dir + '/variant_compare/' + i, name, "all")
                        api_path.sg_variant_compare_impact_bar(self.option("main_id"), self.option("task_id"),
                                                               self.output_dir + '/variant_compare/' + i, name, "snp")
                        api_path.sg_variant_compare_impact_bar(self.option("main_id"), self.option("task_id"),
                                                               self.output_dir + '/variant_compare/' + i, name, "indel")
                    else:
                        api_path.sg_variant_compare_impact_bar(self.option("main_id"), self.option("task_id"),
                                                               self.output_dir + '/variant_compare/' + i, name,
                                                               self.option("variant_type").lower())
                else:
                    self.set_error("结果文件中存在着不明文件")
            api_path.update_variant_compare(list_vcf, self.option("main_id"))
        else:
            api_path.sg_varian_compare_detail(self.option("main_id"),
                                              self.output_dir + '/variant_compare/variant_compare.detail',
                                              "variant_compare")
            api_path.update_variant_compare(self.output_dir+'/variant_compare/variant_compare.filter.vcf',
                                            self.option("main_id"))
            if os.path.exists(self.output_dir + '/variant_compare/variant_compare.eff'):
                api_path.add_sg_variant_compare_effect(self.option("main_id"), self.output_dir +
                                                   '/variant_compare/variant_compare.eff', "variant_compare")
            if os.path.exists(self.output_dir + '/variant_compare/variant_compare.func'):
                api_path.add_sg_variant_compare_func(self.option("main_id"), self.output_dir +
                                                 '/variant_compare/variant_compare.func', "variant_compare")
            if os.path.exists(self.output_dir + '/variant_compare/variant_compare.snp_indel_stat.txt'):
                api_path.add_sg_variant_compare_stat_v2(self.option("main_id"),
                                                        self.output_dir +
                                                        '/variant_compare/variant_compare.snp_indel_stat.txt',
                                                        "variant_compare")
            else:
                api_path.add_sg_variant_compare_stat(self.option("main_id"),
                                                     self.output_dir + '/variant_compare/variant_compare.filter.vcf',
                                                     "variant_compare")
            # api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
            #                                        self.output_dir + '/variant_compare/variant_compare.eff',
            #                                        "variant_compare", "all" if self.option("variant_type") ==
            #                                                         "SNP,INDEL" else self.option("variant_type"))
            if self.option("variant_type") == "SNP,INDEL":
                if os.path.exists(self.output_dir + '/variant_compare/variant_compare.func'):
                    api_path.sg_variant_compare_impact_bar(self.option("main_id"), self.option("task_id"),
                                                           self.output_dir + '/variant_compare/variant_compare.func',
                                                           "variant_compare", "all")
                    api_path.sg_variant_compare_impact_bar(self.option("main_id"), self.option("task_id"),
                                                           self.output_dir + '/variant_compare/variant_compare.func',
                                                           "variant_compare", "snp")
                    api_path.sg_variant_compare_impact_bar(self.option("main_id"), self.option("task_id"),
                                                       self.output_dir + '/variant_compare/variant_compare.func',
                                                       "variant_compare", "indel")
                if os.path.exists(self.output_dir + '/variant_compare/variant_compare.eff'):
                    api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
                                                            self.output_dir + '/variant_compare/variant_compare.eff',
                                                            "variant_compare", "all")
                    api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
                                                        self.output_dir + '/variant_compare/variant_compare.eff',
                                                        "variant_compare", "snp")
                    api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
                                                        self.output_dir + '/variant_compare/variant_compare.eff',
                                                        "variant_compare", "indel")
            else:
                if os.path.exists(self.output_dir + '/variant_compare/variant_compare.func'):
                    api_path.sg_variant_compare_impact_bar(self.option("main_id"), self.option("task_id"),
                                                       self.output_dir + '/variant_compare/variant_compare.func',
                                                       "variant_compare",
                                                       self.option("variant_type").lower())
                if os.path.exists(self.output_dir + '/variant_compare/variant_compare.eff'):
                    api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
                                                    self.output_dir + '/variant_compare/variant_compare.eff',
                                                    "variant_compare",
                                                    self.option("variant_type").lower())

            # api_path.sg_variant_compare_impact_bar(self.option("main_id"), self.option("task_id"),
            #                                        self.output_dir + '/variant_compare/variant_compare.func',
            #                                        "variant_compare", "all" if
            #                                        self.option("variant_type") == "SNP,INDEL" else
            #                                        self.option("variant_type"))
            # if self.option("variant_type") == "SNP,INDEL":
            #     api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
            #                                            self.output_dir + '/variant_compare/variant_compare.eff',
            #                                            "variant_compare", "all")
            #     api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
            #                                            self.output_dir + '/variant_compare/variant_compare.eff',
            #                                            "variant_compare", "snp")
            #     api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
            #                                            self.output_dir + '/variant_compare/variant_compare.eff',
            #                                            "variant_compare", "indel")
            # else:
            #     api_path.add_variant_compare_effect_bar(self.option("main_id"), self.option("task_id"),
            #                                            self.output_dir + '/variant_compare/variant_compare.eff',
            #                                            "variant_compare",
            #                                            self.option("variant_type").lower())

    def file_check(self):    # 检查是否有eff或者func# 文件，在tool中，vcf为空的时候没有这两个表的结果，也就不需要导表了
        files = os.listdir(self.output_dir + "/variant_compare/")
        num = 1  # eff不存在
        for j in files:
            if re.match(".*eff$", j):
                num = 0
                break
        if num == 0:
            return True
        elif num ==1:
            return False

    def run(self):
        self.get_config()
        self.get_group_path()
        self.run_variant_compare()
        super(VariantCompareWorkflow, self).run()

    def end(self):
        gevent.sleep(1)
        # if self.file_check():
        self.set_db()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(VariantCompareWorkflow, self).end()
