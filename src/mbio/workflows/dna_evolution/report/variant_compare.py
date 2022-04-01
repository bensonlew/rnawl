# -*- coding: utf-8 -*-
# __author__ = 'Binbin Zhao'
# modified 20180823

import os
import json
import re
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
            {"name": "analysis_type", "type": "string"},  # 分析类型，用于区分是样本比较还是样本组比较
            {"name": "variant_type", "type": "string",
                "default": "all"},  # 变异类型，用于区分是snp、indel还是all
            # 这里所有需要default值的都需要下面在判断的时候一个个列出，而且需要考虑list当中有的地方是空值，有的地方不空的情况。
            {"name": "is_same", "type": "string"},
            {"name": "vcf_file", "type": "string"},
            {"name": "efftype", "type": "string"},
            # 3_prime_UTR_variant, stop_lost 这个默认是什么。
            {"name": "funtype", "type": "string"}, # High, Low等
            {"name": "dep", "type": "string"},
            {"name": "maf", "type": "string"},
            {"name": "ad", "type": "string"},
            {"name": "miss", "type": "string"},
            {"name": "location", "type": "string"},
            {"name": "update_info", "type": "string"},
            # {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "region", "type": "string"},
            {"name": "project_sn", "type": "string"}
        ]
        self.variant_type = " "  # 用于存储变异类型信息。
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.var_diff = self.add_tool("dna_evolution.variant_compare")
        # self.chromosome = self.add_tool("dna_evolution.chromosome_window")
        self.group_path = ""

    def check_options(self):
        if self.option('analysis_type'):
            if self.option('analysis_type') not in ["1", "2", "3"]:
                raise OptionError("分析方法类型不合法！必须为1 or 2 or 3", code="14900201")
        else:
            raise OptionError('必须提供analysis_type结果表', code="14900202")
        if not self.option("vcf_file"):
            raise OptionError("缺少vcf_file参数", code="14900203")
        if self.option('variant_type'):
            if self.option('variant_type') not in ["snp", "indel", "all"]:
                raise OptionError("变异类型不合法！必须为snp、indel和all", code="14900204")
            else:
                self.variant_type = "SNP,INDEL" if self.option(
                    "variant_type") == "all" else self.option("variant_type").upper()
        if self.option('analysis_type') == "1":
            if not self.option("sample"):
                raise OptionError("缺少sample参数", code="14900205")
        if self.option('analysis_type') == "2":
            if not self.option("group"):
                raise OptionError("缺少group参数", code="14900206")
        if self.option('analysis_type') == "3":
            if not self.option("group"):
                raise OptionError("缺少group参数", code="14900207")
        return True

    def get_config(self):
        """
        将tool所有使用的参数包成一个config文件传给Tool。
        这里还需要补充染色体的信息。
        :return:
        """
        if os.path.exists(os.path.join(self.work_dir, "config.txt")):
            os.remove(os.path.join(self.work_dir, "config.txt"))
        with open(os.path.join(self.work_dir, "config.txt"), "w") as w:
            w.write("Variant Type=" + self.variant_type + "\n")
            w.write("Variant Eff=" + self.option("funtype").upper() + "\n")
            w.write("Variant Ann=" + self.option("efftype") + "\n")
            if self.option("region") is not None:
                w.write("Region=" + self.option("region") + "\n")
            if self.option("analysis_type") == "3":
                sample = json.loads(self.option("sample"))
                dep = json.loads(self.option("dep"))
                genotype = json.loads(self.option("genotype"))
                is_same = json.loads(self.option("is_same"))
                for (x, y, z, m)in zip(sample, dep, genotype, is_same):
                    try:
                        sample1 = x.strip().split("|")[0]
                        sample2 = x.strip().split("|")[1]
                    except BaseException:
                        self.set_error("样本传入格式不正确！", code="14900201")
                        raise Exception("样本传入格式不正确！")
                    try:
                        dep1_min = y.strip().split(
                            "|")[0].strip().split("-")[0]
                        dep1_max = y.strip().split(
                            "|")[0].strip().split("-")[1]
                        dep2_min = y.strip().split(
                            "|")[1].strip().split("-")[0]
                        dep2_max = y.strip().split(
                            "|")[1].strip().split("-")[1]
                        if dep1_min == "":
                            dep1_min = 5
                        elif float(dep1_min) < 0:
                            self.set_error("深度信息必须为正整数！", code="14900202")
                            raise Exception("深度信息必须为正整数！")
                        if dep1_max == "":
                            dep1_max = 100000
                        elif float(dep1_min) < 0 :
                            self.set_error("深度信息必须为正整数！", code="14900203")
                            raise Exception("深度信息必须为正整数！")
                        if dep2_min == "":
                            dep1_min = 5
                        elif float(dep2_min) < 0:
                            self.set_error("深度信息必须为正整数！", code="14900204")
                            raise Exception("深度信息必须为正整数！")
                        if dep2_max == "":
                            dep2_max = 100000
                        elif float(dep2_min) < 0:
                            self.set_error("深度信息必须为正整数！", code="14900205")
                            raise Exception("深度信息必须为正整数！")
                    except BaseException:
                        self.set_error("测序深度信息传入格式不正确！", code="14900206")
                        raise Exception("测序深度信息传入格式不正确！")
                    w.write(
                        "Sample Diff=" +
                        sample1 +
                        "," +
                        str(dep1_min) +
                        "," +
                        str(dep1_max) +
                        "," +
                        z.strip().split("|")[0] +
                        "," +
                        sample2 +
                        "," +
                        str(dep2_min) +
                        "," +
                        str(dep2_max) +
                        "," +
                        z.strip().split("|")[1] +
                        "," +
                        m +
                        "\n")
                group = json.loads(self.option("group"))
                ad = json.loads(self.option("ad"))
                miss = json.loads(self.option("miss"))
                maf = json.loads(self.option("maf"))
                for (x, y, z, m) in zip(group, ad, miss, maf):
                    try:
                        group_split = x.strip().split(":")[0]
                    except BaseException:
                        self.set_error("样本组传入格式不正确！", code="14900207")
                        raise Exception("样本组传入格式不正确！")
                    try:
                        ad_min = y.strip().split("-")[0]
                        ad_max = y.strip().split("-")[1]
                        if ad_min == "":
                            ad_min = 5
                        elif ad_min < 0 and not isinstance(ad_min, int):
                            self.set_error("平均深度信息必须为正整数！", code="14900208")
                            raise Exception("平均深度信息必须为正整数！")
                        if ad_max == "":
                            ad_max = 100000
                        elif ad_max < 0 and not isinstance(ad_max, int):
                            self.set_error("平均深度信息必须为正整数！", code="14900209")
                            raise Exception("平均深度信息必须为正整数！")
                    except BaseException:
                        self.set_error("平均深度信息传入格式不正确！", code="14900210")
                        raise Exception("平均深度信息传入格式不正确！")
                    try:
                        miss_min = z.strip().split("-")[0]
                        miss_max = z.strip().split("-")[1]
                        if miss_min == "":
                            miss_min = 0.05
                        elif float(miss_min) > 1 or float(miss_min) < 0:
                            self.set_error("缺失率必须在0和1之间！", code="14900211")
                            raise Exception("缺失率必须在0和1之间！")
                        if miss_max == "":
                            miss_max = 1
                        elif float(miss_max) > 1 or float(miss_min) < 0:
                            self.set_error("缺失率必须在0和1之间！", code="14900212")
                            raise Exception("缺失率必须在0和1之间！")
                    except BaseException:
                        self.set_error("缺失率传入格式不正确！", code="14900213")
                        raise Exception("缺失率传入格式不正确！")
                    try:
                        maf_min = m.strip().split("-")[0]
                        maf_max = m.strip().split("-")[1]
                        if maf_min == "":
                            maf_min = 0.3
                        elif float(maf_min) > 1 and float(maf_max) < 0:
                            self.set_error("平均频率必须在0和1之间！", code="14900214")
                            raise Exception("平均频率必须在0和1之间！")
                        if maf_max == "":
                            maf_max = 0.3
                        elif float(maf_max) > 1 or float(maf_max) < 0:
                            self.set_error("平均频率必须在0和1之间！", code="14900215")
                            raise Exception("平均频率必须在0和1之间！")
                    except BaseException:
                        self.set_error("平均频率传入格式不正确！", code="14900216")
                        raise Exception("平均频率传入格式不正确！")
                    w.write(
                        "Group Info=" +
                        group_split +
                        "," +
                        str(ad_min) +
                        "," +
                        str(ad_max) +
                        "," +
                        str(miss_min) +
                        "," +
                        str(miss_max) +
                        "," +
                        str(maf_min) +
                        "," +
                        str(maf_max) +
                        "\n")
            elif self.option("analysis_type") == "2":
                group = json.loads(self.option("group"))
                ad = json.loads(self.option("ad"))
                miss = json.loads(self.option("miss"))
                maf = json.loads(self.option("maf"))
                for (x, y, z, m) in zip(group, ad, miss, maf):
                    try:
                        group_split = x.strip().split(":")[0]
                    except BaseException:
                        self.set_error("样本组传入格式不正确！", code="14900217")
                        raise Exception("样本组传入格式不正确！")
                    try:
                        ad_min = y.strip().split("-")[0]
                        ad_max = y.strip().split("-")[1]
                        if ad_min == "":
                            ad_min = 5
                        elif ad_min < 0 and not isinstance(ad_min, int):
                            self.set_error("平均深度信息必须为正整数！", code="14900218")
                            raise Exception("平均深度信息必须为正整数！")
                        if ad_max == "":
                            ad_max = 100000
                        elif ad_max < 0 and not isinstance(ad_max, int):
                            self.set_error("平均深度信息必须为正整数！", code="14900219")
                            raise Exception("平均深度信息必须为正整数！")
                    except BaseException:
                        self.set_error("平均深度信息传入格式不正确！", code="14900220")
                        raise Exception("平均深度信息传入格式不正确！")
                    try:
                        miss_min = z.strip().split("-")[0]
                        miss_max = z.strip().split("-")[1]
                        if miss_min == "":
                            miss_min = 0.05
                        elif float(miss_min) > 1 or float(miss_min) < 0:
                            self.set_error("缺失率必须在0和1之间！", code="14900221")
                            raise Exception("缺失率必须在0和1之间！")
                        if miss_max == "":
                            miss_min = 1
                        elif float(miss_max) > 1 and float(miss_min) < 0:
                            self.set_error("缺失率必须在0和1之间！", code="14900222")
                            raise Exception("缺失率必须在0和1之间！")
                    except BaseException:
                        self.set_error("缺失率传入格式不正确！", code="14900223")
                        raise Exception("缺失率传入格式不正确！")
                    try:
                        maf_min = m.strip().split("-")[0]
                        maf_max = m.strip().split("-")[1]
                        if maf_min == "":
                            maf_min = 0.3
                        elif float(maf_min) > 1 and float(maf_max) < 0:
                            self.set_error("平均频率必须在0和1之间！", code="14900224")
                            raise Exception("平均频率必须在0和1之间！")
                        if maf_max == "":
                            maf_max = 0.3
                        elif float(maf_max) > 1 or float(maf_max) < 0:
                            self.set_error("平均频率必须在0和1之间！", code="14900225")
                            raise Exception("平均频率必须在0和1之间！")
                    except BaseException:
                        self.set_error("平均频率传入格式不正确！", code="14900226")
                        raise Exception("平均频率传入格式不正确！")
                    w.write(
                        "Group Info=" +
                        group_split +
                        "," +
                        str(ad_min) +
                        "," +
                        str(ad_max) +
                        "," +
                        str(miss_min) +
                        "," +
                        str(miss_max) +
                        "," +
                        str(maf_min) +
                        "," +
                        str(maf_max) +
                        "\n")
            elif self.option("analysis_type") == "1":
                sample = json.loads(self.option("sample"))
                dep = json.loads(self.option("dep"))
                genotype = json.loads(self.option("genotype"))
                is_same = json.loads(self.option("is_same"))
                for (x, y, z, m) in zip(sample, dep, genotype, is_same):
                    try:
                        sample1 = x.strip().split("|")[0]
                        sample2 = x.strip().split("|")[1]
                    except BaseException:
                        self.set_error("样本传入格式不正确！", code="14900227")
                        raise Exception("样本传入格式不正确！")
                    try:
                        dep1_min = y.strip().split(
                            "|")[0].strip().split("-")[0]
                        dep1_max = y.strip().split(
                            "|")[0].strip().split("-")[1]
                        dep2_min = y.strip().split(
                            "|")[1].strip().split("-")[0]
                        dep2_max = y.strip().split(
                            "|")[1].strip().split("-")[1]
                        if dep1_min == "":
                            dep1_min = 5
                        elif float(dep1_min) < 0 and not isinstance(dep1_min, int):
                            self.set_error("深度信息必须为正整数！", code="14900228")
                            raise Exception("深度信息必须为正整数！")
                        if dep1_max == "":
                            dep1_max = 100000
                        elif dep1_min < 0 and not isinstance(dep1_min, int):
                            self.set_error("深度信息必须为正整数！", code="14900229")
                            raise Exception("深度信息必须为正整数！")
                        if dep2_min == "":
                            dep1_min = 5
                        elif dep2_min < 0 and not isinstance(dep2_min, int):
                            self.set_error("深度信息必须为正整数！", code="14900230")
                            raise Exception("深度信息必须为正整数！")
                        if dep2_max == "":
                            dep2_max = 100000
                        elif dep2_min < 0 and not isinstance(dep2_min, int):
                            self.set_error("深度信息必须为正整数！", code="14900231")
                            raise Exception("深度信息必须为正整数！")
                    except BaseException:
                        self.set_error("测序深度信息传入格式不正确！", code="14900232")
                        raise Exception("测序深度信息传入格式不正确！")
                    w.write(
                        "Sample Diff=" +
                        sample1 +
                        "," +
                        str(dep1_min) +
                        "," +
                        str(dep1_max) +
                        "," +
                        z.strip().split("|")[0] +
                        "," +
                        sample2 +
                        "," +
                        str(dep2_min) +
                        "," +
                        str(dep2_max) +
                        "," +
                        z.strip().split("|")[1] +
                        "," +
                        m +
                        "\n")

    def run_variant_compare(self):
        self.logger.info("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
        self.logger.info(self.option("project_sn"))
        options = {
            "filter_recode_vcf": self.option("vcf_file"),
            "group_table": self.group_path,
            "variant_compare_config": os.path.join(
                self.work_dir,
                "config.txt"),
        }
        self.var_diff.set_options(options)
        # self.var_diff.on("end", self.run_chromosome_window)
        self.var_diff.on("end", self.set_output, "variant_compare")
        self.var_diff.run()

    # def run_chromosome_window(self):
    #     api_path1 = self.api.api("dna_evolution.chromosome_window")
    #     params_json = {"step_num": 10,
    #                    "variant_type": "all",
    #                    "project_sn": self.option("project_sn"),
    #                    "task_id": self.option("task_id")}
    #     # params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
    #     main_id = api_path1.add_sg_chromosome_window(self.option("project_sn"), params_json, self.option("task_id"),
    #                                                  self.option("main_id"))
    #     pop_table_path = os.path.join(self.var_diff.output_dir, "pop.table")
    #     self.logger.info(pop_table_path)
    #     options = {
    #         "pop_table": pop_table_path,
    #         "step_num": 10,
    #         "variant_type": "all",
    #         "main_id": str(main_id),
    #     }
    #     self.chromosome.set_options(options)
    #     self.chromosome.on("end", self.set_output, "chromosome_window")
    #     # self.chromosome.on("end", self.end)
    #     self.chromosome.run()

    def get_group_path(self):
        """
        group参数是tool的一个参数，这里需要group参数的路径
        :return:
        """
        sample_list = []  # 用于存放用户选择的group文件。
        if self.option("analysis_type") == "1":
            with open(os.path.join(self.work_dir, "config.txt")) as f:
                lines = f.readlines()
                for line in lines:
                    if re.match("Sample Diff=", line):
                        item = re.split("[=,]", line)
                        if not item[1] in sample_list:
                            sample_list.append(item[1])
                        if not item[5] in sample_list:
                            sample_list.append(item[5])
            if os.path.exists(os.path.join(self.work_dir, "group_new.txt")):
                os.remove(os.path.join(self.work_dir, "group_new"))
            with open(os.path.join(self.work_dir, "group_new"), "w") as w:
                for i in sample_list:
                    w.write(i + " " + "default" + "\n")
            self.group_path = os.path.join(self.work_dir, "group_new")
        else:
            if os.path.exists(os.path.join(self.work_dir, "group_new.txt")):
                os.remove(os.path.join(self.work_dir, "group_new"))
            group1 = {}  # 用于存储group name和 sample的对应值。
            with open(os.path.join(self.work_dir, "group_new"), "w") as w:
                for i in json.loads(self.option("group")):
                    if i.strip().split(":")[0] not in group1.keys():
                        group1[i.strip().split(":")[0]] = []
                        group1[i.strip().split(":")[0]
                               ] = i.strip().split(":")[1]
                for j in group1.keys():
                    sample_list1 = group1[j]
                    samples = sample_list1.strip().split(",")
                    for m in samples:
                        w.write(m + " " + j + "\n")
            self.group_path = os.path.join(self.work_dir, "group_new")
            return self.group_path

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'variant_compare':
            self.linkdir(obj.output_dir, 'variant_compare')
        # if event['data'] == 'chromosome_window':
        #     self.linkdir(obj.output_dir, 'chromosome_window')
            self.end()

    def vcf_none(self, vcf_path):
        with open(vcf_path, "r") as f:
            num = 0
            for lines in f:
                if not re.match('#', lines):  # 检查vcf不为空
                    num = 1 # 代表的是不为空
                    break
            if num == 0:
                return False
            else:
                return True

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
        # 导表的主表名字为variant_compare
        api_path = self.api.api("dna_evolution.variant_compare")
        chromosome_window = self.api.api("dna_evolution.chromosome_window")
        params_json = {"step_num": 10,
                       "variant_type": "all",
                       "project_sn": self.option("project_sn"),
                       "task_id": self.option("task_id")}
        main_id = chromosome_window.add_sg_chromosome_window(self.option("project_sn"), params_json, self.option("task_id"),
                                                     self.option("main_id"))
        chromosome_window.add_sg_distribution(self.output_dir+ "/variant_compare/", main_id, "all")
        download_path = self._sheet.output + "/variant_compare/pop.table"
        filter_vcf_path = self._sheet.output + "/variant_compare/pop.filtered.vcf"
        api_path.update_varaint_compare(download_path, filter_vcf_path, self.option("main_id"))
        api_path.add_sg_variant_compare_effect(
            self.option("main_id"),
            self.var_diff.output_dir + "/eff.type")
        api_path.add_variant_compare_effect_bar(
            self.option("main_id"),
            self.option("task_id"),
            self.var_diff.output_dir + "/eff.type",
            self.option("variant_type"))
        api_path.sg_variant_compare_impact(
            self.option("main_id"),
            self.var_diff.output_dir +
            "/function.type")
        api_path.sg_variant_compare_impact_bar(
            self.option("main_id"),
            self.option("task_id"),
            self.var_diff.output_dir +
            "/function.type",
            self.option("variant_type"))
        api_path.sg_varian_compare_detail(
            self.option("main_id"),
            self.var_diff.output_dir +
            "/pop.table",
            self.option("variant_type"))
        if self.option("variant_type") == "all":
            api_path.add_variant_compare_effect_bar(
                self.option("main_id"),
                self.option("task_id"),
                self.var_diff.output_dir + "/eff.type",
                "snp")
            api_path.add_variant_compare_effect_bar(
                self.option("main_id"),
                self.option("task_id"),
                self.var_diff.output_dir + "/eff.type",
                "indel")
            api_path.sg_variant_compare_impact_bar(
                self.option("main_id"),
                self.option("task_id"),
                self.var_diff.output_dir +
                "/function.type",
                "snp")
            api_path.sg_variant_compare_impact_bar(
                self.option("main_id"),
                self.option("task_id"),
                self.var_diff.output_dir +
                "/function.type",
                "indel")
            api_path.sg_varian_compare_detail(
                self.option("main_id"),
                self.var_diff.output_dir +
                "/pop.table",
                "snp")
            api_path.sg_varian_compare_detail(
                self.option("main_id"),
                self.var_diff.output_dir +
                "/pop.table",
                "indel")

    def run(self):
        self.get_config()
        self.get_group_path()
        self.run_variant_compare()
        super(VariantCompareWorkflow, self).run()

    def end(self):
        if self.vcf_none(self.output_dir + "/variant_compare/pop.filtered.vcf"):
            self.set_db()
        else:
            chromosome_window = self.api.api("dna_evolution.chromosome_window")
            self.logger.info("***********************************************************")
            self.logger.info(self.option("project_sn"))
            params_json = {"step_num": 10,
                           "variant_type": "all",
                           "project_sn": self.option("project_sn"),
                           "task_id": self.option("task_id")}
            main_id = chromosome_window.add_sg_chromosome_window(self.option("project_sn"), params_json,
                                                                 self.option("task_id"),
                                                                 self.option("main_id"))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(VariantCompareWorkflow, self).end()
