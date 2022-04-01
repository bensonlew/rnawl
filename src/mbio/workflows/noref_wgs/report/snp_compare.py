# -*- coding: utf-8 -*-
# __author__ = 'Binbin Zhao'
# modified 20190110

import os
import json
import re
import gevent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class SnpCompareWorkflow(Workflow):
    """
    交互分析：样本比较分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SnpCompareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "sample", "type": "string"},
            {"name": "genotype", "type": "string"},
            {"name": "group", "type": "string"},  # group用于存放所有的group信息。
            {"name": "analysis_type", "type": "string"},  # 分析类型，用于区分是样本比较还是样本组比较
            # 这里所有需要default值的都需要下面在判断的时候一个个列出，而且需要考虑list当中有的地方是空值，有的地方不空的情况。
            {"name": "marktype", "type": "string"}, # 是diff还是same
            {"name": "vcf_file", "type": "string"},
            {"name": "dep", "type": "string"},
            {"name": "maf", "type": "string"},
            {"name": "ad", "type": "string"},     # 为什么ad在下面的时候变成int格式了
            {"name": "max_miss", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "tag_file", "type": "string"},  # 为population.tag文件
            {"name": "task_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_sn", "type": "string"},
            {"name": "analysis_model", "type": "string"},  # 用于判断分析的是但条件还是多条件, 为single或者multiple
            {"name": "alle_number", "type": "string"},  # 用于判断分析的是但条件还是多条件, 为single或者multiple
          ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.group_path = ""

    def check_options(self):
        if self.option('analysis_model') not in ["single", "multiple"]:
            raise OptionError("分析模式不合法！必须为single和multiple其中之一。", code="15500107")
        if self.option('analysis_model') == "multiple":
            if self.option('analysis_type'):
                if self.option('analysis_type') not in ["1", "2", "3"]:
                    raise OptionError("分析方法类型不合法！必须为1 or 2 or 3", code="15500108")
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
                    self.set_error("样本传入格式不正确！", code="15500137")
                with open(os.path.join(self.work_dir + "/config_file", sample1 + "_vs_" + sample2 + "_compare"), "w") as w:
                    w.write("{" + "\n")
                    w.write("Alle_Number:" + self.option("alle_number") + "\n")
                    w.write("Model:Single" + "\n")
                    w.write("Selecttion:" + "\n")
                    w.write("\t" + "mode:sample" + "\n")
                    w.write("\t" + "genetype:" + self.option("genotype") + "\n")
                    min_dep = self.option("dep").strip().split(",")[0]
                    max_dep = self.option("dep").strip().split(",")[1]
                    if max_dep == "":
                        max_dep = 100000
                    w.write("\t" + "depth:" + str(min_dep) + "," + str(max_dep) + "\n")
                    w.write("\t" + "marktype:" + self.option("marktype") + "\n")
                    self.logger.info("#############################")
                    self.logger.info("\t" + "Sample:" + sample1 + "," + sample2 + "\n")
                    self.logger.info("#####################%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                    w.write("\t" + "Sample:" + sample1 + "," + sample2 + "\n")
                    w.write("}" + "\n")

        else:
            with open(os.path.join(self.work_dir, "diff.config"), "w") as w:
                w.write("{" + "\n")
                w.write("Alle_Number:" + self.option("alle_number") + "\n")
                w.write("Model:Multiple" + "\n")
                w.write("Selecttion:" + "\n")
                sample = json.loads(self.option("sample"))
                if not sample == [""]:
                    dep = json.loads(self.option("dep"))
                    genotype = json.loads(self.option("genotype"))
                    is_same = json.loads(self.option("marktype"))
                    for (x, y, z, m)in zip(sample, dep, genotype, is_same):
                        try:
                            sample1 = x.strip().split("|")[0]
                            sample2 = x.strip().split("|")[1]
                        except BaseException:
                            self.set_error("样本传入格式不正确！", code="15500138")
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
                                self.set_error("深度信息必须为正整数！", code="15500139")
                            if dep1_max == "":
                                dep1_max = 100000
                            elif float(dep1_max) < 0:
                                self.set_error("深度信息必须为正整数！", code="15500140")
                            if float(dep2_min) < 0:
                                self.set_error("深度信息必须为正整数！", code="15500141")
                            if dep2_max == "":
                                dep2_max = 100000
                            elif float(dep2_max) < 0:
                                self.set_error("深度信息必须为正整数！", code="15500142")
                        except BaseException:
                            self.set_error("测序深度信息传入格式不正确！", code="15500143")
                        try:
                            genetype_1 = z.strip().split("|")[0]
                            genetype_2 = z.strip().split("|")[1]
                        except:
                            self.set_error("基因型信息传入格式不正确！", code="15500144")
                        w.write("sampleinfo:" + "\n")
                        w.write("\t\t" + "sample1_name:" + sample1 + "\n")
                        w.write("\t\t" + "sample1_type:" +  genetype_1 + "\n")
                        w.write("\t\t" + "sample1_depth:" + str(dep1_min) + "," + str(dep1_max) + "\n")
                        w.write("\t\t" + "sample2_name:" + sample2 + "\n")
                        w.write("\t\t" + "sample2_type:" + genetype_2 + "\n")
                        w.write("\t\t" + "sample2_depth:" + str(dep2_min) + "," + str(dep2_max) + "\n")
                        w.write("\t\t" + "sampleresult:" + m + "\n")

                group = json.loads(self.option("group"))
                print group
                ad = json.loads(self.option("ad"))
                miss = json.loads(self.option("max_miss"))
                maf = json.loads(self.option("maf"))
                for (x, y, z, m) in zip(group, ad, miss, maf):
                    try:
                        group_split = x.strip().split(":")[0]
                        print group_split
                        group_sample = x.strip().split(":")[1]
                        print group_sample
                    except BaseException:
                        self.set_error("样本组传入格式不正确！", code="15500145")
                    try:
                        ad_min = y.strip().split("-")[0]
                        ad_max = y.strip().split("-")[1]
                        if ad_min < 0 and not isinstance(ad_min, int):
                            self.set_error("平均深度信息必须为正整数！", code="15500146")
                        if ad_max < 0 and not isinstance(ad_max, int):
                            self.set_error("平均深度信息必须为正整数！", code="15500147")
                        elif ad_max == "":
                            ad_max = 100000
                    except BaseException:
                        self.set_error("平均深度信息传入格式不正确！", code="15500148")
                    if float(z) > 1 or float(z) < 0:
                        self.set_error("最大缺失率必须在0和1之间！", code="15500149")
                    try:
                        maf_min = m.strip().split("-")[0]
                        maf_max = m.strip().split("-")[1]
                        if float(maf_min) > 1 and float(maf_max) < 0:
                            self.set_error("平均频率必须在0和1之间！", code="15500150")
                        if float(maf_max) > 1 or float(maf_max) < 0:
                            self.set_error("平均频率必须在0和1之间！", code="15500151")
                    except BaseException:
                        self.set_error("平均频率传入格式不正确！", code="15500152")
                    w.write("groupinfo:" + "\n")
                    w.write("\t\t" + "group_name:" + group_split + "\n")
                    w.write("\t\t" + "group_sample:" + group_sample + "\n")
                    w.write("\t\t" + "group_depth:" + str(ad_min) + "," + str(ad_max) + "\n")
                    w.write("\t\t" + "group_af:" + maf_min + "," + maf_max +  "\n")
                    w.write("\t\t" + "group_miss:" + z + "\n")

    def run_snp_compare(self):
        if self.option("analysis_model") == "single":
            self.snp_compare_tools = []
            config_file = os.listdir(self.work_dir + "/config_file")
            for config in config_file:
                snp_compare = self.add_tool("noref_wgs.snp_compare")
                snp_compare.set_options(
                    {
                        "vcf_file": self.option("vcf_file"),
                        "tag_file": self.option("tag_file"),
                        "analysis_name": config,
                        "config_file": os.path.join(self.work_dir + "/config_file", config),
                        "analysis_model": self.option("analysis_model")
                    }
                )
                self.snp_compare_tools.append(snp_compare)
            for j in range(len(self.snp_compare_tools)):
                self.snp_compare_tools[j].on("end", self.set_output, "snp_compare_s")
            if self.snp_compare_tools:
                if len(self.snp_compare_tools) > 1:
                    self.on_rely(self.snp_compare_tools, self.end)
                elif len(self.snp_compare_tools) == 1:
                    self.snp_compare_tools[0].on('end', self.end)
            else:
                self.set_error("snp_compare_tools为空！", code="15500153")
            for tool in self.snp_compare_tools:
                self.logger.info("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                self.logger.info(tool)
                gevent.sleep(1)
                tool.run()

        else:

            self.snp_compare = self.add_tool("noref_wgs.snp_compare")
            options = {
                "vcf_file": self.option("vcf_file"),
                "tag_file": self.option("tag_file"),
                "analysis_name": "snp_compare",
                "config_file": os.path.join(
                    self.work_dir,
                    "diff.config"),
                "analysis_model": self.option("analysis_model")
            }
            self.snp_compare.set_options(options)
            self.snp_compare.on("end", self.set_output, "snp_compare_m")
            self.snp_compare.run()


    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'snp_compare_s':
            self.linkdir(obj.output_dir, self.output_dir + "/snp_compare")
        elif event['data'] == 'snp_compare_m':
            self.linkdir(obj.output_dir, self.output_dir + "/snp_compare")
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
        api_path = self.api.api("noref_wgs.snp_compare")
        if self.option("analysis_model") == "single":
            list_vcf = []
            vcf_path = self._sheet.output + "/snp/"
            path = os.listdir(self.output_dir+"/snp_compare/")
            for i in path:
                name = i.strip().split(".")[0]
                if re.match(".*xls$", i):
                    api_path.add_snp_compare_detail(self.output_dir + '/snp_compare/' + i,
                                                    self.option("main_id"), name)
                    self.logger.info("_______________-----------------------------------")
                    self.logger.info(i)

                elif re.match(".*results$", i):
                    api_path.add_sg_snp_compare_stat(self.output_dir + '/snp_compare/' + name + ".results", self.option("main_id"), name)
                elif re.match(".*vcf$", i):
                    list_vcf.append(vcf_path + i)
                else:
                    self.set_error("结果文件中存在着不明文件", code="15500154")
            api_path.update_snp_compare(list_vcf, self.option("main_id"))
        else:
            api_path.add_snp_compare_detail(self.output_dir+'/snp_compare/snp_compare.table.xls', self.option("main_id"), "snp_compare")
            api_path.add_sg_snp_compare_stat(self.output_dir+'/snp_compare/snp_compare.results', self.option("main_id"), "snp_compare")
            api_path.update_snp_compare(self.output_dir+'/snp_compare/snp_compare.vcf', self.option("main_id"))  # 此处需要修改

    def file_check(self):    #检查是否有results文件，在tool中，vcf为空的时候没有results结果，也就不需要导表了
        files = os.listdir(self.output_dir + "/snp_compare/")
        num = 1  # results不存在
        for j in files:
            if re.match(".*results$", j):
                num = 0
                break
        if num == 0:
            return True
        elif num ==1:
            return False

    def run(self):
        self.get_config()
        self.run_snp_compare()
        super(SnpCompareWorkflow, self).run()

    def end(self):
        if self.file_check():
            self.set_db()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SnpCompareWorkflow, self).end()
