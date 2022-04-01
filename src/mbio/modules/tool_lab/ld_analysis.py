# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modify:20180612
# last modified by binbinzhao@20190328

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import json


class LdAnalysisModule(Module):
    """
     连锁不平衡接口module。
     传入的分组情况以字典的形式。如下：
     {"a": "1,2,4", "b": "25,6,7", "c": "23,34,56"}
    """
    def __init__(self, work_id):
        super(LdAnalysisModule, self).__init__(work_id)
        options = [
            {"name": "vcf_file", "type": 'infile', "format": "dna_gmap.vcf"},
            {"name": "min_dp", "type": "string"},
            {"name": "max_dp", "type": "string"},
            {"name": "max_missing", "type": "string", "default": "0.3"},
            {"name": "min_maf", "type": "string", "default": "0.05"},
            {"name": "max_maf", "type": "string", "default": "1"},
            {"name": "group_dict", "type": "string"},
            {"name": "recode", "type": "bool"},
            {"name": "task_id", "type": "string"},
            {"name": "project_sn", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            # {"name": "graph_path", "type": "string"}
        ]
        self.add_option(options)
        self.vcftools_filter_tools = []
        self.ld_decay_tools = []
        self.group_sample = {}

    def check_options(self):
        a = vars(self.option("vcf_file"))
        if not self.option("vcf_file"):
            raise OptionError("请输入vcf.file文件", code="24900601")
        if not self.option("max_missing"):
            raise OptionError("请输入缺失率", code="24900602")
        else:
            if float(self.option("max_missing")) < 0 or float(self.option("max_missing")) > 1:
                raise OptionError("缺失率应该为0到1之间的数字", code="24900603")
        if not self.option("min_maf"):
            raise OptionError("请输入次要等位基因最小值", code="24900604")
        else:
            if float(self.option("min_maf")) > 1:
                raise OptionError("次要等位基因频率最大值为1", code="24900605")
        if not self.option("max_maf"):
            raise OptionError("请输入次要等位基因最大值", code="24900606")
        else:
            if float(self.option("max_maf")) > 1:
                raise OptionError("次要等位基因频率最大值为1", code="24900607")
        # self.logger.info("___________________________")
        # self.logger.info(self.option("min_dp"))
        # print "___________________________"
        # print self.option("min_dp")
        # if self.option("min_dp"):
        #     if not float(self.option("min_dp")) > 0:
        #         raise OptionError("平均测序深度必须为正整数")
        # if self.option("max_dp"):
        #     if not float(self.option("max_dp")) > 0:
        #         raise OptionError("平均测序深度必须为正整数")
        return True

    def run_vcftools_filter(self):
        # group_list_ = {"group_list1": "Tibetan_1,Tibetan_2,Tibetan_3,Tibetan_4,Tibetan_5",
        #                              "group_list2": "XSBN_1,XSBN_2,XSBN_3,XSBN_4,XSBN_5",
        #                              "group_list3": "RED_1,RED_2,RED_3,RED_4,RED_5"}
        group_list_ = json.loads(self.option("group_dict"))
        if not os.path.exists(self.work_dir + "/group_list"):
            os.mkdir(self.work_dir + "/group_list")
        else:
            os.system("rm -r {}".format(self.work_dir + "/group_list"))
            os.mkdir(self.work_dir + "/group_list")
        for group_name in group_list_.keys():  # 这里的self.output_dir指的是tool的dir还是module的，应该是module的。
            path = os.path.join(self.work_dir, "group_list")
            if os.path.exists(os.path.join(path,  group_name)):
                os.remove(os.path.join(path,  group_name))
            with open(os.path.join(path,  group_name), "w") as w:
                list = group_list_[group_name].strip().split(",")
                # list = group_list_[group_name]
                self.group_sample[str(group_name)] = len(list)
                for i in list:
                    w.write(i + "\n")
        self.logger.info("self.group_sample：{}".format(self.group_sample))
        for g in os.listdir(os.path.join(self.work_dir, "group_list")):
            keep = os.path.join(os.path.join(self.work_dir, "group_list"), g)
            vcftools_filter = self.add_tool("dna_evolution.vcftools_filter")
            options = ({
                "vcf_path": self.option("vcf_file").prop["path"],
                "recode": True,  # 未设置传参，直接写死True
                "max_missing": float(self.option('max_missing')),
                "keep": keep,
                "min_maf": self.option("min_maf"),
                "group_name": g,
                "max_maf": float(self.option('max_maf'))
            })
            if self.option('min_dp'):
                options["minDP"] = int(self.option('min_dp'))
            if self.option('max_dp'):
                options["maxDP"] = int(self.option('max_dp'))
            vcftools_filter.set_options(options)
            self.vcftools_filter_tools.append(vcftools_filter)
        vcftools_filter1 = self.add_tool("dna_evolution.vcftools_filter")  # 需要多跑一个“all”的tool。
        options = ({
            "vcf_path": self.option("vcf_file").prop["path"],
            "recode": True,  # 未设置传参，直接写死True
            "max_missing": float(self.option('max_missing')),
            "min_maf": self.option("min_maf"),
            "group_name": "Overall",
            "max_maf": float(self.option('max_maf'))
        })
        vcftools_filter1.set_options(options)
        if self.option('min_dp'):
            options["minDP"] = int(self.option('min_dp'))
        if self.option('max_dp'):
            options["maxDP"] = int(self.option('max_dp'))
        self.vcftools_filter_tools.append(vcftools_filter1)
        for j in range(len(self.vcftools_filter_tools)):
            self.vcftools_filter_tools[j].on("end", self.set_output, 'vcf_filter_dir')
        if self.vcftools_filter_tools:
            if len(self.vcftools_filter_tools) > 1:
                self.on_rely(self.vcftools_filter_tools, self.ld_decay_run)
            elif len(self.vcftools_filter_tools) == 1:
                self.vcftools_filter_tools[0].on('end', self.ld_decay_run)
        else:
            raise Exception("vcftools_filter_tools列表为空！")
        for tool in self.vcftools_filter_tools:
            gevent.sleep(1)
            tool.run()

    def ld_decay_run(self):
        if len(self.vcftools_filter_tools) > 1:
            result_path = self.output_dir + "/vcf_filter_dir"
        else:
            result_path = self.vcftools_filter_tools[0].output_dir
        for i in os.listdir(result_path):
            if not self.check_vcf_none(os.path.join(result_path, i)):  # 不为空以及不是一个样本
                ld_decay = self.add_tool("dna_evolution.ld_decay")
                options = {
                    "vcf_file": os.path.join(result_path, i),
                    "maf": float(self.option("min_maf")),
                    "miss": float(self.option("max_missing")),
                    "group_name": i.split(".")[0]
                }
                ld_decay.set_options(options)
                self.ld_decay_tools.append(ld_decay)
            else:
                continue
        for j in range(len(self.ld_decay_tools)):
            self.ld_decay_tools[j].on("end", self.set_output, 'ld_decay_dir')
        if self.ld_decay_tools:
            if len(self.ld_decay_tools) > 1:
                self.on_rely(self.ld_decay_tools, self.merge_file)
            elif len(self.ld_decay_tools) == 1:
                self.ld_decay_tools[0].on('end', self.merge_file)
        else:
            raise Exception("ld_decay_tools列表为空！")
        for tool in self.ld_decay_tools:
            gevent.sleep(1)
            tool.run()
    #
    # def check_end(self, file_dir):
    #     """
    #     当有一个分组的结果为空的时候，流程就截止
    #     :return:
    #     """
    #     for m in os.listdir(file_dir):
    #         self.logger.info("3333333:{}".format(os.path.join(file_dir, m)))
    #         if self.check_vcf_none(os.path.join(file_dir, m)):
    #             return True
    #     return False

    def check_vcf_none(self, vcf_path):
        """
        默认vcf不是空的以及是一个样本的分组
        :param vcf_path:
        :return:
        """
        with open(vcf_path, "r") as f:
            num = 0
            for lines in f:
                if re.match('#CHROM', lines):  # 检查分组中样本不是一个
                    if len(lines.strip().split('\t')) < 11:
                        num = 0
                        break
                if not re.match('#', lines):  # 检查vcf不为空
                    num = 1
                    break
            if num == 0:
                return True
            else:
                return False

    def merge_file(self):  # 用于ld_decay_run生成的文件写入一个新的文件内。
        path1 = os.path.join(self.output_dir, "ld_decay_dir")
        if os.listdir(path1) == "":
            raise Exception("缺少.stat.gz文件")
        else:
            if os.path.exists(os.path.join(self.work_dir, "gro_list")):
                os.remove(os.path.join(self.work_dir, "gro_list"))
            with open(os.path.join(self.work_dir, "gro_list"), "w") as w:
                w.write("popid" + "\t" + "file" + "\tsampleNum\n")
                for i in os.listdir(path1):
                    self.logger.info("group_name:{}".format(i.split(".")[0]))
                    if i.split(".")[0] in ['all', 'Overall']:
                        sample_num = 0
                        for key in self.group_sample.keys():
                            sample_num += self.group_sample[key]
                    else:
                        try:
                            sample_num = self.group_sample[i.split(".")[0]]
                        except:
                            raise Exception("获取分组中样本数目失败！")
                    w.write(i.split(".")[0] + "\t" +
                            os.path.join(os.path.join(self.output_dir, "ld_decay_dir"), i) + "\t" + str(sample_num) + "\n")
        self.ld_draw = self.add_tool("dna_evolution.ld_draw")
        options = {
            "gro_list": os.path.join(self.work_dir, "gro_list")
        }
        self.ld_draw.set_options(options)
        self.ld_draw.on("end", self.set_output, 'ld_draw_dir')
        gevent.sleep(1)
        self.ld_draw.run()

    # def new_group_list(self):
    #     path2 = os.path.join(self.output_dir, "ld_graph_dir")
    #     if os.listdir(path2) == "":
    #         raise Exception("缺少绘图文件，请核实")
    #     else:
    #         if os.path.exists(os.path.join(self.work_dir, "new_gro_list")):
    #             os.remove(os.path.join(self.work_dir, "new_gro_list"))
    #         with open(os.path.join(self.work_dir, "new_gro_list"), "w") as w:
    #             w.write("popid" + "\t" + "file" + "\n")
    #             for i in os.listdir(path2):
    #                 w.write(i.split(".")[0] + "\t" + os.path.join(os.path.join(self.output_dir, "ld_graph_dir"), i) +
    #                         "\n")
    #     if self.option("main_id"):
    #         self.set_db()
    #     self.end()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'vcf_filter_dir':
            self.linkdir(obj.output_dir, 'vcf_filter_dir')
        elif event['data'] == 'ld_decay_dir':
            self.linkdir(obj.output_dir, 'ld_decay_dir')
        elif event['data'] == 'ld_draw_dir':
            self.linkdir(obj.output_dir, 'ld_draw_dir')
            # self.new_group_list()
            self.end()

        else:
            pass

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

        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    # def set_db(self):
    #     """
    #     保存结果到mongo
    #     """
    #     if self.option("main_id"):
    #         # graph_path = self.option("graph_path") + "ld_analysis/ld_draw_dir/ld.png"
    #         graph_path = self._sheet.output + "/ld_draw_dir/ld.png"
    #         api_ld_analysis = self.api.api("dna_evolution.ld_analysis")
    #         api_ld_analysis.sg_ld_detail(self.option("main_id"), self.work_dir + "/gro_list")
    #         # api_ld_analysis.add_ld_curve(self.work_dir + "/new_gro_list", self.option("main_id"),
    #         #  self.option("task_id"))
    #         api_ld_analysis.update_ld_analysis(graph_path, self.option("main_id"))
    #     self.end()

    def run(self):
        super(LdAnalysisModule, self).run()
        self.run_vcftools_filter()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(LdAnalysisModule, self).end()
