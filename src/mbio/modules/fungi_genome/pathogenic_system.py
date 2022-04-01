# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180316
# 增加dfvf 和 signalp


from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class PathogenicSystemModule(Module):
    def __init__(self, work_id):
        super(PathogenicSystemModule, self).__init__(work_id)
        options = [
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "analysis", "type": "string", "default": "all"}  # 挑选指定分析进行，默认进行所有
        ]
        self.add_option(options)
        self.anno_modules = []
        self.anno_dir={}
        self.tcdb = self.add_module("annotation.tcdb_dna_anno")
        self.card = self.add_module("annotation.card_dna_anno")
        #self.vfdb = self.add_module("annotation.vfdb_dna_anno")
        self.dfvf = self.add_module("annotation.dfvf_anno")
        self.phi = self.add_module("annotation.phi_dna_anno")
        self.tmhmm = self.add_tool("align.tmhmm_anno")
        self.signalp = self.add_tool("align.signalp_anno")

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="22101301")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="22101302")
        self.analysis_type = list(set(self.option("analysis").split(',')))
        if "all" in self.analysis_type:
        #    self.analysis_type = ["card", "vfdb", "tcdb", "phi", "tmhmm"]
            self.analysis_type = ["card", "dfvf", "tcdb", "phi", "tmhmm","signalp"]
        for one_analysis in self.analysis_type:
            if not one_analysis in ["card", "dfvf", "tcdb", "phi", "tmhmm","signalp"]:
                raise OptionError("分析项%s不被支持", variables=(one_analysis), code="22101303")
        return True

    def run_anno(self, run_tool):
        one_run = eval('self.' + run_tool)
        opt_dic = {
            "query": self.option("query"),
            "sample": self.option("sample")
        }
        if run_tool == "signalp":
            opt_dic["type"] = "euk"
            opt_dic["out_format"] = "short"
        if run_tool == "card":
            opt_dic["category"] = "drug_class" # 细菌升级添加 by ysh in 20190410
        one_run.set_options(opt_dic)
        self.anno_modules.append(one_run)
        self.anno_dir[run_tool] = one_run.output_dir

    def get_all_run(self):
        for i in self.analysis_type:
            self.run_anno(i)
        if len(self.anno_modules) > 1:
            self.on_rely(self.anno_modules, self.set_output)
            self.logger.info(self.anno_modules)
        else:
            self.anno_modules[0].on('end', self.set_output)
        for module in self.anno_modules:
            self.logger.info(module)
            module.run()

    def set_output(self):
        for i in self.analysis_type:
            self.logger.info(self.anno_dir[i])
            self.linkdir(self.anno_dir[i], i.upper())
        rm_list = ["gene_tcdb_anno.xls","tcdb_align_table.xls","top1_tcdb_align_table.xls"]
        for i in rm_list:
            tmp_path = self.output_dir +'/TCDB/'+ i
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        self.end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(PathogenicSystemModule, self).run()
        self.get_all_run()

    def end(self):
        super(PathogenicSystemModule, self).end()
