# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.11

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import math


class TraitAnnovarAgent(Agent):
    """
    遗传图谱：性状分析
    """
    def __init__(self, parent):
        super(TraitAnnovarAgent, self).__init__(parent)
        options = [
            {"name": "trait_file", "type": "infile", "format": "dna_gmap.trait"},  # 性状文件
            {"name": "sample_list", "type": "infile", "format": "dna_gmap.path"},  # 样本list
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("trait_file").is_set:
            raise OptionError("请设置性状文件trait_file", code="34801901")
        if not self.option("sample_list").is_set:
            raise OptionError("请设置参与性状分析的样本list:sample_list", code="34801902")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(TraitAnnovarAgent, self).end()


class TraitAnnovarTool(Tool):
    def __init__(self, config):
        super(TraitAnnovarTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.rscript = "program/R-3.3.3/bin/Rscript"
        self.trait_annovar = self.config.PACKAGE_DIR + "/dna_gmap/trit-annovar.pl"
        self.trait_annovar_new = self.config.PACKAGE_DIR + "/dna_gmap/trait_analysis.R"
        # if self.option("trait_file").prop["is_excel"]:
        #     self.new_trait = self.option("trait_file").prop["path"]
        # else:
        #     self.new_trait = os.path.join(self.work_dir, "new_trait.xls")
        #     self.option("trait_file").get_xls_trait(self.new_trait)

    def run_trait_annovar(self):
        """
        trait-annovar.pl
        """
        cmd = "{} {} -trit {}".format(self.perl_path, self.trait_annovar, self.new_trait)
        cmd += " -select {} -out {}".format(self.option("sample_list").prop["path"], self.work_dir)
        command = self.add_command("trait_annovar", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("trait-annovar.pl运行完成")
        else:
            self.set_error("性状分析运行失败，请检查", code="34801901")
            self.set_error("trait-annovar.pl运行失败，请检查", code="34801907")

    def run_trait_annovar_new(self):
        """
        trait_annovar.R
        """
        cmd = "{} {} --infile {}".format(self.rscript, self.trait_annovar_new, self.option("trait_file").prop["path"])
        cmd += " --outfile {}".format(os.path.join(self.work_dir, "annovar.trait.xls"))
        command = self.add_command("new_trait_annovar", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("trait_annovar.R运行完成")
        else:
            self.set_error("性状分析运行失败，请检查", code="34801902")
            self.set_error("trait_annovar.R运行失败，请检查", code="34801908")

    def get_trait_bar(self):
        """
        得到画柱状图的数据
        """
        if not os.path.exists(self.work_dir + "/bar_dir"):
            os.mkdir(self.work_dir + "/bar_dir")
        with open(self.option("trait_file").prop["path"], "r") as f:
            lines = f.readlines()
            head = lines[0].strip().split("\t")
            trit_dict = {}
            for i in range(1, len(head)):
                trit_dict[i] = []
            for line in lines[1:]:
                item = line.strip().split("\t")
                for i in range(1, len(item)):
                    try:
                        if item[i] != "NaN":
                            trit_dict[i].append(float(item[i]))
                    except:
                        self.logger.info("样本{}的性状{}是空的{}".format(item[0], head[i], item[i]))
            for i in range(1, len(head)):
                with open(self.work_dir + "/bar_dir/" + head[i] + ".bar.txt", "w") as w:
                    trit_dict[i].sort()
                    min = math.floor(trit_dict[i][0])
                    max = math.ceil(trit_dict[i][-1])
                    step = float(max - min) / 20
                    for j in range(1, 20):
                        min_ = min + step * (j - 1)
                        max_ = min + step * j
                        sample_num = 0
                        for v in trit_dict[i]:
                            if min <= v < max_:
                                sample_num += 1
                                trit_dict[i].remove(v)
                            else:
                                break
                        w.write(str(min_) + "-" + str(max_) + "\t" + str(sample_num) + "\n")
                    w.write(str(max_) + "-" + str(max) + "\t" + str(len(trit_dict[i])) + "\n")

    def set_output(self):
        anno_sample = os.path.join(self.work_dir, "annovar.trait.xls")
        # anno_trit = os.path.join(self.work_dir, "annovar.trit.xls")
        anno_trit = self.option("trait_file").prop["path"]
        f1 = os.path.join(self.output_dir, "annovar.trait.xls")
        f2 = os.path.join(self.output_dir, "annovar.sample.xls")
        if os.path.exists(f1):
            os.remove(f1)
        os.link(anno_sample, f1)
        if os.path.exists(f2):
            os.remove(f2)
        os.link(anno_trit, f2)

    def set_db(self):
        """
        将结果导入mongo数据库
        """
        self.logger.info("将结果导入mongo数据库")
        trait_api = self.api.api("dna_gmap.trait_annovar")
        feature_id = self.option("main_id")
        task_id = self.option("task_id")
        anno_trait = os.path.join(self.output_dir, "annovar.trait.xls")
        anno_sample = os.path.join(self.output_dir, "annovar.sample.xls")
        bar_dir = os.path.join(self.work_dir, "bar_dir")
        trait_api.add_sg_feature_annovar(task_id, feature_id, anno_trait)
        trait_api.add_sg_feature_samples(feature_id, anno_sample)
        trait_api.add_sg_feature_bar(feature_id, task_id, bar_dir)

    def run(self):
        super(TraitAnnovarTool, self).run()
        # self.run_trait_annovar()
        self.run_trait_annovar_new()
        self.get_trait_bar()
        self.set_output()
        if self.option("main_id"):
            self.set_db()
        self.end()
