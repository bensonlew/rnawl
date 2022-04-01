# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0822

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import math


class TraitAnnovarAgent(Agent):
    """
    群体进化，GWAS关联分析“性状分布统计表”的tool
    """
    def __init__(self, parent):
        super(TraitAnnovarAgent, self).__init__(parent)
        options = [
            {"name": "upload_trait_path", "type": "infile", "format": "dna_gmap.trait"}  # 性状文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("upload_trait_path").is_set:
            raise OptionError("请设置性状文件upload_trait_path")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(TraitAnnovarAgent, self).end()


class TraitAnnovarTool(Tool):
    def __init__(self, config):
        super(TraitAnnovarTool, self).__init__(config)
        # self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.rscript = "program/R-3.3.3/bin/Rscript"
        self.trait_annovar = self.config.PACKAGE_DIR + "/dna_gmap/trait_analysis.R"

    def run_trait_annovar(self):
        """
        trait_annovar.R
        """
        cmd = "{} {} --infile {}".format(self.rscript, self.trait_annovar, self.option("upload_trait_path").prop["path"])
        cmd += " --outfile {}".format(os.path.join(self.work_dir, "annovar.trait.xls"))
        command = self.add_command("trait_annovar", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("trait_annovar.R运行完成")
        else:
            self.set_error("性状分析运行失败，请检查")
            raise Exception("trait_annovar.R运行失败，请检查")

    def get_trait_bar(self):
        """
        得到画柱状图的数据；组距100等分或者1000等分，目前100等分。
        """
        if not os.path.exists(self.work_dir + "/bar_dir"):
            os.mkdir(self.work_dir + "/bar_dir")
        with open(self.option("upload_trait_path").prop["path"], "r") as f:
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
                    step = float(max - min) / 100
                    for j in range(1, 100):
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
        anno_trit = self.option("upload_trait_path").prop["path"]
        f1 = os.path.join(self.output_dir, "annovar.trait.xls")
        f2 = os.path.join(self.output_dir, "annovar.sample.xls")
        if os.path.exists(f1):
            os.remove(f1)
        os.link(anno_sample, f1)
        if os.path.exists(f2):
            os.remove(f2)
        os.link(anno_trit, f2)

    def run(self):
        super(TraitAnnovarTool, self).run()
        self.run_trait_annovar()
        self.get_trait_bar()
        self.set_output()
        # self.set_db()
        self.end()
