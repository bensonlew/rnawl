#-*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# 20190930

import os
import re
import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_anno import CommonAnno
from mbio.packages.bac_comp_genome.common_function import link_file



class BlastclustAgent(Agent):
    """
    可以用于长序列聚类
    """

    def __init__(self, parent):
        super(BlastclustAgent, self).__init__(parent)
        options = [
            {"name": "dir", "type": "string", "format": "sequence.fasta"},  # 序列的文件夹路径
            {"name": "indentity", "type": "string", "default": "50"},
            {"name": "coverage", "type": "string", "default": "0.5"},
            {"name": "method", "type": "string", "default": "prephage"},#prephage,island
            {"name": "sample_list", "type": "infile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("dir"):
            raise OptionError("必须设置参数genome_name")
        if not self.option("indentity"):
            raise OptionError("必须设置参数indentity")
        if not self.option("coverage"):
            raise OptionError("必须设置参数coverage")
        if not self.option("method"):
            raise OptionError("必须设置参数method")
        else:
            if self.option("method") not in ["prephage", "island"]:
                raise OptionError("请提供正确method！")

    def set_resource(self):
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(BlastclustAgent, self).end()


class BlastclustTool(Tool):
    def __init__(self, config):
        super(BlastclustTool, self).__init__(config)
        self.blastclust = "bioinfo/metaGenomic/blast-2.2.18/bin/blastclust"

    def run_blastclust(self):
        """
        根据序列用blastclust进行序列聚类
        :return:
        """
        list = []
        cmd = ''
        for i in os.listdir(self.option("dir")):
            if re.search('.fna',i):
                list.append(self.option("dir") + "/" + i)
        des=' '.join(list)
        if self.option("method") in ['prephage']:
            os.system("cat {} > {}".format(des, self.work_dir + "/all.prephage.fna"))
            cmd = "{} -i {} -o {} -L {} -S {} -e F -a 5 -p F".format(self.blastclust,self.work_dir + "/all.prephage.fna",self.work_dir + "/all.blastclust.xls",
            self.option("coverage"),self.option("indentity"))
        elif self.option("method") in ['island']:
            os.system("cat {} > {}".format(des, self.work_dir + "/all.island.fna"))
            cmd = "{} -i {} -o {} -L {} -S {} -e F -a 5 -p F".format(self.blastclust, self.work_dir + "/all.island.fna",
                                                                self.work_dir + "/all.blastclust.xls",
                                                                self.option("coverage"), self.option("indentity"))
        self.logger.info(cmd)
        command = self.add_command("run_blastclust", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.stat_summary()
            self.logger.info("run_blastclust文件运行完成")
        else:
            self.set_error("run_blastclust文件运行出错!")

    def stat_blastclust(self, file, out):
        with open (file, "r") as f,open (out, "w") as g:
            if self.option("method") in ['prephage']:
                g.write("Type\tsample\tPh_id\n")
            elif self.option("method") in ['island']:
                g.write("Type\tsample\tGI_id\n")
            lines = f.readlines()
            n=0
            for line in lines:
                n += 1
                lin = line.strip().split(" ")
                for i in lin:
                    if self.option("method") in ['prephage']:
                        ph, sample=i.split("__")
                        g.write("Ph_Cluster"+ str(n) + "\t" + sample + "\t" + ph + "\n")
                    elif self.option("method") in ['island']:
                        island, sample= i.split("__")
                        g.write("GI_Cluster" + str(n) + "\t" + sample + "\t" + island + "\n")

    def stat_summary(self):
        self.sample_list = self.get_sample(self.option("sample_list").prop['path'])
        if self.option("method") in ['prephage']:
            self.stat_blastclust(self.work_dir + "/all.blastclust.xls", self.work_dir + "/all.prephage_type.xls")
            anno = CommonAnno()
            anno.anno_abund(self.work_dir + "/all.prephage_type.xls", "Type", "Ph_id",
                            self.work_dir + "/all.prephage_abund.txt")
            anno.fill_data(self.work_dir + "/all.prephage_abund.txt",self.sample_list, self.work_dir + "/all.prephage_abund.xls", "0", "Type")
            anno.anno_ano_genelist(self.work_dir + "/all.prephage_type.xls", "Type", "Ph_id",
                                   self.work_dir + "/all.prephage_genelist.txt")
            anno.fill_data(self.work_dir + "/all.prephage_genelist.txt", self.sample_list,
                           self.work_dir + "/all.prephage_genelist.xls", "-", "Type")
        elif self.option("method") in ['island']:
            self.stat_blastclust(self.work_dir + "/all.blastclust.xls", self.work_dir + "/all.island_type.xls")
            anno = CommonAnno()
            anno.anno_abund(self.work_dir + "/all.island_type.xls", "Type", "GI_id",
                            self.work_dir + "/all.island_abund.txt")
            anno.fill_data(self.work_dir + "/all.island_abund.txt", self.sample_list,
                           self.work_dir + "/all.island_abund.xls", "0", "Type")
            anno.anno_ano_genelist(self.work_dir + "/all.island_type.xls", "Type", "GI_id",
                                   self.work_dir + "/all.island_genelist.txt")
            anno.fill_data(self.work_dir + "/all.island_genelist.txt", self.sample_list,
                           self.work_dir + "/all.island_genelist.xls", "-", "Type")


    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("正在设置结果文件目录")
        if self.option("method") in ['prephage']:
            for i in ["all.prephage_abund.xls", "all.prephage_genelist.xls"]:
                link_file(self.work_dir + "/" + i, self.output_dir + "/" + i)
        elif self.option("method") in ['island']:
            for i in ["all.island_abund.xls", "all.island_genelist.xls"]:
                link_file(self.work_dir + "/" + i, self.output_dir + "/" + i)

    def get_sample(self,file):
        list = []
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                list.append(lin[0])
        return list

    def run(self):
        super(BlastclustTool, self).run()
        self.run_blastclust()
        self.set_output()
        self.end()