# -*- coding:utf-8 -*-
# __author__ = 'chenyanyan'
# last_modifiy:2016.09.26

from __future__ import division
from biocluster.agent import Agent
from biocluster.tool import Tool
import shutil
import os
import glob
from biocluster.core.exceptions import OptionError
import json
import time


class StarIndexAgent(Agent):
    """
    star 比对工具
    """

    def __init__(self, parent):
        super(StarIndexAgent, self).__init__(parent)
        options = [

            {"name": "ref_genome_custom", "type": "infile", "format": "sequence.fasta"},  # 用户上传参考基因组文件
            {"name": "ref_genome", "type": "string"},  # 参考基因组模式选项 用户自定义、选择已有生物物种
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因组的gtf文件 ，gtf文件和fasta文件要配套使用
            {"name": "seq_method", "type": "string"},
            {"name": "is_indexed", "type": "bool", "default": False},
            {"name": "star_index1", "type": "outfile", "format": "align.star.star_index"}

        ]
        self.add_option(options)
        self.step.add_steps('star')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.star.start()
        self.step.update()

    def step_end(self):
        self.step.star.finish()
        self.step.update()

    def check_options(self):

        """
        检查参数设置
        """
        if self.option("ref_genome") == "customer_mode" and not self.option("ref_genome_custom").is_set:
            raise OptionError("请上传自定义参考基因组")

    def set_resource(self):
        self._cpu = 4
        if self.option("ref_genome") == "customer_mode":
            if self.option("ref_genome_custom").prop["size"] / 1024 / 1024 < 512:
                self._memory = '20G'
            elif self.option("ref_genome_custom").prop["size"] / 1024 / 1024 < 1024:
                self._memory = '40G'
            else:
                self._memory = '100G'
        else:
            self._memory = '100G'   # 设置资源大小

    def end(self):
        super(StarIndexAgent, self).end()    # 继承超类的end方法


class StarIndexTool(Tool):

    def __init__(self, config):
        super(StarIndexTool, self).__init__(config)
        self.star_path = "bioinfo/rna/star-2.5/bin/Linux_x86_64/"  # 设置star的路径
        self.samtools_path = self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.3.1/'

    def star_index1(self, genomeDir, ref_fa):
        """
        step1:第一步建索引；用star建立参考基因组的索引，当用户不上传参考基因组时，该步骤省略，直接调用已有的序列文件
        """
        if self.option("ref_genome") == "customer_mode":
            if self.option("ref_genome_custom").prop["size"] / 1024 / 1024 < 512:
                ram = str(20 * 1024 * 1024 * 1024)
            elif self.option("ref_genome_custom").prop["size"] / 1024 / 1024 < 1024:
                ram = str(40 * 1024 * 1024 * 1024)
            else:
                ram = str(100 * 1024 * 1024 * 1024)
        else:
            ram = str(100 * 1024 * 1024 * 1024)
        cmd = "{}STAR --runMode genomeGenerate --limitGenomeGenerateRAM {} --genomeDir {} --genomeFastaFiles {} --runThreadN 10".format(self.star_path, ram, genomeDir, ref_fa)
        self.logger.info("使用star建立参考基因组索引")
        command = self.add_command("star_index1", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("成功构建参考序列索引index1！")
        else:
            ram = str(100 * 1024 * 1024 * 1024)
            cmd = "{}STAR --runMode genomeGenerate --limitGenomeGenerateRAM {} --genomeDir {} --genomeFastaFiles {} --runThreadN 10".format(self.star_path, ram, genomeDir, ref_fa)
            # command.rerun()
            command = self.add_command("star_index1", cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("成功构建参考序列索引index1！")
            else:
                self.set_error("构建索引出错!")
                raise Exception("运行star出错")

    def run(self):
        """
        运行
        """
        super(StarIndexTool, self).run()

        if self.option("ref_genome") == "customer_mode" and self.option("ref_genome_custom").is_set:  # 自定义模式的时候需要建索引
            self.logger.info("在参考基因组为自定义模式下运行star！")
            ref_fa = self.option("ref_genome_custom").prop["path"]
            if not os.path.exists("ref_star_index1"):
                os.mkdir("ref_star_index1")
            genomeDir_path1 = os.path.join(self.work_dir, "ref_star_index1")   # 准备第一次建索引的路径
            if self.option("is_indexed") is False:
                self.star_index1(genomeDir_path1, ref_fa)  # 第一步：建索引，传入第一步索引的文件夹（此时是空文件夹）
            self.option("star_index1", genomeDir_path1)
        else:  # 参考基因组来自数据库
            ref_genome_json = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.json"
            with open(ref_genome_json, "r") as f:
                ref_dict = json.loads(f.read())
            rel_index = ref_dict[self.option("ref_genome")]["dna_index"]
            abs_index = self.config.SOFTWARE_DIR +  "/database/Genome_DB_finish/" +  rel_index
            self.logger.info(abs_index)
            self.option("star_index1", abs_index)
        self.end()
