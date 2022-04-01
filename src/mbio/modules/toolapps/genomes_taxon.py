# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.19

from biocluster.module import Module
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir


class GenomesTaxonModule(Module):
    """
    所有样品的物种注释分类
    """
    def __init__(self, work_id):
        super(GenomesTaxonModule, self).__init__(work_id)
        options = [
            {"name": "type", "type": "string", "default": "GTDB"},  # GTDB or custom
            {'name': '16s', 'type': 'infile', "format": "sequence.fasta_dir"},  ##基因组的16s序列
            {'name': 'house', 'type': 'infile', "format": "sequence.fasta_dir"},  ##基因组的看家基因序列
            {"name": "genomes", "type": "infile", "format": "sequence.fasta_dir"},  ##样品的基因组文件
            {"name": "ref_genomes", "type": "infile", "format": "sequence.fasta_dir"},  ##自定义数据库的基因组文件
            {"name": "ref_s16", "type": "outfile", "format": "sequence.fasta"},
            {'name': 'custom_fa', 'type': 'infile', "format": "sequence.fasta"},  ##custom时的数据库序列
            {"name": "cus_table", "type": "infile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.modules =[]

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("genomes").is_set:
            raise OptionError("必须设置参数genomes的目录！")

    def run_genomepredict(self):
        if self.option("type") in ['GTDB']:
            for file in os.listdir(self.option("16s").prop['path']):
                sample = file.split(".16S.ffn")[0]
                genome_ani = self.add_module('toolapps.genome_ani')
                opts = {
                    'type': self.option("type"),
                    'input': self.option("16s").prop['path'] +"/"+file,
                    'house': self.option("house").prop['path'] +"/"+sample + ".core_gene.fa",
                    'genome': self.option("genomes").prop['path'] +"/"+sample + ".fasta",
                    'sample': sample,
                }
                self.logger.info(opts)
                genome_ani.set_options(opts)
                self.modules.append(genome_ani)
        elif self.option("type") in ['custom']:
            for file in os.listdir(self.option("16s").prop['path']):
                sample = file.split(".16S.ffn")[0]
                genome_ani = self.add_module('toolapps.genome_ani')
                opts = {
                    'type': self.option("type"),
                    'input': self.option("16s").prop['path'] +"/"+file,
                    'house': self.option("house").prop['path'] +"/"+sample + ".core_gene.fa",
                    'genome': self.option("genomes").prop['path'] +"/"+sample + ".fasta",
                    'ref_s16': self.option("ref_s16"),
                    'ref_genomes': self.option("ref_genomes"),
                    'custom_fa': self.option("custom_fa"),
                    'cus_table': self.option("cus_table"),
                    'sample': sample,
                }
                genome_ani.set_options(opts)
                self.modules.append(genome_ani)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        elif len(self.modules) == 1:
            self.modules[0].on("end", self.set_output)
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(GenomesTaxonModule, self).run()
        self.run_genomepredict()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if os.path.exists(self.output_dir+"/anno"):
            shutil.rmtree(self.output_dir+"/anno")
        os.mkdir(self.output_dir+"/anno")
        if os.path.exists(self.output_dir+"/16s_blast"):
            shutil.rmtree(self.output_dir+"/16s_blast")
        os.mkdir(self.output_dir+"/16s_blast")
        if os.path.exists(self.output_dir+"/tree"):
            shutil.rmtree(self.output_dir+"/tree")
        os.mkdir(self.output_dir+"/tree")
        if len(self.modules) > 1:
            des = []
            for module in self.modules:
                if module.option("s16_table").is_set:
                    s16 = os.path.basename(module.option("s16_table").prop['path'])
                    os.link(module.option("s16_table").prop['path'], self.output_dir+"/16s_blast/"+s16)
                for file in os.listdir(module.output_dir):
                    if re.search(".taxon.xls", file):
                        des.append(module.output_dir+"/"+file)
                    elif re.search(".anno.xls", file):
                        os.link(module.output_dir+"/"+file, self.output_dir+"/anno/"+file)
                    elif re.search(".nwk", file):
                        os.link(module.output_dir+"/"+file, self.output_dir+"/tree/"+file)
            os.system("cat {} >{}".format(" ".join(des), self.output_dir + "/all.taxon.xls"))
        elif len(self.modules) == 1:
            if self.modules[0].option("s16_table").is_set:
                s16 = os.path.basename(self.modules[0].option("s16_table").prop['path'])
                os.link(self.modules[0].option("s16_table").prop['path'], self.output_dir + "/16s_blast/" + s16)
            for file in os.listdir(self.modules[0].output_dir):
                if re.search(".taxon.xls", file):
                    if os.path.exists(self.output_dir + "/all.taxon.xls"):
                        os.remove(self.output_dir + "/all.taxon.xls")
                    os.link(self.modules[0].output_dir + "/" + file, self.output_dir + "/all.taxon.xls")
                elif re.search(".anno.xls", file):
                    if os.path.exists(self.output_dir + "/anno/" + file):
                        os.remove(self.output_dir + "/anno/" + file)
                    os.link(self.modules[0].output_dir + "/" + file, self.output_dir + "/anno/" + file)
                elif re.search(".nwk", file):
                    if os.path.exists(self.output_dir + "/tree/" + file):
                        os.remove(self.output_dir + "/tree/" + file)
                    os.link(self.modules[0].output_dir + "/" + file, self.output_dir + "/tree/" + file)
        self.end()

    def end(self):
        super(GenomesTaxonModule, self).end()