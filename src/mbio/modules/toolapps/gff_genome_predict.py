# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.14

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir


class GffGenomePredictModule(Module):
    """
    所有基因组的结果统计
    """
    def __init__(self, work_id):
        super(GffGenomePredictModule, self).__init__(work_id)
        options = [
            {"name": "gffs", "type": "infile", "format": "toolapps.gff_dir"},
            {"name": "genomes", "type": "infile", "format": "sequence.fasta_dir"},
            {'name': 'map_info', 'type': 'infile', "format": "sequence.profile_table"},
            {'name': 'stat', 'type': 'outfile', "format": "sequence.profile_table"},
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

    def get_info(self):
        dict ={}
        with open(self.option("map_info").prop['path']) as f:
            lines = f.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                dict[lin[0]]=[lin[1],lin[2]]
        return dict

    def run_genomepredict(self):
        self.gff_dict = self.get_info()
        for sample in self.gff_dict.keys():
            genome_pre = self.add_tool('toolapps.gff_genome_predict')
            opts = {
                'input_genome': self.option("genomes").prop['path'] + "/" + self.gff_dict[sample][0],
                'gff': self.option("gffs").prop['path'] + "/" + self.gff_dict[sample][1],
                'sample': sample,
            }
            genome_pre.set_options(opts)
            self.modules.append(genome_pre)
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
        super(GffGenomePredictModule, self).run()
        self.run_genomepredict()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if os.path.exists(self.output_dir+"/house"):
            shutil.rmtree(self.output_dir+"/house")
        os.mkdir(self.output_dir+"/house")
        if os.path.exists(self.output_dir+"/s16"):
            shutil.rmtree(self.output_dir+"/s16")
        os.mkdir(self.output_dir+"/s16")
        if len(self.modules) > 1:
            des = []
            for module in self.modules:
                if module.option("s16").is_set:
                    file = os.path.basename(module.option("s16").prop['path'])
                    os.link(module.option("s16").prop['path'], self.output_dir + "/s16/" + file)
                file2 = os.path.basename(module.option("house").prop['path'])
                os.link(module.option("house").prop['path'], self.output_dir + "/house/" + file2)
                des.append(module.option("stat").prop['path'])
            if os.path.exists(self.output_dir + "/all.stat.xls"):
                os.remove(self.output_dir + "/all.stat.xls")
            os.system("cat {} >{}".format(" ".join(des), self.output_dir + "/all.stat.xls"))
            self.option("stat", self.output_dir + "/all.stat.xls")
        elif len(self.modules) == 1:
            if self.modules[0].option("s16").is_set:
                file = os.path.basename(self.modules[0].option("s16").prop['path'])
                os.link(self.modules[0].option("s16").prop['path'], self.output_dir + "/s16/" + file)
            file2 = os.path.basename(self.modules[0].option("house").prop['path'])
            os.link(self.modules[0].option("house").prop['path'], self.output_dir + "/house/" + file2)
            if os.path.exists(self.output_dir + "/all.stat.xls"):
                os.remove(self.output_dir + "/all.stat.xls")
            os.link(self.modules[0].option("stat").prop['path'], self.output_dir + "/all.stat.xls")
            self.option("stat", self.output_dir + "/all.stat.xls")
        self.end()

    def end(self):
        super(GffGenomePredictModule, self).end()