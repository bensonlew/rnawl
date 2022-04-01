# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180416

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class CombineCircosModule(Module):
    def __init__(self, work_id):
        super(CombineCircosModule, self).__init__(work_id)
        options = [
            {"name": "assemble", "type": "infile", "format": "sequence.fasta"},  # 拼接序列文件
            {"name": "location", "type": "string", "default": "Scaffold"},  # 位置信息
            {"name": "trna", "type": "infile", "format": "gene_structure.gff3"},  # trna.gff文件
            {"name": "rrna", "type": "infile", "format": "gene_structure.gff3"},  # rrna.gff文件
            {"name": "gene", "type": "infile", "format": "gene_structure.gff3"},  # gene.gff文件
            {"name": "specimen_id", "type": "string", "default": ""},  # 样品名
            {"name": "anno_cog", "type": "infile", "format": "sequence.profile_table"}  # cog注释结果文件
        ]
        self.add_option(options)
        self.module_list = []

    def check_options(self):
        if not self.option("assemble").is_set:
            raise OptionError("必须设置参数拼接序列文件", code="21401301")
        if not self.option("trna").is_set:
            raise OptionError("必须设置trna.gff文件", code="21401302")
        if not self.option("rrna").is_set:
            raise OptionError("必须设置rrna.gff文件", code="21401303")
        if not self.option("gene").is_set:
            raise OptionError("必须设置gene.gff文件", code="21401304")
        if not self.option("anno_cog").is_set:
            raise OptionError("必须设置cog注释结果文件", code="21401305")
        return True

    def run_circos(self):
        self.location = []
        self.location1 = self.option("location").split(",")
        for one in self.location1:
            if one.startswith('Chromosome'):
                self.location.append(one)
            elif one.startswith('chromosome'):
                self.location.append(one)
            elif one.startswith('Plasmid'):
                if os.path.exists(self.work_dir + '/' + one + '.gff'):
                    self.location.append(one)
        self.circos_out = {}
        if self.option("location") != "Scaffold":
            for one in self.location:
                if not os.path.exists(self.work_dir + '/' + one + '.fa'):
                    self.set_error("不存在%s的拼接文件", variables=(one), code="21401301")
                self.circos = self.add_module("bacgenome.circos")
                self.circos.set_options({
                    "assemble": self.work_dir + '/' + one + '.fa',
                    "location": one,
                    "trna": self.option("trna"),
                    "rrna": self.option("rrna"),
                    "gene": self.work_dir + '/' + one + '.gff',
                    "anno_cog": self.option("anno_cog"),
                    "specimen_id": self.option("specimen_id"),
                })
                self.circos_out[one] = self.circos.output_dir
                self.module_list.append(self.circos)
        else:
            self.circos = self.add_module("bacgenome.circos")
            self.circos.set_options({
                "assemble": self.option("assemble"),
                "location": "Scaffold",
                "trna": self.option("trna"),
                "rrna": self.option("rrna"),
                "gene": self.option("gene"),
                "anno_cog": self.option("anno_cog"),
                "specimen_id": self.option("specimen_id"),
            })
            self.circos_out["Scaffold"] = self.circos.output_dir
            self.module_list.append(self.circos)
        if len(self.module_list) > 1:
            self.on_rely(self.module_list, self.set_output)
            self.logger.info(self.module_list)
        else:
            self.module_list[0].on('end', self.set_output)
        for module in self.module_list:
            module.run()

    def run_split(self):
        if self.option("location") != "Scaffold":
            self.option("assemble").split_single_seq(self.work_dir)
            self.option("gene").split_gff_by_name(self.work_dir)
        self.run_circos()

    def set_output(self):
        if self.option("location") != "Scaffold":
            for one in self.location:
                if os.path.exists(self.output_dir + '/' + one):
                    shutil.rmtree(self.output_dir + '/' + one)
                shutil.copytree(self.circos_out[one], self.output_dir + '/' + one)
        else:
            if os.path.exists(self.output_dir + "/Scaffold"):
                shutil.rmtree(self.output_dir + "/Scaffold")
            shutil.copytree(self.circos_out["Scaffold"], self.output_dir + "/Scaffold")
        self.end()

    def run(self):
        super(CombineCircosModule, self).run()
        self.run_split()

    def end(self):
        super(CombineCircosModule, self).end()

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
