# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180408

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class CircosModule(Module):
    def __init__(self, work_id):
        super(CircosModule, self).__init__(work_id)
        options = [
            {"name": "assemble", "type": "infile", "format": "sequence.fasta"},  # 序列文件
            {"name": "location", "type": "string", "default": "Scaffold"},  # 位置信息
            {"name": "trna", "type": "infile", "format": "gene_structure.gff3"},  # trna.gff文件
            {"name": "rrna", "type": "infile", "format": "gene_structure.gff3"},  # rrna.gff文件
            {"name": "gene", "type": "infile", "format": "gene_structure.gff3"},  # gene.gff文件
            {"name": "specimen_id", "type": "string", "default": ""},  # 样品名
            {"name": "anno_cog", "type": "infile", "format": "sequence.profile_table"},  # cog注释结果文件
            {"name": "labs", "type": "string", "default": "cog;ncrna"}, #分号分割
        ]
        self.add_option(options)
        self.pre_file = self.add_tool("bacgenome.circos_pre")
        self.plot_circos = self.add_tool("bacgenome.plot_circos")

    def check_options(self):
        if not self.option("assemble").is_set:
            raise OptionError("必须设置参数拼接序列文件", code="21401101")
        if not self.option("trna").is_set:
            raise OptionError("必须设置trna.gff文件", code="21401102")
        if not self.option("rrna").is_set:
            raise OptionError("必须设置rrna.gff文件", code="21401103")
        if not self.option("gene").is_set:
            raise OptionError("必须设置gene.gff文件", code="21401104")
        if not self.option("anno_cog").is_set:
            raise OptionError("必须设置cog注释结果文件", code="21401105")
        return True

    def run_pre(self):
        self.pre_file.set_options({
            "assemble": self.option("assemble"),
            "location": self.option("location"),
            "trna": self.option("trna"),
            "rrna": self.option("rrna"),
            "gene": self.option("gene"),
            "anno_cog": self.option("anno_cog")
        })
        self.pre_file.on('end', self.run_plot)
        self.pre_file.run()

    def run_plot(self):
        opts = {
            "k": self.pre_file.work_dir + "/karyotype.txt",
            "c": self.pre_file.work_dir + "/sense_strand_cog.txt",
            "t": self.pre_file.work_dir + "/temp.txt",
            "ac": self.pre_file.work_dir + "/antisense_strand_cog.txt",
            "pgc": self.pre_file.work_dir + "/positive_gc_count.txt",
            "ngc": self.pre_file.work_dir + "/negative_gc_count.txt",
            "pgs": self.pre_file.work_dir + "/positive_gc_skew.txt",
            "ngs": self.pre_file.work_dir + "/negative_gc_skew.txt",
            "labs" : self.option("labs")
        }
        if not self.option('location').startswith('p'):
            opts['f'] = self.pre_file.work_dir + "/ncRNA.txt"
            opts['type'] = 1
        else:
            opts['type'] = 2
        self.plot_circos.set_options(opts)
        self.plot_circos.on('end', self.set_output)
        self.plot_circos.run()

    def set_output(self):
        # if self.option("specimen_id").startswith('Chr'):
        #     name = self.option("specimen_id") + '_' + 'chromosome' + i[3:] + '_circos'
        # elif self.option("specimen_id").startswith('p'):
        #     name = self.option("specimen_id") + '_' + 'plasmid' + i[1:] + '_circos'
        # else:
        #     name = self.option("specimen_id") + '_whole_genome_circos'
        # if os.path.exists(self.output_dir + '/' + name + '.png'):
        #     os.remove(self.output_dir + '/' + name + '.png')
        # os.link(self.plot_circos.output_dir + '/circos.png', self.output_dir + '/' + name + '.png')
        # if os.path.exists(self.output_dir + '/' + name + '.svg'):
        #     os.remove(self.output_dir + '/' + name + '.svg')
        # os.link(self.plot_circos.output_dir + '/circos.svg', self.output_dir + '/' + name + '.svg')
        self.linkdir(self.plot_circos.output_dir, "")
        self.linkdir(self.pre_file.output_dir, "pre_file")
        self.end()

    def run(self):
        super(CircosModule, self).run()
        self.run_pre()

    def end(self):
        super(CircosModule, self).end()

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
