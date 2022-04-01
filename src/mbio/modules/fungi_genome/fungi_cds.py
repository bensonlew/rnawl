# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modify:2018.5.26

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError
import datetime
from mbio.packages.fungi_genome.get_genemark_key import update_gm_key_pip

class FungiCdsModule(Module):
    def __init__(self, work_id):
        super(FungiCdsModule, self).__init__(work_id)
        options = [
            {"name": "masked", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "cegma_gff", "type":"string", "default": "" },
            {"name": "sample", "type": "string","default":""}, # 样品名称
            {"name": "ref_protein","type": "infile", "format":"sequence.fasta"} #参考蛋白

        ]
        self.add_option(options)

        self.genemark = self.add_tool("predict.genemark_es")
        self.snap1 = self.add_tool("predict.gff_to_pyuhmm")
        self.snap2 = self.add_tool("predict.gff_to_pyuhmm")
        self.maker1 = self.add_tool("predict.maker")
        self.augustus = self.add_tool("predict.augustus")
        self.maker2 = self.add_tool("predict.maker")
        self.format_stat = self.add_tool("fungi_genome.gff_faa_fnn_format")
        self.rm_include = self.add_tool("fungi_genome.gene_rm_include")
        self.augustus_time = ''
        #self.tool1 = []
        #self.genemark_gmhmm_mod = ""
        # self.snap1_pyu_hmm = ""
        # self.maker1
        # self.snap2_pyu_hmm = ""
        #
        # self.maker_gff = ""


    def check_options(self):
        if not self.option("masked").is_set:
            raise OptionError("必须设置参数masked", code="22100701")
        if not self.option("ref_protein").is_set:
            raise OptionError("必须设置参数ref_protein", code="22100702")
        if not self.option("cegma_gff"):
            raise OptionError("必须设置参数cegma_gff", code="22100703")
        if not self.option("sample"):
            raise OptionError("必须设置参数sample", code="22100704")
        return True

    def run_genemark(self):
        self.genemark.set_options({
            "input_genome":self.option("masked")
        })
        self.genemark.on("end",self.run_snap1)
        #self.genemark_gmhmm_mod = os.path.join(self.genemark.output_dir,"gmhmm.mod")
        #self.logger.info(self.genemark_gmhmm_mod)
        self.genemark.run()

        #self.tool1.append(self.genemark)

    def run_snap1(self):
        self.snap1.set_options({
            "cegma_gff" : self.option("cegma_gff"),
            "fasta" : self.option("fasta"),
            "type" : "cegma"
        })
        self.snap1.on('end',self.run_maker1)
        self.snap1.run()
        #self.tool1.append(self.snap1)


    def run_maker1(self):
        self.maker1.set_options({
            "input_genome": self.option("masked"),
            "base" : self.option("sample"),
            "pyu_hmm" : os.path.join(self.snap1.output_dir,"pyu.hmm"),
            "gmhmm_mod" : os.path.join(self.genemark.output_dir,"gmhmm.mod"),
            "ref_protein" : self.option("ref_protein"),
            "get_fasta" : "F"

        })
        self.maker1.on('end',self.run_snap2)
        self.maker1.run()


    def run_snap2(self):
        self.maker_gff = os.path.join(self.maker1.output_dir,"{}.all.gff".format(self.option("sample")))
        self.snap2.set_options({
            "maker_gff": self.maker_gff,
            "type" : "maker"
        })
        self.snap2.on("end",self.run_augustus)
        self.snap2.run()


    def run_augustus(self):
        self.augustus_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        self.augustus.set_options({
            "input_genome": self.option("masked"),
            "species" : self.option("sample") + self.augustus_time,   # 解决存在相同样本名的问题
            "trainingset" : self.maker_gff
        })
        self.augustus.on("end",self.run_maker2)
        self.augustus.run()

    def run_maker2(self):
        self.maker2.set_options({
            #"input_genome": self.option("masked"),
            "input_genome": self.option("fasta"),
            "base" : self.option("sample"),
            "pyu_hmm" : os.path.join(self.snap2.output_dir,"pyu.hmm"),
            "gmhmm_mod" : os.path.join(self.genemark.output_dir,"gmhmm.mod"),
            "species" : self.option("sample") + self.augustus_time,
            "get_fasta" : "T"
        })
        self.maker2.on("end", self.run_rm_include)
        self.maker2.run()

    def run_rm_include(self):
        self.rm_include.set_options({
            "gff": "{1}/{0}.all.gff".format(self.option("sample"),self.maker2.output_dir),
            "faa": "{1}/{0}.proteins.fasta".format(self.option("sample"),self.maker2.output_dir),
            "sample": self.option("sample")
        })
        self.rm_include.on("end", self.run_format_stat)
        self.rm_include.run()

    def run_format_stat(self):
        self.format_stat.set_options({
            "gff": "{1}/{0}.all.gff".format(self.option("sample"),self.maker2.output_dir),
            "faa" : "{1}/{0}_rm_include.faa".format(self.option("sample"),self.rm_include.output_dir),
            "ffn" : "{1}/{0}.transcripts.fasta".format(self.option("sample"),self.maker2.output_dir),
            "genome":self.option("masked"),
            "sample": self.option("sample")

        })
        self.format_stat.on("end",self.set_output)
        self.format_stat.run()

    def set_output(self):
        new_files = "{1}/{0}_CDS.gff {1}/{0}_CDS.faa {1}/{0}_CDS.fnn".format(self.option("sample"),self.output_dir).split(" ")
        old_files = "{1}/{0}_CDS.gff {1}/{0}_CDS.faa {1}/{0}_CDS.fnn".format(self.option("sample"),self.format_stat.output_dir).split(" ")
        for file in new_files:
            if os.path.exists(file) :
                os.remove(file)
            os.link(old_files[new_files.index(file)],file)

        new_files2 = "{1}/{0}_CDS_length.xls {1}/{0}_CDS_statistics.xls".format(self.option("sample"),self.output_dir).split(" ")
        old_files2 = "{1}/{0}_CDS_length.xls {1}/{0}_CDS_statistics.xls".format(self.option("sample"),self.format_stat.output_dir).split(" ")
        for file in new_files2:
            if os.path.exists(file) :
                os.remove(file)
            os.link(old_files2[new_files2.index(file)],file)

        self.end()


    def run(self):
        super(FungiCdsModule, self).run()
        retstr = update_gm_key_pip(self.work_dir)  #20191118
        self.logger.info('update_gm_key_pip RESULT: %s'%retstr)
        self.run_genemark()

        # self.run_genemark()
        # self.run_snap1()
        # self.run_maker1()
        # self.run_augustus()
        # self.run_snap2()
        #
        # self.on_rely(self.tool1, self.maker1)
        # self.maker1.on('end',self.augustus)
        # self.augustus.on('end',self.self.snap2)
        # self.snap2('end',self.maker2)
        # for i in self.tool1:
        #     i.run()
        # self.maker1.run()
        # self.augustus.run()
        # self.snap2.run()
        # self.maker2.run()


    def end(self):
        super(FungiCdsModule, self).end()
