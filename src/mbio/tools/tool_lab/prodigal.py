# -*- coding: utf-8 -*-
# copy from prdict/prodigal


import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class ProdigalAgent(Agent):
    """
    prodigal 进行基因预测
    """

    def __init__(self, parent):
        super(ProdigalAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "sample", "type": "string","default":"out"},
            {"name": "trans_table","type":"string", "defalut":"11"},  #翻译的系统
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(ProdigalAgent, self).end()


class ProdigalTool(Tool):
    def __init__(self, config):
        super(ProdigalTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.prodigal_path = self.config.SOFTWARE_DIR +"/bioinfo/metaGenomic/Prodigal-2.6.3/prodigal"

    def run_prodigal(self):
        cmd = "{0} -i {2}  -o {1}.gff  -f gff -a {1}.faa -d {1}.ffn ".format(self.prodigal_path,self.option('sample'),self.genome_fasta)
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("prodigal运行完成")
        except subprocess.CalledProcessError:
            self.set_error("prodigal运行出错")
        num = 0
        with open(self.work_dir + '/' + self.option('sample') + '.gff',"r") as g:
            data = g.readlines()
            for i in data:
                if i.startswith("#") or i.strip() == "":
                    pass
                else:
                    num +=1
        return num
    def get_info(self):
        if os.path.exists(self.output_dir + '/' + self.option('sample')):
            os.remove(self.output_dir + '/' + self.option('sample'))
        os.mkdir(self.output_dir + '/' + self.option('sample'))
        all_gene_len = 0
        with open(self.output_dir + '/' + self.option('sample') + '/' + "Prodigal_Gene prediction_detail.xls","w") as f, open(self.work_dir + '/' + self.option('sample') + '.gff',"r") as v, open(self.option("input_genome").prop['path'],"r") as b:
            data1 = v.readlines()
            data2 = b.read()
            num = 0
            f.write("Gene ID\tSample Name\tLocation\tStrand\tStart\tEnd\tGene Len（bp）\tInitiator Codon\tTerminator Codon\n")
            for i in data1:
                if i.startswith("#"):
                    pass
                else:
                    num += 1
                    gene_id = "ORF" + str(num).zfill(4)
                    tmp = i.strip().split("\t")
                    all_gene_len += abs(int(tmp[3])-int(tmp[2]))
                    f.write(gene_id + "\t" + self.option('sample') + "\t" + tmp[0] + "\t" + tmp[5] + "\t" + tmp[2] + "\t" + tmp[3]
                            + "\t" + str(abs(int(tmp[3])-int(tmp[2]))) + "\t" + str(tmp[2].split("=")[-1]) + "\t" + str(tmp[3].split("=")[-1]) + "\n")


        all_len,gc_percent = self.get_gc(self.option('input_genome').prop["path"])
        with open(self.output_dir + '/' + self.option('sample') + '/' + "Prodigal_Gene prediction_stat.xls","w") as h:
            h.write("Sample Name\tGene No.\tGene Total Len（bp）\tGene Average Len（bp）\tGC Content  (%)\tGene/Genome (%)\n")
            h.write(self.option('sample') + "\t" + str(num) + "\t" + str(all_gene_len) + "\t" + str(all_gene_len/float(num))
                    + "\t" + str(gc_percent) + "\t" + str(all_gene_len/float(all_len)) + "\n")

    def set_output(self):
        gff_path = self.output_dir + '/' + self.option('sample') + '.gff'
        faa_path = self.output_dir + '/' + self.option('sample') + '.faa'
        ffn_path = self.output_dir + '/' + self.option('sample') + '.ffn'
        if os.path.exists(gff_path):
            os.remove(gff_path)
        if os.path.exists(faa_path):
            os.remove(faa_path)
        if os.path.exists(ffn_path):
            os.remove(ffn_path)
        os.link(self.work_dir + '/' + self.option('sample') + '.gff', gff_path)
        os.link(self.work_dir + '/' + self.option('sample') + '.faa', faa_path)
        os.link(self.work_dir + '/' + self.option('sample') + '.ffn', ffn_path)

    def get_gc(self,genome):
        GC_len=0; N_len=0; AT_len=0
        with open(genome,"r") as f:
            for line in f:
                if line.startswith(">"):
                    pass
                else:
                    value = line.strip()
                    GC_len += (value.count('g') + value.count('G') + value.count('c') + value.count('C'))
                    N_len += (value.count('n') + value.count('N'))
                    AT_len += (value.count('a') + value.count('A') + value.count('t') + value.count('T'))
        all_len = GC_len + N_len + AT_len
        gc_percent = round(GC_len / float(all_len) *100,2)
        return all_len,gc_percent

    def run(self):
        super(ProdigalTool, self).run()
        num = self.run_prodigal()
        if num > 0:
            self.get_info()
            self.set_output()
        self.end()
