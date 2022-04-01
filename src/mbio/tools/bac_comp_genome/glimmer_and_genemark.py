# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 20190918


import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import sort_seq
from mbio.packages.ref_rna.trans_step import step_count  # 对fasta序列做长度分布图输入表


class GlimmerAndGenemarkAgent(Agent):
    """
    #比较基因组进行基因预测
    Glimmer3 进行细菌基因预测
    GenemarkS 进行质粒基因预测（plas）
    """

    def __init__(self, parent):
        super(GlimmerAndGenemarkAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "genome_type", "type": "string", "default": ""},  # scaffold类型是细菌序列还是质粒序列("plas")
            {"name": "seq", "type": "outfile", "format": "sequence.fasta"},  # 预测基因DNA序列
            {"name": "prot_seq", "type": "outfile", "format": "sequence.fasta"},  # 预测基因prot序列
            {"name": "orf_prefix", "type": "string", "default": "\"\""},  # ORF，自定义前缀
            {"name": "gene_statistics", "type": "outfile", "format": "sequence.profile_table"},  # 预测基因统计文件
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式
            {"name": "software", "type":"string","default":"glimmer"},  #选择的软件： glimmer， genemark
            {"name": "genome_name", "type" : "string"}, #基因组名称用于以区别组装序列的名称
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(GlimmerAndGenemarkAgent, self).end()


class GlimmerAndGenemarkTool(Tool):
    def __init__(self, config):
        super(GlimmerAndGenemarkTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = self.option('genome_name')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bac_comp_genome/"
        self.glimmer_path = "/bioinfo/Genomic/Sofware/glimmer3.02/scripts/"
        self.transeq = "/bioinfo/seq/EMBOSS-6.6.0/emboss/transeq"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/lib"
        self.set_environ(LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)
        self.genemark_path = self.config.SOFTWARE_DIR + \
                             "/bioinfo/Genomic/Sofware/genemark_suite_linux_64/gmsuite/gmsn.pl"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/lib"
        self.set_environ(LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)

    def run_glimmer(self):
        """
        glimmer软件进行预测
        :return:
        """
        sort_seq(self.genome_fasta, self.work_dir + "/seq.fasta")
        cmd = "{}g3-iterated.csh {} {}".format(self.glimmer_path, self.work_dir + "/seq.fasta", self.sample_name)
        command = self.add_command("genepredict", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Glimmer运行完成")
            self.run_glimmer2ffn()
        else:
            self.set_error("Glimmer运行出错!")

    def run_glimmer2ffn(self):
        """
        根据预测结果提取出gff文件和fnn序列
        :return:
        """
        predict_file = self.work_dir + "/" + self.sample_name + ".predict"
        cmd = '{} {}glimmer2ffn.pl {} {} {} 0 {}'.format(self.perl_path, self.perl_script, self.work_dir + "/seq.fasta",
                                                         predict_file, self.option("orf_prefix"), self.sample_name)
        command = self.add_command("run_glimmer2ffn", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn运行完成")
            if os.path.getsize(self.work_dir + "/" + self.sample_name + ".fna"):
                self.run_transeq()
            else:
                self.end()
        else:
            self.set_error("提取fnn运行出错!")

    def run_transeq(self):
        """
        对核酸文件进行翻译
        :return:
        """
        nul_seq = self.work_dir + "/" + self.sample_name + ".fna"
        port_seq = self.work_dir + "/" + self.sample_name + ".faa"
        cmd = '{} -trim -table 11 -sequence {} -outseq {}'.format(self.transeq, nul_seq, port_seq)
        command = self.add_command("transeq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("翻译蛋白序列运行完成")
            self.run_get_information()
        else:
            self.set_error("翻译蛋白序列运行出错!")

    def run_get_information(self):
        """
        根据核酸文件和蛋白文件完成统计
        :return:
        """
        if not os.path.exists(self.work_dir + "/seq.fasta"):
            os.link(self.genome_fasta, self.work_dir + "/seq.fasta")
        else:
            os.remove(self.work_dir + "/seq.fasta")
            os.link(self.genome_fasta, self.work_dir + "/seq.fasta")
        nul_seq = self.work_dir + "/" + self.sample_name + ".fna"
        if os.path.exists(nul_seq):
            nul_seq = nul_seq
        else:
            if os.path.exists(self.work_dir + "/" + self.sample_name + ".fnn"):
                os.link(self.work_dir + "/" + self.sample_name + ".fnn", self.work_dir + "/" + self.sample_name + ".fna")
                nul_seq = self.work_dir + "/" + self.sample_name + ".fna"
            else:
                self.set_error("未能成功的生成核酸结果，请检查!")
        port_seq = self.work_dir + "/" + self.sample_name + ".faa"
        cmd = '{} {}trim_gene_info.pl {} {} {} {} {}'.format(self.perl_path, self.perl_script, self.work_dir + "/seq.fasta",
                                                             nul_seq, port_seq, self.option("orf_prefix"),
                                                             self.output_dir)
        command = self.add_command("get_information", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("统计预测结果运行完成")
        else:
            self.set_error("统计预测结果运行出错!")

    def run_genemark(self):
        """
        对质粒进行基因预测（genemark）
        :return:
        """
        os.system("cp ~/.gm_key ./")
        cmd = '{} {} {} --output {} --fnn --faa'.format(self.perl_path, self.genemark_path, self.genome_fasta,
                                                        self.sample_name)
        command = self.add_command("run_genemark", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("genemark运行完成")
            self.run_get_information()
            self.genemark_predict()
        else:
            self.set_error("genemark运行出错!")

    def genemark_predict(self):
        """
        生成标准格式的gff文件
        :return:
        """
        self.logger.info("开始对predict.gff进行标准化")
        predict_file = self.work_dir + "/" + self.sample_name
        predict_put = self.output_dir + "/" + self.sample_name + "_CDS.gff"
        gene_id = ""
        sequence_id = ""
        sample = ""
        scaff = ""
        with open(predict_file, "r") as f, open(predict_put, "wb") as w:
            w.write("Gene id\tSequence id\tStart\tEnd\tStrand\tGene Length(bp)\tProtein Length\n")
            for line in f:
                if re.match(r"FASTA definition line", line):
                    sample_scaff = line.strip().split(' ')[-1]
                    scaff = sample_scaff.strip()
                    li = sample_scaff.split('_')
                    # scaff = li.pop()
                    # scaff = li.pop().split('Plasmid')
                    # scaff = "p" + scaff[-1]
                    sample = "_".join(li)
                    nu = 0
                if re.match("^\s+[0-9]+.*$", line):
                    info = re.subn("\s+", "\t", line.strip('\n'))[0]
                    arr = info.lstrip('\t').split('\t')
                    if len(self.option("orf_prefix").strip("\"")) != 0:
                        gene_id = self.option("orf_prefix") + "0" * (4 - len(arr[0])) + arr[0]
                        nu = nu + 1
                        arr[0] = str(nu)
                        sequence_id = scaff + "_ORF" + "0" * (4 - len(arr[0])) + arr[0]
                    else:
                        gene_id = sample + "_ORF" + "0" * (4 - len(arr[0])) + arr[0]
                        nu = nu + 1
                        arr[0] = str(nu)
                        sequence_id = scaff + "_ORF" + "0" * (4 - len(arr[0])) + arr[0]
                    arr[0] = gene_id + "\t" + sequence_id
                    strand = arr[1]
                    arr[1] = arr[2].split('>')[-1].split('<')[-1]
                    arr[2] = arr[3].split('>')[-1].split('<')[-1]
                    arr[3] = strand
                    arr[5] = str(int(arr[4]) / 3 - 1)
                    w.write("\t".join(arr) + "\n")
                else:
                    pass

    def set_output(self):
        self.logger.info("开始设置结果文件目录")
        gli_path = self.work_dir + "/" + self.sample_name + ".predict.gff"
        if self.option("software") == "glimmer" and os.path.exists(gli_path):
            if os.path.exists(self.output_dir + "/" + self.sample_name + "_CDS.gff"):
                os.remove(self.output_dir + "/" + self.sample_name + "_CDS.gff")
            os.link(gli_path, self.output_dir + "/" + self.sample_name + "_CDS.gff")
        nul_path = self.output_dir + "/" + self.sample_name + ".fna"
        port_path = self.output_dir + "/" + self.sample_name + ".faa"
        self.option('seq',nul_path)
        self.option('prot_seq', port_path)
        self.option('gff', self.output_dir + "/" + self.sample_name + "_CDS.gff")

    def run(self):
        self.logger.info("开始运行tool！")
        super(GlimmerAndGenemarkTool, self).run()
        if self.option("software") == "glimmer":
            self.run_glimmer()
        else:
            self.run_genemark()
        self.set_output()
        self.end()
