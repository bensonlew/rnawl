#-*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# 20190903

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna.trans_step import step_count  # 对fasta序列做长度分布图输入表


class GenePredictAgent(Agent):
    """
    Glimmer3 进行细菌基因预测
    GenemarkS 进行质粒基因预测（plas）
    """

    def __init__(self, parent):
        super(GenePredictAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "genome_type", "type": "string", "default": ""},  # scaffold类型是细菌序列还是质粒序列("plas")
            {"name": "genome_gff", "type": "infile", "format": "gene_structure.gff3"},  #gff文件
            {"name": "gene_tag", "type": "string", "default": ""}, #基因前缀
            {"name": "seq", "type": "outfile", "format": "sequence.fasta"},  # 预测基因DNA序列
            {"name": "prot_seq", "type": "outfile", "format": "sequence.fasta"},  # 预测基因prot序列
            {"name": "orf_prefix", "type": "string", "default": "\"\""},  # ORF，自定义前缀
            {"name": "gene_statistics", "type": "outfile", "format": "sequence.profile_table"},  # 预测基因统计文件
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")
        if not self.option('genome_type'):
            raise OptionError("必须设置参数genome_type")
        if not self.option("genome_gff").is_set:
            raise OptionError("必须设置参数genome_gff")

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(GenePredictAgent, self).end()


class GenePredictTool(Tool):
    def __init__(self, config):
        super(GenePredictTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.glimmer_path = "/bioinfo/Genomic/Sofware/glimmer3.02/scripts/"
        self.transeq = "/bioinfo/seq/EMBOSS-6.6.0/emboss/transeq"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/lib"
        self.set_environ(LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)
        self.genemark_path = self.config.SOFTWARE_DIR + \
                             "/bioinfo/Genomic/Sofware/genemark_suite_linux_64/gmsuite/gmsn.pl"
        self.choose_seq = self.config.PACKAGE_DIR + "/bac_comp_genome/"
        self.python_path = "/miniconda2/bin/python"


    def get_gene_from_gff(self):
        """
        从上传的gff文件中提取出序列文件和整理文gff
        :return:
        """
        gff_path = self.option('genome_gff').prop['path']
        cmd = '{} {}get_gff_from_gff.py -i {} -tag {} -type {} -genome {} -out {}'.format(self.python_path, self.package_path,
                                    gff_path, self.option('gene_tag'), 'cds', self.option('gene_type'), self.work_dir + '/gene_out')
        self.logger.info(cmd)
        command = self.add_command('get_gff', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("get_gff运行完成")
        else:
            self.set_error("get_gff运行出错!")

    def run_choose_seq(self):
        """
        从原始gff文件中生成faa、fnn和新的gff文件(增加起始密码子和终止密码子)
        :return:
        """

        faa_name = self.option('sample') + '.gene.faa'
        faa_path = os.path.join(self.db_path, self.option('sample') + faa_name)
        if os.path.exists(self.work_dir + '/out_seq'):
            shutil.rmtree(self.work_dir + '/out_seq')
        os.mkdir(self.work_dir + '/out_seq')
        if self.option('genome_type') == 'complete' or self.option('genome_type') == 'chromosome':
            gff_path = os.path.join(self.work_dir, self.option('genome') + file)
            fna_path = os.path.join(self.db_path, self.option('genome') + '.fna')
        else:
            gff_path = os.path.join(self.work_dir, self.option('sample') + file)
            fna_path = os.path.join(self.db_path, self.option('sample') + '.fna')
        cmd = '{} {}get_fa_bygff.pl {} {} {} {} {}'.format(self.perl_path, self.package_path,
                    fna_path, faa_path, gff_path, self.work_dir + '/out_seq', self.option('genome_type'))
        self.logger.info(cmd)
        command_name = "gff" + '_' + str(gff_num)
        command = self.add_command(command_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取序列完成运行完成")
        else:
            self.set_error("提取序列运行出错!")

    def run_glimmer(self):
        cmd = "{}g3-iterated.csh {} {}".format(self.glimmer_path, self.genome_fasta, self.sample_name)
        command = self.add_command("genepredict", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Glimmer运行完成")
            self.run_glimmer2ffn()
        else:
            self.set_error("Glimmer运行出错!", code="33300501")

    def run_glimmer2ffn(self):
        predict_file = self.work_dir + "/" + self.sample_name + ".predict"
        cmd = '{} {}glimmer2ffn.pl {} {} {} 0 {}'.format(self.perl_path, self.perl_script, self.genome_fasta,
                                                         predict_file, self.option("orf_prefix"), self.sample_name)
        command = self.add_command("run_glimmer2ffn", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn运行完成")
            if os.path.getsize(self.work_dir + "/" + self.sample_name + ".fnn"):
                self.run_transeq()
            else:
                self.end()
        else:
            self.set_error("提取fnn运行出错!", code="33300502")

    def run_transeq(self):
        nul_seq = self.work_dir + "/" + self.sample_name + ".fnn"
        port_seq = self.work_dir + "/" + self.sample_name + ".faa"
        cmd = '{} -trim -table 11 -sequence {} -outseq {}'.format(self.transeq, nul_seq, port_seq)
        command = self.add_command("transeq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("翻译蛋白序列运行完成")
            self.run_get_information()
        else:
            self.set_error("翻译蛋白序列运行出错!", code="33300503")

    def run_get_information(self):
        nul_seq = self.work_dir + "/" + self.sample_name + ".fnn"
        port_seq = self.work_dir + "/" + self.sample_name + ".faa"
        cmd = '{} {}trim_gene_info.pl {} {} {} {} {}'.format(self.perl_path, self.perl_script, self.genome_fasta,
                                                             nul_seq, port_seq, self.option("orf_prefix"),
                                                             self.output_dir)
        command = self.add_command("get_information", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("统计预测结果运行完成")

        else:
            self.set_error("统计预测结果运行出错!", code="33300504")

    def run_genemark(self):
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
            self.set_error("genemark运行出错!", code="33300505")

    def genemark_predict(self):
        predict_file = self.work_dir + "/" + self.sample_name
        predict_put = self.output_dir + "/" + self.sample_name + ".predict.gff"
        gene_id = ""
        sequence_id = ""
        sample = ""
        scaff = ""
        with open(predict_file, "r") as f, open(predict_put, "wb") as w:
            w.write("Gene id\tSequence id\tStart\tEnd\tStrand\tGene Length(bp)\tProtein Length\n")
            for line in f:
                if re.match("^FASTA definition line:\s(\w+)$", line):
                    sample_scaff = line.strip().split(' ')[-1]
                    li = sample_scaff.split('_')
                    scaff = li.pop()
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
                        sequence_id = scaff + "_" + self.option("orf_prefix") + "0" * (4 - len(arr[0])) + arr[0]
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
        gli_path = self.work_dir + "/" + self.sample_name + ".predict.gff"
        if self.option("genome_type") == "" and os.path.exists(gli_path):
            if os.path.exists(self.output_dir + "/" + self.sample_name + ".predict.gff"):
                os.remove(self.output_dir + "/" + self.sample_name + ".predict.gff")
            os.link(gli_path, self.output_dir + "/" + self.sample_name + ".predict.gff")
        nul_path = self.output_dir + "/" + self.sample_name + ".fnn"
        port_path = self.output_dir + "/" + self.sample_name + ".faa"
        statistics_path = self.output_dir + "/gene_statistics.xls"
        if os.path.getsize(nul_path):
            self.option('seq', nul_path)
            self.option('prot_seq', port_path)
            self.logger.info(step_count(self.option('seq').prop['path'], self.work_dir + "/fnn_stat.xls", 11, 100,
                                        self.output_dir + "/length_distribute.txt"))
            step_count(self.option('seq').prop['path'], self.work_dir + "/fnn_stat.xls", 11, 100,
                       self.output_dir + "/length_distribute.txt")
        self.option('gene_statistics', statistics_path)
        self.option('gff', self.output_dir + "/" + self.sample_name + ".predict.gff")

    def run(self):
        super(GenePredictTool, self).run()
        if self.option("genome_type") == "":
            self.run_glimmer()
        else:
            self.run_genemark()
        self.set_output()
        self.end()





