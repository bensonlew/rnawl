# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import ConfigParser
import unittest
import shutil
import subprocess
import glob
from Bio import SeqIO
from skbio.parse.sequences import parse_fasta


class NovelMirnaNewAgent(Agent):
    """
    新miRNA预测
    """

    def __init__(self, parent):
        super(NovelMirnaNewAgent, self).__init__(parent)
        options = [
            {"name": "category", "type": "string", "default": ""},  # 物种分类，Animal or Plant
            {"name": "species", "type": "string", "default": ""},  # 具体物种
            {"name": "input_genome", "type": "infile", "format": "small_rna.fasta"},  # 参考基因组文件
            {"name": "input_fa", "type": "infile", "format": "small_rna.fasta"},  # 输入序列文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "mapping_arf", "type": "infile", "format": "small_rna.common"},  # clean reads mapping到参考基因组的结果
            {"name": "mrd", "type": "infile", "format": "small_rna.common"},  # 已知miRNA定量结果文件
            {"name": "known_miRNA_express", "type": "infile", "format": "small_rna.common"},  # 已知miRNA定量结果
            {"name": "clean_fa", "type": "infile", "format": "small_rna.fasta"},  # 用来mapping到参考基因组的fasta序列
            {"name": "novel_count", "type": "outfile", "format": "small_rna.common"},  # 新miRNA表达count表
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"},  # 鉴定完新miRNA的过滤文件
            {'default': '', 'type': 'string', 'name': 'index'},  # bowtie索引
            {'default': 'mireap', 'type': 'string', 'name': 'software'},  # 植物预测miRNA所用软件
            {'default': 1, 'type': 'float', 'name': 'rpm'},  # mirdeep-p2参数，用于测序数据过滤
        ]
        self.add_option(options)
        self.step.add_steps("novel_mirna")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.novel_mirna.start()
        self.step.update()

    def stepfinish(self):
        self.step.novel_mirna.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if self.option("category").lower() not in ["animal", "plant"]:
            raise OptionError("物种分类只能为Animal/Plant")
        if self.option("species") == "null":
            raise OptionError("必须指定具体物种")
        if not self.option("input_fa").is_set:
            raise OptionError("必须提供质控后的FASTA文件")
        if not self.option("mapping_arf").is_set:
            raise OptionError("必须提供clean reads mapping到参考基因组的结果文件")
        if not self.option("config").is_set:
            raise OptionError("必须提供配置文件")
        with open(self.option("input_fa").prop["path"], "r") as f:
            line = f.readline()
            if not re.compile(r'^>\S+_x\d+').match(line):
                raise OptionError("质控后的序列文件格式不对")

        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '30G'

    def end(self):
        super(NovelMirnaNewAgent, self).end()


class NovelMirnaNewTool(Tool):
    def __init__(self, config):
        super(NovelMirnaNewTool, self).__init__(config)
        self.mirdp2_bash = "bioinfo/miRNA/miRDP2-v1.1.4/miRDP2-v1.1.4_pipeline_modify.bash"
        self.mirdp2_scripts = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/miRDP2-v1.1.4/scripts/"
        self.mirdeep = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mirdeep2/"
        self.mireap = self.config.SOFTWARE_DIR + "/bioinfo/miRNA/mireap_0.2/bin/"
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.1.2/'
        self.perl = 'program/perl-5.24.0/bin/'
        self.perl_old = '/usr/bin/'
        self.rnafold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Progs/'
        self.parse_arf = self.config.PACKAGE_DIR + "/small_rna/"
        self.FFW1 = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mireap_0.2/lib/'
        self.RNA = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Perl/blib/lib/'
        self.lib_rna = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/ViennaRNA-1.8.4/Perl/blib/arch/auto/RNA/'
        self.essential_lib = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/lib/site_perl/5.24.0/'
        self.squid = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/squid-1.9g/'
        self.randfold = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/mirdeep2/essentials/randfold-2.0/'
        self.set_environ(LD_LIBRARY_PATH=self.lib_rna)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.perl)
        # self.set_environ(PATH=self.perl_old)
        self.set_environ(PERL5LIB=self.FFW1)
        self.set_environ(PERL5LIB=self.RNA)
        if self.config.SOFTWARE_DIR == "/mnt/ilustre/users/sanger-dev/app":
            self.set_environ(PERL5LIB=self.essential_lib)
        self.set_environ(PATH=self.parse_arf)
        self.set_environ(PATH=self.mirdeep)
        self.set_environ(PATH=self.randfold)
        self.set_environ(PATH=self.squid)
        self.set_environ(PATH=self.mireap)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.rnafold)
        self.set_environ(PATH=self.mirdp2_bash)
        self.set_environ(PATH=self.mirdp2_scripts)

    def run(self):
        """
        运行
        :return:
        """
        super(NovelMirnaNewTool, self).run()
        self.get_mirna()
        if self.option("category").lower() == "animal":
            self.run_remove_white_space_in_id()
            self.run_replace_nonatcgun()
            self.run_mirdeep2()
            self.parse_mirdeep_result()
        elif self.option("category").lower() == "plant" and self.option("software") == "mireap":
            self.parse_mapping_file()
            self.run_mireap()
            self.parse_mireap_results()
            self.rna_fold()
            self.ps_2_pdf()
            self.set_mireap_results()
        else:
            self.run_mirdp2()
            self.parse_mirdp2_results()
            self.set_mirdp2_results()
        self.end()

    def get_mirna(self):
        species = self.option("species").split(",")
        self.mature_fa = os.path.join(self.work_dir, species[0] + "_mature.fa")
        self.hairpin_fa = os.path.join(self.work_dir, species[0] + "_hairpin.fa")
        mature_i = self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.mature.fa".format(species[0], species[0])
        hairpin_i = self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.hairpin.fa".format(species[0],
                                                                                                   species[0])
        if os.path.exists(self.mature_fa):
            os.remove(self.mature_fa)
        if os.path.exists(self.hairpin_fa):
            os.remove(self.hairpin_fa)
        os.link(mature_i, self.mature_fa)
        os.link(hairpin_i, self.hairpin_fa)

        self.other_mature_fa = os.path.join(self.work_dir, "other_mature.fa")
        ## 预测新miRNA时，把其他所有的物种作为参考
        with open(self.other_mature_fa, "w") as w1:
            mature_fa_file = []
            for i in species:
                mature_fa_file.extend(
                    [self.config.SOFTWARE_DIR + "/database/mirbase/species/{}/{}.mature.fa".format(i, i)])
                self.logger.info(mature_fa_file)
            other_fa_file = glob.glob(self.config.SOFTWARE_DIR + "/database/mirbase/species/*/*.mature.fa")
            for file in other_fa_file:
                if file not in mature_fa_file:
                    with open(file, "r") as r1:
                        lines = r1.readlines()
                        w1.writelines(lines)

    def parse_mirdeep_result(self):
        self.logger.info("处理mirdeep2预测结果")
        file = os.path.basename(glob.glob(self.work_dir + "/miRNAs_expressed_all_samples_*.csv")[0])
        id = re.search("miRNAs_expressed_all_samples_(.*).csv", file).group(1)
        cmd = "{}perl {}convent_mirdeep_result_modify.pl -id {} -o mirdeep -config {}".format(self.perl, self.mirdeep,
                                                                                              id, self.option(
                "config").prop["path"])
        command = self.add_command("parse_mirdeep_result", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("parse_mirdeep_result完成")
        else:
            self.set_error("parse_mirdeep_result出错")

        infile1 = self.work_dir + "/mirdeep/predicted_total_sumary_table.xls"
        with open(infile1, "r") as f, open(self.work_dir + "/novel_miR_count.xls", "w") as w:
            for line in f:
                if re.search(r'novel_miRNAs', line):
                    continue
                elif re.search(r'miRBase_miRNAs', line):
                    break
                elif re.search(r'#provisional_id', line):
                    items = line.strip().split("\t")
                    w.write("miRNA_ID" + "\t" + "\t".join(items[2:]) + "\n")
                else:
                    items = line.strip().split("\t")
                    w.write(items[0] + "\t" + "\t".join(items[2:]) + "\n")
        '''
        cmd1 = """sed 1d %s | sed '$d' | sed '$d' | sed 's/#provisional_id/miRNA_ID/' |cut -f 1,3- > novel_miR_count.xls""" % (infile1)
        try:
            subprocess.check_output(cmd1, shell=True)
            self.logger.info('生成novel_miR_count文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成novel_miR_count文件失败')
            self.set_error("生成novel_miR_count文件失败")
        '''

        infile2 = self.work_dir + "/mirdeep/novo_mature_normalization.xls"
        cmd2 = """sed 's/#provisional_id/miRNA_ID/' %s > novel_miR_norm.xls""" % (infile2)
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info('生成novel_miR_norm文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成novel_miR_norm文件失败')
            self.set_error("生成novel_miR_norm文件失败")

        mirna_detail = os.path.join(self.work_dir, "novel_mirna_detail")
        result_bed = glob.glob(self.work_dir + "/result_*.bed")
        precursors_str = glob.glob(self.work_dir + "/mirdeep_runs/run_*/tmp/precursors.str")[0]
        novel_precursors_str = os.path.join(self.work_dir, "precursors.str")
        novel_mature_fa = {}
        novel_precursor_fa = {}
        precursors_str_dict = {}
        for seq_record in SeqIO.parse(self.work_dir + "/mirdeep/novo_mature_seq.fa", "fasta"):
            novel_mature_fa[seq_record.id] = seq_record.seq
        for seq_record in SeqIO.parse(self.work_dir + "/mirdeep/novo_precursor_seq.fa", "fasta"):
            novel_precursor_fa[seq_record.id] = seq_record.seq
        with open(precursors_str, "r") as f:
            seq_id = {}
            id = ""
            sequence = ""
            for line in f:
                if line.startswith(">"):
                    id = line.strip().split(" ")[0].split(">")[1]
                    seq_id[id] = line.strip()
                    precursors_str_dict[id] = ""
                else:
                    if id != "":
                        precursors_str_dict[id] += line
        with open(result_bed[0], "r") as f1, open(mirna_detail, "w") as w1, open(novel_precursors_str, "w") as w2:
            head_list = ["miRNA_id", "miRNA_seq", "miRNA_len", "pre_name", "pre_seq", "pre_len", "pre_position",
                         "score"]
            w1.write("\t".join(head_list) + "\n")
            for line in f1:
                items = line.strip().split("\t")
                if re.match(r'novel', items[3]):
                    id = items[3].split(":")[1]
                    pos = items[0] + "(" + items[5] + "):" + items[6] + "-" + items[7]
                    score = items[4]
                    if id in novel_mature_fa.keys():
                        mature_seq = novel_mature_fa[id]
                        precursor_seq = novel_precursor_fa[id]
                        w1.write("\t".join(
                            [id, str(mature_seq), str(len(mature_seq)), id, str(precursor_seq), str(len(precursor_seq)),
                             pos, score]) + "\n")
                        w2.write(">" + id + "\n" + precursors_str_dict[id])

        mirna_count = os.path.join(self.work_dir, "novel_miR_count.xls")
        mirna_norm = os.path.join(self.work_dir, "novel_miR_norm.xls")
        mirna_detail = os.path.join(self.work_dir, "novel_mirna_detail")
        mrd = glob.glob(self.work_dir + "/mirdeep_runs/run_*/output.mrd")[0]
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_count.xls"))
        os.link(mirna_count, os.path.join(self.output_dir, "novel_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_norm.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        os.link(mirna_norm, os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_detail.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        os.link(mirna_detail, os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        if os.path.exists(os.path.join(self.output_dir, "miRBase.mrd")):
            os.remove(os.path.join(self.output_dir, "miRBase.mrd"))
        os.link(mrd, os.path.join(self.output_dir, "miRBase.mrd"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursors.str")):
            os.remove(os.path.join(self.output_dir, "novel_precursors.str"))
        os.link(novel_precursors_str, os.path.join(self.output_dir, "novel_precursors.str"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mature_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_mature_seq.fa"))
        os.link(self.work_dir + "/mirdeep/novo_mature_seq.fa", os.path.join(self.output_dir, "novel_mature_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursor_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        os.link(self.work_dir + "/mirdeep/novo_precursor_seq.fa",
                os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        if os.path.exists(self.output_dir + "/structure_pdf"):
            shutil.rmtree(self.output_dir + "/structure_pdf")
        os.mkdir(self.output_dir + "/structure_pdf")
        self.option("novel_count", os.path.join(self.output_dir, "novel_mirna_count.xls"))
        pdf_files = glob.glob(self.work_dir + "/mirdeep/pdfs_*/*.pdf")
        for file in pdf_files:
            base_name = os.path.basename(file).split(".pdf")[0]
            if base_name in novel_mature_fa.keys():
                os.link(file, self.output_dir + "/structure_pdf/" + os.path.basename(file))
        ## 整理鉴定完novel miRNA剩余的reads
        mrd = os.path.join(self.output_dir, "miRBase.mrd")
        filter_fa = os.path.join(self.output_dir, "filtered.fa")
        input_fa = self.option("input_fa").prop["path"]
        a = dict()
        with open(mrd, "r") as f1:
            for line in f1:
                if re.compile(r'\S+_\d+_x\d+\s+').match(line):
                    a[line.strip().split()[0]] = 1
        with open(filter_fa, "w") as w:
            for seq_id, seq_sequence in parse_fasta(input_fa):
                if seq_id not in a:
                    w.write('>%s\n%s\n' % (seq_id, seq_sequence))
        self.option("filter_fa", filter_fa)
        ## 去掉PDF转换为PNG步骤
        '''
        self.logger.info("开始PDF转换为PNG")
        pdf_paths = [os.path.join(self.output_dir + '/structure_pdf', pdf) for pdf in
                     os.listdir(self.output_dir + '/structure_pdf') if re.match(r'^\S+\.pdf$', pdf.strip())]
        self.logger.info(pdf_paths)
        if len(pdf_paths) < 1:
            self.logger.info('未生成pdf文件：路径为 %s' % self.output_dir + '/structure_pdf')
        else:
            for i in range(0, len(pdf_paths)):
                pdf_paths_new = pdf_paths[i]
                re_string = re.compile('.pdf$')
                png_paths = re_string.sub(".png", pdf_paths_new)
                self.convert_pdf_to_png(i+1, pdf_paths_new, png_paths)
        '''

    def convert_pdf_to_png(self, num, olds, news):
        self.image_magick = '/program/ImageMagick/bin/convert'
        cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + olds + ' ' + news
        self.logger.info(cmd)
        command = self.add_command('convert_{}'.format(num), cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            pass
        else:
            self.set_error("PDF转换PNG出错!")

    def parse_mapping_file(self):
        self.logger.info("开始对mapping arf文件进行处理")
        cmd = "perl {}parse_arf.pl -rfam {} -mrd {} -arf {}".format(self.parse_arf,
                                                                    self.option("clean_fa").prop["path"],
                                                                    self.option("mrd").prop["path"],
                                                                    self.option("mapping_arf").prop["path"])
        command = self.add_command("parse_arf", cmd)
        command.software_dir = "/usr/bin/"
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("对mapping arf文件处理完成")
        else:
            self.set_error("对mapping arf文件处理出错")

    def parse_mireap_results(self):
        self.logger.info("开始对mireap结果文件进行处理")
        cmd = "{}perl {}parse_mireap_result_modify.pl -exp {} -aln {} -gff {} -config {} -o {}".format(
            self.perl, self.mireap, self.option("known_miRNA_express").prop["path"],
            os.path.join(self.work_dir, "mireap-Nov.aln"),
            os.path.join(self.work_dir, "mireap-Nov.gff"), self.option("config").prop["path"],
            os.path.join(self.work_dir, "mireap"))
        command = self.add_command("parse_mireap_results", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("对mireap结果文件处理完成")
        else:
            self.set_error("对mireap结果文件处理出错")

    def run_remove_white_space_in_id(self):
        input_genome = self.option("input_genome").prop["path"]
        self.parse_genome = os.path.join(self.work_dir, os.path.basename(input_genome))
        cmd = "perl {}remove_white_space_in_id.pl {} > {}".format(self.mirdeep, input_genome, self.parse_genome)
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('remove_white_space_in_id计算完成')
        except subprocess.CalledProcessError:
            self.logger.info('remove_white_space_in_id计算失败')
            self.set_error("remove_white_space_in_id运行失败")

    def run_replace_nonatcgun(self):
        temp_mature_fa = os.path.join(self.work_dir, "temp_mature.fa")
        temp_hairpin_fa = os.path.join(self.work_dir, "temp_hairpin.fa")
        temp_other_mature_fa = os.path.join(self.work_dir, "temp_other_mature.fa")
        cmd1 = "perl {}remove_nonatcg_in_seq.pl {} > {}".format(self.mirdeep, self.mature_fa, temp_mature_fa)
        self.logger.info(cmd1)
        try:
            subprocess.check_output(cmd1, shell=True)
            self.logger.info('replace_nonatcgun计算完成')
        except subprocess.CalledProcessError:
            self.set_error("replace_nonatcgun运行失败")
        cmd2 = "perl {}remove_nonatcg_in_seq.pl {} > {}".format(self.mirdeep, self.hairpin_fa, temp_hairpin_fa)
        self.logger.info(cmd2)
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info('replace_nonatcgun计算完成')
        except subprocess.CalledProcessError:
            self.set_error("replace_nonatcgun运行失败")
        cmd3 = "perl {}remove_nonatcg_in_seq.pl {} > {}".format(self.mirdeep, self.other_mature_fa,
                                                                temp_other_mature_fa)
        self.logger.info(cmd3)
        try:
            subprocess.check_output(cmd3, shell=True)
            self.logger.info('replace_nonatcgun计算完成')
        except subprocess.CalledProcessError:
            self.set_error("replace_nonatcgun运行失败")
        self.mature_fa = temp_mature_fa
        self.hairpin_fa = temp_hairpin_fa
        self.other_mature_fa = temp_other_mature_fa

    def run_mirdeep2(self):
        """
        动物novel miRNA预测
        """
        tmp_files = glob.glob(self.work_dir + "/miRNAs_expressed_all_samples_*.csv") + glob.glob(
            self.work_dir + "/result_*.bed") + glob.glob(self.work_dir + "/mirdeep_runs/run_*/tmp/precursors.str")
        for file in tmp_files:
            os.remove(file)
        self.logger.info("开始进行novel miRNA预测")
        input_fa = self.option("input_fa").prop["path"]
        mapping_arf = self.option("mapping_arf").prop["path"]
        cmd = "{}perl {}miRDeep2_modify.pl {} {} {} {} {} {}".format(self.perl, self.mirdeep, input_fa,
                                                                     self.parse_genome, mapping_arf,
                                                                     self.mature_fa, self.other_mature_fa,
                                                                     self.hairpin_fa)
        command = self.add_command("mirdeep2", cmd)
        # command.software_dir = "/usr/bin/"
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行mirdeep脚本完成")
        else:
            self.set_error("运行mirdeep脚本出错")

    def run_mireap(self):
        """
        植物novel miRNA预测
        """
        self.logger.info("开始进行novel miRNA预测")
        map = self.work_dir + "/map.txt"
        input_genome = self.option("input_genome").prop["path"]
        cmd = "perl {}mireap.pl -i {} -m {} -r {} -t Nov".format(self.mireap,
                                                                 os.path.join(self.work_dir, "filtered.fa"), map,
                                                                 input_genome)
        command = self.add_command("mireap", cmd)
        command.software_dir = "/usr/bin/"
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行mireap脚本完成")
        else:
            self.set_error("运行mireap脚本出错")

    def run_mirdp2(self):
        """
        植物novel miRNA预测
        """
        input_genome = self.option("input_genome").prop["path"]
        mirdp2_output = self.work_dir + "/mirdp2"
        os.mkdir(mirdp2_output)
        cmd = "{} -g {} -x {} -i {} -o {} -R {}".format(self.mirdp2_bash, input_genome, self.option("index"),
                                                            self.option("input_fa").prop["path"], mirdp2_output,self.option("rpm"))
        command = self.add_command("mirdp2", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行mirdp2脚本完成")
        else:
            self.set_error("运行mirdp2脚本出错")

    def parse_mirdp2_results(self):
        self.logger.info("开始对mirdp2结果文件进行处理")
        prediction_file = glob.glob(self.work_dir + "/mirdp2/*/*filter_P_prediction")[0]
        mirdp2_result = os.path.join(self.work_dir, "mirdp2_result")
        os.mkdir(mirdp2_result)
        cmd = "{}perl {}parse_mirdp2_result.pl -exp {} -prediction {} -config {} -o {}".format(
            self.perl, self.mirdp2_scripts, self.option("known_miRNA_express").prop["path"],
            prediction_file, self.option("config").prop["path"], mirdp2_result)
        command = self.add_command("parse_mirdp2_results", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("对mirdp2结果文件处理完成")
        else:
            self.set_error("对mirdp2结果文件处理出错")

        with open(prediction_file, "r") as f:
            for line in f:
                pre_miRNA = str(line.strip().split("\t")[3])[0:12]
                ps = pre_miRNA + "_ss.ps"
                pdf = pre_miRNA + ".pdf"
                cmd1 = "ps2pdf {} {}".format(ps, pdf)
                try:
                    subprocess.check_output(cmd1, shell=True)
                    self.logger.info('ps2pdf转换完成')
                except:
                    self.set_error("ps2pdf转换失败")

    def set_mirdp2_results(self):
        self.logger.info("开始设置mirdp2结果目录")
        ## novel miRNA二级结构图片
        if os.path.exists(os.path.join(self.output_dir, "structure_pdf")):
            shutil.rmtree(os.path.join(self.output_dir, "structure_pdf"))
        if os.path.exists(os.path.join(self.work_dir, "structure_ps")):
            shutil.rmtree(os.path.join(self.work_dir, "structure_ps"))
        os.mkdir(os.path.join(self.output_dir, "structure_pdf"))
        os.mkdir(os.path.join(self.work_dir, "structure_ps"))
        pdf = glob.glob(self.work_dir + "/*.pdf")
        ps = glob.glob(self.work_dir + "/*_ss.ps")
        for file in pdf:
            shutil.move(file, os.path.join(self.output_dir, "structure_pdf"))
        for file in ps:
            shutil.move(file, os.path.join(self.work_dir, "structure_ps"))

        novel_mature_fa = {}
        novel_precursor_fa = {}
        mirna_info = self.work_dir + "/mirdp2_result/novel_miR_mature_infor.xls"
        mirna_detail = os.path.join(self.work_dir, "novel_mirna_detail.xls")
        for seq_record in SeqIO.parse(self.work_dir + "/mirdp2_result/novel_miR_mature.fa", "fasta"):
            novel_mature_fa[seq_record.id] = seq_record.seq
        for seq_record in SeqIO.parse(self.work_dir + "/mirdp2_result/novel_miR_pre.fa", "fasta"):
            novel_precursor_fa[seq_record.id] = seq_record.seq
        with open(mirna_info, "r") as f1, open(mirna_detail, "w") as w1:
            head_list = ["miRNA_id", "miRNA_seq", "miRNA_len", "pre_name", "pre_seq", "pre_len", "pre_position",
                         "energy"]
            w1.write("\t".join(head_list) + "\n")
            headline = f1.readline()
            for line in f1:
                items = line.strip().split("\t")
                mat_id = items[0]
                pre_id = items[5]
                pos = items[6] + "(" + items[9] + "):" + items[7] + "-" + items[8]
                energy = items[11]
                mature_seq = novel_mature_fa[mat_id]
                precursor_seq = novel_precursor_fa[pre_id]
                w1.write("\t".join(
                    [mat_id, str(mature_seq), str(len(mature_seq)), pre_id, str(precursor_seq), str(len(precursor_seq)),
                     pos, energy]) + "\n")

        prediction_file = glob.glob(self.work_dir + "/mirdp2/*/*filter_P_prediction")[0]
        filter_fa = os.path.join(self.output_dir, "filtered.fa")
        input_fa = self.option("input_fa").prop["path"]
        a = dict()
        with open(prediction_file, "r") as f1:
            for line in f1:
                reads = line.strip().split("\t")[11].split(";")
                for read in reads:
                    if read not in a:
                        a[read] = 1
        with open(filter_fa, "w") as w:
            for seq_id, seq_sequence in parse_fasta(input_fa):
                if seq_id not in a:
                    w.write('>%s\n%s\n' % (seq_id, seq_sequence))
        self.option("filter_fa", filter_fa)
        if os.path.exists(os.path.join(self.output_dir, "novel_precursors.str")):
            os.remove(os.path.join(self.output_dir, "novel_precursors.str"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_pre.tr", os.path.join(self.output_dir, "novel_precursors.str"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mature_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_mature_seq.fa"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_mature.fa", os.path.join(self.output_dir, "novel_mature_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursor_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_pre.fa", os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_count.xls"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_count.xls", os.path.join(self.output_dir, "novel_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_norm.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        os.link(self.work_dir + "/mirdp2_result/novel_miR_norm.xls", os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_detail.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        os.link(mirna_detail, os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        self.option("novel_count", os.path.join(self.output_dir, "novel_mirna_count.xls"))

    def rna_fold(self):
        self.logger.info("开始进行二级结构预测")
        pdfs = self.work_dir + "/mireap/pdfs"
        if os.path.exists(pdfs):
            shutil.rmtree(pdfs)
        else:
            os.mkdir(pdfs)
        cmd = "{}RNAfold < {} > {}".format(self.rnafold, os.path.join(self.work_dir, "mireap/novel_miR_pre.fa"),
                                           os.path.join(pdfs, "NovmiR.str"))
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('rnafold计算完成')
        except subprocess.CalledProcessError:
            self.logger.info('rnafold计算失败')
            self.set_error("rnafold运行失败")

    def ps_2_pdf(self):
        ps = glob.glob(self.work_dir + "/*.ps")
        for file in ps:
            file1 = file + ".pdf"
            cmd1 = "ps2pdf {} {}".format(file, file1)
            try:
                subprocess.check_output(cmd1, shell=True)
                self.logger.info('ps2pdf转换完成')
            except:
                self.set_error("ps2pdf转换失败")

    def set_mireap_results(self):
        self.logger.info("开始设置mireap结果目录")
        ## novel miRNA二级结构图片
        if os.path.exists(os.path.join(self.output_dir, "structure_pdf")):
            shutil.rmtree(os.path.join(self.output_dir, "structure_pdf"))
        os.mkdir(os.path.join(self.output_dir, "structure_pdf"))
        pdf = glob.glob(self.work_dir + "/*.pdf")
        for file in pdf:
            os.link(file, os.path.join(self.output_dir, "structure_pdf", os.path.basename(file)))

        novel_mature_fa = {}
        novel_precursor_fa = {}
        mirna_info = self.work_dir + "/mireap/novel_miR_mature_infor.xls"
        mirna_detail = os.path.join(self.work_dir, "novel_mirna_detail.xls")
        for seq_record in SeqIO.parse(self.work_dir + "/mireap/novel_miR_mature.fa", "fasta"):
            novel_mature_fa[seq_record.id] = seq_record.seq
        for seq_record in SeqIO.parse(self.work_dir + "/mireap/novel_miR_pre.fa", "fasta"):
            novel_precursor_fa[seq_record.id] = seq_record.seq
        with open(mirna_info, "r") as f1, open(mirna_detail, "w") as w1:
            head_list = ["miRNA_id", "miRNA_seq", "miRNA_len", "pre_name", "pre_seq", "pre_len", "pre_position",
                         "energy"]
            w1.write("\t".join(head_list) + "\n")
            headline = f1.readline()
            for line in f1:
                items = line.strip().split("\t")
                mat_id = items[0]
                pre_id = items[5]
                pos = items[6] + "(" + items[9] + "):" + items[7] + "-" + items[8]
                energy = items[11]
                mature_seq = novel_mature_fa[mat_id]
                precursor_seq = novel_precursor_fa[pre_id]
                w1.write("\t".join(
                    [mat_id, str(mature_seq), str(len(mature_seq)), pre_id, str(precursor_seq), str(len(precursor_seq)),
                     pos, energy]) + "\n")

        mrd = os.path.join(self.work_dir, "mireap-Nov.aln")
        filter_fa = os.path.join(self.output_dir, "filtered.fa")
        input_fa = self.option("input_fa").prop["path"]
        a = dict()
        with open(mrd, "r") as f1:
            for line in f1:
                if len(line.strip().split(" ")) >= 3:
                    if re.compile(r'S.*_\d+\s+').match(line.strip().split(" ")[1]):
                        mirna_id = line.strip().split(" ")[1] + "_" + line.strip().split(" ")[2]
                        a[mirna_id] = 1
        with open(filter_fa, "w") as w:
            for seq_id, seq_sequence in parse_fasta(input_fa):
                if seq_id not in a:
                    w.write('>%s\n%s\n' % (seq_id, seq_sequence))
        self.option("filter_fa", filter_fa)
        if os.path.exists(os.path.join(self.output_dir, "novel_precursors.str")):
            os.remove(os.path.join(self.output_dir, "novel_precursors.str"))
        os.link(self.work_dir + "/mireap/novel_miR_pre.tr", os.path.join(self.output_dir, "novel_precursors.str"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mature_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_mature_seq.fa"))
        os.link(self.work_dir + "/mireap/novel_miR_mature.fa", os.path.join(self.output_dir, "novel_mature_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_precursor_seq.fa")):
            os.remove(os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        os.link(self.work_dir + "/mireap/novel_miR_pre.fa", os.path.join(self.output_dir, "novel_precursor_seq.fa"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_count.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_count.xls"))
        os.link(self.work_dir + "/mireap/novel_miR_count.xls", os.path.join(self.output_dir, "novel_mirna_count.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_norm.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        os.link(self.work_dir + "/mireap/novel_miR_norm.xls", os.path.join(self.output_dir, "novel_mirna_norm.xls"))
        if os.path.exists(os.path.join(self.output_dir, "novel_mirna_detail.xls")):
            os.remove(os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        os.link(mirna_detail, os.path.join(self.output_dir, "novel_mirna_detail.xls"))
        self.option("novel_count", os.path.join(self.output_dir, "novel_mirna_count.xls"))

    def parse_config(self, file=None, section=None, name=None):
        config = ConfigParser.ConfigParser()
        config.read(file)
        return config.get(section, name)


class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "NovelMirna_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.novel_mirna_new",
            "instant": False,
            "options": dict(
                # species="sly",
                # category="Animal",
                # input_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/Ensemble_release_36/dna/Solanum_lycopersicum.SL2.50.dna_sm.toplevel.fa",
                # input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/filtered.fa",
                # config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/qc_file.config",
                # mapping_arf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/reads_vs_genome.arf",
                # mrd="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/miRBase.mrd",
                # known_miRNA_express="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/miRNAs_expressed_all_samples_total.csv",
                # index="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/Ensemble_release_36/dna/Solanum_lycopersicum.SL2.50.dna_sm.toplevel_index",
                # species="sly",
                # category="Plant",
                # software="mireap",
                # input_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/Ensemble_release_36/dna/Solanum_lycopersicum.SL2.50.dna_sm.toplevel.fa",
                # input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/filtered.fa",
                # clean_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/uniq.fasta",
                # config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/qc_file.config",
                # mapping_arf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/reads_vs_genome.arf",
                # mrd="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/miRBase.mrd",
                # known_miRNA_express="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/miRNAs_expressed_all_samples_total.csv",
                # index="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/Ensemble_release_36/dna/Solanum_lycopersicum.SL2.50.dna_sm.toplevel_index",
                species="ath",
                category="Plant",
                software="mirdp2",
                input_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/TAIR10_ensembl_44/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
                input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/filtered.fa",
                clean_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/uniq.fasta",
                config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/qc_file.config",
                mapping_arf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/reads_vs_genome.arf",
                mrd="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/miRBase.mrd",
                known_miRNA_express="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/miRNAs_expressed_all_samples_total.csv",
                index="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/TAIR10_ensembl_44/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel",
                rpm=0
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
