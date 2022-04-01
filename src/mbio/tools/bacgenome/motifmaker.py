# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os, shutil
import re
import pandas as pd

class MotifmakerAgent(Agent):
    """
    三代数据甲基化位点预测
    用的是motifmaker软件识别位点
    bam文件
    """

    def __init__(self, parent):
        super(MotifmakerAgent, self).__init__(parent)
        options = [
            {"name": "input", "type": "infile", 'format': "align.bwa.bam"}, # 输入bam文件
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "ref_database", "type": "infile", "format": "sequence.fasta"}, # 输入参考序列
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"},## 输出结果文件
            {"name": "seq_type", "type": "string", 'default': "bam"},  # 文件类型，pacbio一般为bam文件；nanopore一般为fastq文件
            {"name": "identify", "type": "string", 'default': "m6A,m4C"}, # Specific modifications to identify (comma-separated list). Currrent options are m6A, m4C, m5C_TET. Using --control overrides this option. (default: m6A,m4C)
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("input"):
            raise OptionError("必须设置输入bam文件")
        if not self.option("sample"):
            raise OptionError("必须设置输入样本名称")
        if not self.option("ref_database"):
            raise OptionError("必须设置输入参考序列文件")
        return True

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = 9
        self._memory = '120G'

    def end(self):
        super(MotifmakerAgent, self).end()


class MotifmakerTool(Tool):
    def __init__(self, config):
        super(MotifmakerTool, self).__init__(config)
        self._version = "1.0"
        self.pbalign = '/bioinfo/Genomic/Sofware/smrtlink_10.1.0.1/smrtlink/smrtcmds/bin/pbmm2'
        self.pbindex = '/bioinfo/Genomic/Sofware/smrtlink_10.1.0.1/smrtlink/smrtcmds/bin/pbindex'
        self.ipdsummary = "/bioinfo/Genomic/Sofware/smrtlink_10.1.0.1/smrtlink/smrtcmds/bin/ipdSummary"
        self.motifmaker = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/MultiMotifMaker-master/artifacts/MultiMotifMaker.jar"
        self.java = "/program/sun_jdk1.8.0/bin/java"
        self.path = self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/bin:"+ self.config.SOFTWARE_DIR + "/program/Python/bin:"
        self.set_environ(PATH=self.path)
        self.samtools = "/program/Python/bin/samtools"
        self.seq = self.option("input").prop["path"]
        self.ref = self.option("ref_database").prop["path"]
        self.align = self.work_dir + "/align.bam"
        self.gff = self.work_dir + "/ref.gff"
        self.csv = self.work_dir + "/ref.csv"
        self.perl_script = self.config.PACKAGE_DIR + '/bacgenome/get_methy.pl '
        self.perl_path = "/program/perl-5.24.0/bin/perl"

    def run(self):
        super(MotifmakerTool, self).run()
        self.run_mapping()
        self.run_index()
        self.run_index2()
        self.run_predict()
        self.run_motif()
        self.run_motif2()
        self.get_motif_detail()
        self.set_output()
        self.end()

    def run_mapping(self):
        """
        mapping数据库得到align结果
        :return:
        """
        # 先建参考的minimap2索引
        self.logger.info("开始建参考序列的minimap2的索引！")
        self.ref_index = self.work_dir + "/" + os.path.basename(self.ref).split(".")[0]+".mmi"
        cmd = '{} index {} {}'.format(self.pbalign, self.ref, self.ref_index)
        command_index = self.add_command('minimap2_index', cmd).run()
        self.wait(command_index)
        if command_index.return_code == 0:
            self.logger.info("minimap2_index成功!")
        else:
            self.set_error("minimap2_index失败!")

        self.logger.info("开始做mapping！")
        cmd = '{} align {} {} {} -j 9 '.format(self.pbalign, self.ref_index, self.seq,self.align)
        command = self.add_command('mapping', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("mapping成功!")
        else:
            self.set_error("mapping失败!")

    def run_index(self):
        """
        对bam文件建索引
        :return:
        """
        self.logger.info("开始做index！")
        cmd = '{} {} '.format(self.pbindex, self.align)
        command = self.add_command('index', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("建索引成功！")
        else:
            self.set_error("建索引失败！")

    def run_index2(self):
        """
        对fasta文件建索引
        :return:
        """
        self.logger.info("开始做index！")
        cmd = '{} faidx {} '.format(self.samtools, self.ref)
        command = self.add_command('fasta_index', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fasta建索引成功！")
        else:
            self.set_error("fasta建索引失败！")

    def run_predict(self):
        """
        得到注释的gff、csv等文件
        :return:
        """
        self.logger.info("开始做predict！")
        cmd = '{} {} --reference {} --gff {} --identify {} --numWorkers 8 --methylFraction'.format(self.ipdsummary,
                                                                                                   self.align,
                                                                                                   self.ref,
                                                                                                   self.gff,
                                                                                                   self.option("identify"))
        command = self.add_command('predict', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("甲基化预测成功")
        else:
            self.set_error("甲基化预测失败")

    def run_motif(self):
        """
        运行motifmaker软件
        :return:
        """
        self.motif_gff = self.work_dir + "/motifs.csv"
        if os.path.exists(self.motif_gff):
            os.remove(self.motif_gff)
        cmd = '{} -jar {} find -f {} -g {} -o {} -t 9'.format(self.java, self.motifmaker, self.ref, self.gff,
                                                              self.motif_gff)
        command = self.add_command('motifmaker', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("motifmaker生成成功")
        else:
            self.set_error("motifmaker生成失败")

    def run_motif2(self):
        """
        motifmaker软件运行之后没有motif信息，需要提取
        :return:
        """
        self.motif_gff2 = self.work_dir + "/motifs.gff"
        cmd = '{} -jar {} reprocess -f {} -g {} -m {} -o {}'.format(self.java, self.motifmaker, self.ref, self.gff,
                                                              self.motif_gff,self.motif_gff2)
        command = self.add_command('motifmaker_get_information', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("motifmaker提取信息成功")
        else:
            self.set_error("motifmaker提取信息失败")

    def get_motif_detail(self):
        """
        根据gff结果，整理成页面所需要的结果
        :return:
        """
        if os.path.exists(self.output_dir + "/" +self.option("sample")):
            shutil.rmtree(self.output_dir + "/" +self.option("sample"))
        os.mkdir(self.output_dir + "/" +self.option("sample"))
        motif_detail = self.output_dir + "/" + self.option("sample") + "/{}.motif_detail.xls".format(self.option(
            "sample"))
        #self.motif_gff2 = self.gff
        if not os.path.exists(self.motif_gff2):
            self.set_error("motif_gff文件生成失败".format(self.motif_gff2))
        cmd = '{} {} {} {} {} '.format(self.perl_path, self.perl_script, self.motif_gff2, self.option('sample'),motif_detail)
        command = self.add_command('get_motif_detail', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("get_motif_detail文件生成成功")
        else:
            self.set_error("get_motif_detail文件生成失败")

    def set_output(self):
        """
        设置结果文件
        :return:
        """
        self.logger.info("开始设置结果文件!")
        motif_csv = pd.read_csv(self.motif_gff)
        link_motif_csv = self.output_dir + "/" +self.option("sample") + "/{}.motifs.csv".format(self.option("sample"))
        link_motif_gff = self.output_dir + "/" +self.option("sample") + "/{}.motif.gff".format(self.option(
            "sample"))
        if len(motif_csv) > 1:
            if os.path.exists(link_motif_csv):
                os.remove(link_motif_csv)
            if os.path.exists(self.motif_gff):
                self.logger.info(self.motif_gff)
                shutil.copyfile(self.motif_gff, link_motif_csv)
            if os.path.exists(self.motif_gff2):
                self.logger.info(self.motif_gff2)
                shutil.copyfile(self.motif_gff2, link_motif_gff)
