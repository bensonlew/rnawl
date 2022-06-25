#!/usr/bin/python
# -*- coding: utf-8 -*-
#!/usr/bin/python
# -*- coding: utf-8 -*-
# __author__ : konghualei
# last modify: 20170301

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.packages.ref_rna.express.single_sample import *
from mbio.files.sequence.fastq import FastqFile
import shutil
import os
import re
import subprocess

class RsemAgent(Agent):
    """调用Rsem1软件计算基因表达量(对转录本计算表达量), bowtie2方法建索引、比对"""
    def __init__(self, parent):
        super(RsemAgent, self).__init__(parent)
        options = [
            {"name": "fq_type", "type": "string"}, #PE OR SE
            {"name": "ref_denovo","type":"string", "default": "ref"}, #ref or denovo，判断是有参还是无参
            {"name": "transcript_fa", "type": "infile", "format": "sequence.fasta"},  # 转录本fasta文件
            # {"name": "Rsem1_fa", "type": "infile", "format": "sequence.fasta"},  # trinit.fasta文件
            {"name": "fa_build", "type": "outfile", "format": "sequence.fasta"},  # 转录本fasta构建好索引文件
            {"name": "fq_l", "type": "string", "default": "sequence.fastq"},  # PE测序，包含所有样本的左端fq文件的文件夹  不压缩的fq文件
            {"name": "fq_r", "type": "string", "default": "sequence.fastq"},  # PE测序，包含所有样本的左端fq文件的文件夹
            {"name": "fq_s", "type": "string", "default": "sequence.fastq"},  # SE测序，包含所有样本的fq文件的文件夹
            # {"name": "bam", "type": "infile", "format": "align.bwa.bam, ref_rna.assembly.bam_dir"}, # bam文件输入
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf" }, # gtf/gff文件(目前是merged.gtf文件-提取gene2transcript map文件)
            #{"name": "gtf_type", "type": "string", "default": "ref"}, # "ref" "merged_stringtie" "merged_cufflinks"
            {"name": "sample_name", "type":"string"}, # 样本名称
            {"name": "cpu", "type": "int", "default": 8},  # 设置CPU
            {"name": "max_memory", "type": "string", "default": "100G"}, # 设置内存
            {"name": "only_bowtie_build", "type": "bool", "default": False},  # 为true时该tool只建索引
            {"name": "bowtie_build_Rsem1", "type": "bool", "default": False},  # 为true时该tool需要建索引
            {"name": "strand_specific", "type": "bool", "default": False},  # added by gdq
            {"name": "strand_dir", "type": "string", "default": "forward"},  # added by gdq
        ]
        self.add_option(options)
        self._memory_increase_step = 20
        self.step.add_steps("Rsem1")
        self.on("start", self.stepstart)
        self.on("end", self.stepfinish)

    def stepstart(self):
        self.step.Rsem1.start()
        self.step.update()

    def stepfinish(self):
        self.step.Rsem1.finish()
        self.step.update()

    def check_options(self):
        """重写参数检测函数"""
        self.logger.info("输出ref_gtf格式")
        self.logger.info(self.option('ref_gtf').format)
        if not self.option("transcript_fa").is_set:
            raise OptionError("请输入参考基因组的fa文件")
        if self.option("fq_type") not in ['PE', 'SE']:
            raise OptionError("请输入单端或双端!")
        if not self.option("ref_gtf").is_set:
            raise OptionError("请输入参考基因组gtf文件!")

    def set_resource(self):
        self._cpu = 8
        self._memory = '40G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(RsemAgent, self).end()

class RsemTool(Tool):
    def __init__(self, config):
        super(RsemTool, self).__init__(config)
        # self.version = "1.0.1"
        #self.Rsem1_path =self.config.SOFTWARE_DIR + '/bioinfo/rna/Rsem1-1.2.31/bin/'
        self.Rsem1 =  "/bioinfo/rna/scripts/align_and_estimate_abundance.pl"
        self.Rsem1_path = 'bioinfo/rna/RSEM-1.2.31/bin/'
        self.bowtie_path = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie2-2.2.9/'
        #self.star_build_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/star-2.5/bin/Linux_x86_64/"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.perl = self.config.SOFTWARE_DIR + '/miniconda2/bin'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.set_environ(PATH=self.Rsem1_path)
        self.set_environ(PATH=self.bowtie_path)
        self.set_environ(PATH=self.perl)

    def ref_Rsem1_build(self):
        self.Rsem1_index_path = os.path.join(self.work_dir, "Rsem1_index")
        gtf_format = os.path.basename(self.option('ref_gtf').prop['path'])
        gene2transcript = os.path.join(self.work_dir, 'gene2transcript')
        if re.search(r"gtf$",gtf_format):
            gtf(self.option('ref_gtf').prop['path'], gene2transcript)
            self.logger.info("对gtf文件提取gene_transcript成功！")
        elif re.search(r'gff$|gff3$', gtf_format):
            gff3(infile = self.option('ref_gtf').prop['path'], outfile = gene2transcript, out_path=self.work_dir)
            self.logger.info("对gtf文件提取gene_transcript成功！")
        if os.path.exists(gene2transcript):
            self.logger.info("建立gene_id和transcript_id对应关系成功！")
        #if self.option("transcript_fa").is_set:
        cmd = self.Rsem1_path + "rsem-prepare-reference -p 8 --transcript-to-gene-map {} {} {} --bowtie2 --bowtie2-path {} ".format(gene2transcript, self.option("transcript_fa").prop['path'],self.Rsem1_index_path, self.bowtie_path)
        self.logger.info(cmd)
        bowtie2_cmd = self.add_command("bowtie2_build", cmd, ignore_error=True).run()
        self.wait(bowtie2_cmd)
        if bowtie2_cmd.return_code == 0:
            self.logger.info("%s运行完成" % bowtie2_cmd)
            #self.option("fa_build").set_path(self.Rsem1_index_path)
        elif bowtie2_cmd.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("%s运行出错!" % bowtie2_cmd)

    def ref_Rsem1_run(self, fasta_build):
        data = FastqFile()
        if self.option("fq_type") == "PE":
            if self.option('strand_specific'):
                if self.option('strand_dir') == "forward":
                    forward_prob = 0
                else:
                    forward_prob = 1
            else:
                forward_prob = 0.5
            cmd = self.Rsem1_path + "rsem-calculate-expression " \
                                    "--paired-end " \
                                    "-p 8 " \
                                    "{} {} {} {} " \
                                    "--bowtie2 " \
                                    "--bowtie2-path {} " \
                                    "--forward-prob {}".format(self.option("fq_l"),
                                                               self.option("fq_r"),
                                                               fasta_build,
                                                               self.output_dir +"/"+self.option("sample_name"),
                                                               self.bowtie_path,
                                                               forward_prob)
        elif self.option("fq_type") == "SE":
            cmd = self.Rsem1_path + "rsem-calculate-expression -p 8 {} {} {} --bowtie2 --bowtie2-path {}".format(self.option("fq_s"),\
                   fasta_build, self.output_dir+"/"+self.option("sample_name"), self.bowtie_path)
        self.logger.info(cmd)
        exp_cmd = self.add_command("express_run", cmd, ignore_error=True).run()
        self.wait(exp_cmd)
        if exp_cmd.return_code == 0:
            self.logger.info("{}运行成功！".format(exp_cmd))
        elif exp_cmd.return_code in [1, -9, 255]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("{}运行失败".format(exp_cmd))

    def denovo_bowtie_build(self):
        Rsem1_fasta = self.work_dir + '/' + os.path.basename(self.option('transcript_fa').prop['path'])
        if os.path.exists(Rsem1_fasta):
            os.remove(Rsem1_fasta)
            os.link(self.option('transcript_fa').prop['path'], Rsem1_fasta)
        else:
            os.link(self.option('transcript_fa').prop['path'], Rsem1_fasta)
        cmd = self.Rsem1 + ' --transcripts %s --seqType fq --single test.fq --est_method  Rsem1 --output_dir %s --trinity_mode --aln_method bowtie2 --prep_reference' % (Rsem1_fasta, self.work_dir)
        self.logger.info('开始运行bowtie2建索引')
        bowtie_cmd = self.add_command('bowtie_build', cmd, ignore_error=True).run()
        self.wait(bowtie_cmd)
        if bowtie_cmd.return_code == 0:
            self.logger.info("%s运行完成" % bowtie_cmd)
            self.option('fa_build', Rsem1_fasta)
        elif bowtie_cmd.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("%s运行出错!" % bowtie_cmd)

    def denovo_run_Rsem1(self, Rsem1_fasta):
        if self.option('fq_type') == 'SE':
            sample = os.path.basename(self.option('fq_s')).split('_sickle_s.fastq')[0]
            Rsem1_cmd = self.Rsem1 + ' --transcripts %s --seqType fq --single %s --est_method  Rsem1 --output_dir %s --thread_count 6 --trinity_mode --aln_method bowtie2 --output_prefix %s' % (Rsem1_fasta, self.option('fq_s'), self.work_dir, sample)
        else:
            sample = os.path.basename(self.option('fq_l')).split('_sickle_l.fastq')[0]
            Rsem1_cmd = self.Rsem1 + ' --transcripts %s --seqType fq --right %s --left %s --est_method  Rsem1 --output_dir %s --thread_count 6 --trinity_mode --aln_method bowtie2 --output_prefix %s' % (Rsem1_fasta, self.option('fq_r'), self.option('fq_l'), self.work_dir, sample)

        self.logger.info("开始运行_Rsem1_cmd")
        cmd = self.add_command("Rsem1_cmd", Rsem1_cmd, ignore_error=True).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("%s运行完成" % cmd)
        elif cmd.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("%s运行出错!" % cmd)

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        #for root, dirs, files in os.walk(self.output_dir):
        #    for names in files:
        #        os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir)
        try:
            for f in results:
                if re.search(r'results$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                    # os.system(""" sed -i "s/gene://g" %s """ % (file_path))
                    # os.system(""" sed -i "s/transcript://g" %s """ % (file_path))
            self.logger.info("设置Rsem1分析结果目录成功")
        except Exception as e:
            self.logger.info("设置Rsem1分析结果目录失败{}".format(e))

    def run(self):
        super(RsemTool, self).run()
        if self.option("only_bowtie_build") == True:
            if self.option("ref_denovo") == "ref":
                self.ref_Rsem1_build()
            if self.option("ref_denovo") == "denovo":
                self.denovo_bowtie_build()
            self.logger.info("Rsem1建立索引成功！")
        else:
            if self.option("bowtie_build_Rsem1") == False:
                if self.option("ref_denovo") == "ref":
                    self.ref_Rsem1_build()
                    # self.ref_Rsem1_run(self.option("fa_build").prop['path'])
                    self.ref_Rsem1_run(self.Rsem1_index_path)
                if self.option("ref_denovo") == "denovo":
                    self.denovo_bowtie_build()
                    self.denovo_run_Rsem1(self.option('fa_build').prop['path'])
                self.logger.info("Rsem1计算表达量成功！")
            else:
                """输入建好索引的文件，直接进行表达量的计算"""
                if self.option("ref_denovo") == "ref":
                    self.ref_Rsem1_run(self.option("transcript_fa").prop["path"])
                if self.option("ref_denovo") == "denovo":
                    self.denovo_run_Rsem1(self.option('transcript_fa').prop['path'])
        self.set_output()
        self.end()
