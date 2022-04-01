#!/usr/bin/python
# -*- coding:utf-8 -*-
# __author__:khl
# testkhl
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.packages.ref_rna.express.single_sample import *
from mbio.packages.ref_rna.express.set_strand import set_strand
from mbio.packages.ref_rna.express.express_distribution import *
import shutil
import os
import re
from mbio.files.meta.otu.group_table import *


class FeaturecountsAgent(Agent):

    def __init__(self, parent):
        super(FeaturecountsAgent, self).__init__(parent)
        options = [
            {"name": "fq_type", "type": "string","default": "PE"},  # PE OR SE
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因组的gtf文件
            {"name": "bam", "type": "infile", "format": "align.bwa.bam,align.bwa.bam_dir"},  # 样本比对后的bam文件或文件夹
            {"name": "strand_specific", "type": "bool", "default": False},  # PE测序，是否链特异性, 默认是0, 无特异性
            {"name": "strand_dir", "type": "string","default": "None"},  # "forward", "reverse" 默认不设置此参数
            {"name": "count", "type":"outfile", "format": "rna.express_matrix"},  # featurecounts 输出结果总结
            {"name": "fpkm", "type": "outfile", "format": "rna.express_matrix"},  # fpkm表达量表
            {"name": "tpm", "type": "outfile", "format": "rna.express_matrix"},  # tpm表达量表
            {"name": "out_file", "type": "outfile", "format": "rna.express_matrix"}, # featureCounts软件直接生成的表达量文件
            {"name": "gtf", "type": "string"},  # 该物种的参考基因组文件路径  gtf格式而非gff格式
            {"name": "cpu", "type": "int", "default": 10},  # 设置CPU
            {"name": "is_duplicate", "type": "bool"},  # 是否有生物学重复
            {"name": "edger_group", "type":"infile", "format":"sample.group_table"},
            {"name": "max_memory", "type": "string", "default": "100G"},  #设置内存
            {"name": "exp_way", "type": "string","default": "fpkm"},  # fpkm水平表达量  fpkm tpm all
            {"name": "all_gene_list", "type": "outfile", "format": "rna.express_matrix"},  #所有基因list列表，提供给下游差异分析module
        ]
        self.add_option(options)
        self.step.add_steps("featurecounts")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.featurecounts.start()
        self.step.update()

    def stepfinish(self):
        self.step.featurecounts.finish()
        self.step.update()

    def check_options(self):
        if not self.option('fq_type'):
            raise OptionError('必须设置测序类型：PE OR SE')
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内')
        if not self.option("strand_specific"):
            if not self.option("strand_dir"):
                raise OptionError("链特异性时需要选择正链或者负链")
        if not self.option("ref_gtf").is_set:
            raise OptionError("需要输入gtf文件")
        #if self.option("feature_id") not in ["gene_id","exon","both"]:
        #    raise OptionError("计算基因或者外显子的count值")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '30G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "summary", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            #[r"featurecounts+\.summary", "summary", "summary 记录"],
            #[r"featurecounts+\^sample\d", "","gene count 表"]
        ])
        super(FeaturecountsAgent, self).end()


class FeaturecountsTool(Tool):

    def __init__(self, config):
        super(FeaturecountsTool, self).__init__(config)
        # self._version = '1.0.1'
        self.featurecounts_path = 'bioinfo/align/subread-1.5.0/bin/featureCounts'
        self.set_environ(PATH = self.featurecounts_path)
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '
        self.parse_featurecounts_perl = os.path.join(Config().SOFTWARE_DIR, 'bioinfo/rna/scripts/parse_featurecount_all.pl ')
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:$PATH"
        self._r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path1 = "/program/R-3.3.1/bin/Rscript "
        self.python_path = self.config.SOFTWARE_DIR+"/program/Python/bin:$PATH"
        self.set_environ(PATH=self.python_path)
        self.distribution_path = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/denovo_rna/express'

    def featurecounts_run(self):
        self.logger.info("开始进行表达量分析！")
        self.new_gtf = self.option('ref_gtf').prop['path']
        self.logger.info(self.new_gtf)
        if os.path.isdir(self.option('bam').prop['path']):
            file_name=[]
            self.samples=[]
            for files in os.listdir(self.option('bam').prop['path']):
                file_name.append(os.path.join(self.option('bam').prop['path'],files))
                self.samples.append(files.split('.bam')[0])
            self.input_bam_file=" ".join(file_name)
            self.logger.info(self.samples)
            self.output_name="-".join(self.samples)
            # self._out_dir= os.path.join(self.output_dir, "-".join(self.samples))
            self._out_dir= os.path.join(self.output_dir, "featureCounts_result")
        else:
            self.output_name=os.path.basename(self.option("bam").prop['path']).split('.bam')[0]
            self._out_dir=os.path.join(self.output_dir,"temp_sample")
            self.input_bam_file=self.option('bam').prop['path']
        if self.option('fq_type') == 'PE':
            self.logger.info("双端测序")
            if self.option("strand_specific"):
                self.strand_specific, self.strand_dir = set_strand(method='featurecounts', strand_specific = True, strand_dir = self.option('strand_dir'))
                cmd = self.featurecounts_path + " -T %s -p -a %s -t exon -M -s %s -o %s %s" % (self.option('cpu'), self.new_gtf,
                     self.strand_dir, self._out_dir, self.input_bam_file)
            else:
                    cmd = self.featurecounts_path + " -T %s -p -a %s -t exon -M -s %s -o %s %s" % (self.option('cpu'), self.new_gtf,
                         0, self._out_dir, self.input_bam_file)
            self.logger.info(cmd)
        else:
            cmd = self.featurecounts_path + " -T %s -a %s -t exon -M -o %s %s" % (self.option('cpu'), self.new_gtf,
                         self._out_dir, self.input_bam_file)
        self.logger.info("开始运行featureCounts计算表达量")
        featurecounts_cmd = self.add_command("featurecounts", cmd).run()
        self.wait(featurecounts_cmd)
        if featurecounts_cmd.return_code == 0:
            self.logger.info("%s运行完成" % featurecounts_cmd)
        else:
            self.set_error("%s运行出错" % cmd)

    def fpkm_tpm(self):
        """计算fpkm，tpm的表达量"""
        if os.path.exists(self._out_dir):
            # remove_header(file_path=self._out_dir, file_name="vs".join(os.path.basename(self._out_dir).split("-")))
            remove_header(file_path=self._out_dir, file_name="featureCounts_remove_header")
        else:
            remove_header(file_path=self._out_dir, file_name=self.output_name)
        self.count_path = os.path.join(self.output_dir,"count.xls")
        self.gene_length_path = os.path.join(self.output_dir,"gene_length.xls")
        self.logger.info(self._out_dir)
        self.logger.info(self.count_path)
        self.logger.info(self.gene_length_path)
        # prepare(input_file=os.path.join(self.output_dir, "vs".join(os.path.basename(self._out_dir).split("-"))), gtf_file=self.new_gtf, count_matrix=self.count_path, gene_length=self.gene_length_path)
        prepare(input_file=os.path.join(self.output_dir, "featureCounts_remove_header"),
                gtf_file=self.new_gtf, count_matrix=self.count_path, gene_length=self.gene_length_path)
        fpkm_tpm_cmd= self.perl_path+self.parse_featurecounts_perl+self.count_path +" " +self.gene_length_path + " "+"fpkm_tpm"
        _fpkm_tpm_cmd = self.add_command("fpkm.tpm", fpkm_tpm_cmd).run()
        self.wait(_fpkm_tpm_cmd)
        if _fpkm_tpm_cmd.return_code ==0:
            self.logger.info("计算fpkm, tpm表达量成功!\n{}".format(_fpkm_tpm_cmd))
        else:
            self.set_error("计算fpkkm, tpm表达量失败!\n{}".format(_fpkm_tpm_cmd))

    def get_distribution(self, old_fpkm, out_path, exp_way=None):
        if not os.path.exists(out_path):
            os.mkdir(out_path)
            self.logger.info("文件夹{}生成".format(out_path))
        distribution(rfile=out_path+"/express_distribution.r",input_matrix=old_fpkm,outputfile=out_path,filename="gene")
        gcmd=self.r_path1+out_path+"/express_distribution.r"
        cmd1 = self.add_command("gene_{}_cmd".format(exp_way), gcmd).run()
        self.wait(cmd1)
        if cmd1.return_code == 0:
            self.logger.info("表达量分布图{}的数据分析成功".format(exp_way))
        else:
            self.set_error("表达量分布图{}的数据分析出错".format(exp_way))

    def sample_distribution(self):
        #样本count  distribution
        self.get_distribution(self.output_dir + "/count.xls", self.work_dir + "/count", "count")
        if self.option("exp_way") == "fpkm":
            self.get_distribution(self.output_dir+"/fpkm_tpm.fpkm.xls",self.work_dir+"/fpkm","fpkm")
            self.logger.info("计算fpkm成功")
        if self.option("exp_way") == "tpm":
            self.get_distribution(self.output_dir+"/fpkm_tpm.tpm.xls",self.work_dir+"/tpm","tpm")
            self.logger.info("计算tpm成功")
        if self.option("exp_way") == "all":
            self.get_distribution(self.output_dir+"/fpkm_tpm.fpkm.xls",self.work_dir+"/fpkm","fpkm")
            self.logger.info("计算fpkm成功")
            self.get_distribution(self.output_dir+"/fpkm_tpm.tpm.xls",self.work_dir+"/tpm","tpm")
            self.logger.info("计算tpm成功")

    def group_distribution(self, old_fpkm, new_fpkm,old_count, new_count, sample_group_info, outputfile, filename,_type=None): # add by khl
        """
        :param _type: fpkm或者tpm
        """
        group_express(old_fpkm = old_fpkm, new_fpkm = new_fpkm, old_count = old_count, new_count = new_count, sample_group_info = sample_group_info, \
                        outputfile = outputfile, filename = filename)
        # cmd = self.r_path1 + " {}".format(rfile)
        # cmd1 = self.add_command("{}".format(filename+_type), cmd).run()
        # self.logger.info(cmd)
        # self.wait()

        group_fpkm_distribution = self.work_dir + "/group_{}_distribution".format(_type)
        if not os.path.exists(group_fpkm_distribution):
            os.mkdir(group_fpkm_distribution)

        distribution(rfile='./group_{}_distribution.r'.format(_type), input_matrix=new_fpkm, outputfile=group_fpkm_distribution,
                     filename=filename)
        cmd = self.r_path1 + " group_{}_distribution.r".format(_type)
        cmd_run = self.add_command("{}".format((filename + "{}".format(_type)).lower()), cmd).run()

        group_count_distribution = self.work_dir + "/group_count_distribution"
        if not os.path.exists(group_count_distribution):
            os.mkdir(group_count_distribution)
            distribution(rfile='./group_count_distribution.r', input_matrix=new_count, outputfile=group_count_distribution,
                     filename=filename)
            cmd1 = self.r_path1 + " group_count_distribution.r"
            cmd1_run = self.add_command("{}".format((filename + "count").lower()), cmd1).run()
        self.wait(cmd_run)
        if cmd_run.return_code == 0:
            self.logger.info("计算group{}密度分布成功".format(_type))
        else:
            self.logger.info("计算group{}密度分布失败".format(_type))

    def group_detail(self):
        g = GroupTableFile()
        if self.option("edger_group").is_set:
            g.set_path(self.option("edger_group").prop['path'])
            sample_group_info = g.get_group_spname()
            self.logger.info("打印sample_group_info信息")
            self.logger.info(sample_group_info)
            outputfile = self.output_dir
            results=os.listdir(outputfile)
            new_group_path = self.work_dir + "/group"
            if not os.path.exists(new_group_path):
                os.mkdir(new_group_path)
            group_fpkm = new_group_path+"/fpkm"
            group_tpm = new_group_path + "/tpm"
            if not os.path.exists(group_fpkm):
                os.mkdir(group_fpkm)
            if not os.path.exists(group_tpm):
                os.mkdir(group_tpm)
            old_count = outputfile+"/count.xls"

            self.group_distribution(old_fpkm=outputfile+"/fpkm_tpm.fpkm.xls", new_fpkm=new_group_path+"/group.fpkm.xls",
                                    old_count=old_count, new_count=outputfile+"/fpkm_count.xls",
                            sample_group_info=sample_group_info, outputfile=group_fpkm,filename='GroupGenes', _type='fpkm')
            self.group_distribution(old_fpkm=outputfile+"/fpkm_tpm.tpm.xls",new_fpkm=new_group_path+"/group.tpm.xls",
                                    old_count = old_count, new_count = outputfile+"/tpm_count.xls",
                            sample_group_info=sample_group_info, outputfile=group_tpm,filename='GroupGenes', _type='tpm')

            self.logger.info("计算基因转录本group成功！")
        else:
            raise Exception("有生物学重复时，请设置样本生物学分组信息！")

    def set_output(self):
        self.logger.info("设置结果目录")
        if os.path.exists(self.count_path):
            self.option("count").set_path(self.count_path)
        else:
            raise Exception("没有生成表达量表")
        for files in os.listdir(self.work_dir):
            if files.endswith("fpkm.xls"):
                if os.path.exists(os.path.join(self.output_dir,files)):
                    os.remove(os.path.join(self.output_dir,files))
                os.link(os.path.join(self.work_dir,files), os.path.join(self.output_dir,files))
                self.option("fpkm").set_path(os.path.join(self.output_dir, files))
            if files.endswith("tpm.xls"):
                if os.path.exists(os.path.join(self.output_dir,files)):
                    os.remove(os.path.join(self.output_dir,files))
                os.link(os.path.join(self.work_dir,files), os.path.join(self.output_dir,files))
                self.option("tpm").set_path(os.path.join(self.output_dir, files))
        self.option("out_file").set_path(self._out_dir)  #featurecounts软件生成文件
        # all_gene_list(file_path=os.path.join(self.output_dir, "vs".join(os.path.basename(self._out_dir).split("-"))), all_gene_list_path = os.path.join(self.output_dir, "all_gene_list"))
        all_gene_list(file_path=os.path.join(self.output_dir, "featureCounts_remove_header"),
                      all_gene_list_path=os.path.join(self.output_dir, "all_gene_list"))
        self.option("all_gene_list").set_path(os.path.join(self.output_dir, "all_gene_list"))
        self.logger.info("生成gene_list成功！")
        self.logger.info("设置count表达量成功")

    def run(self):
        super(FeaturecountsTool, self).run()
        self.featurecounts_run()
        self.fpkm_tpm()
        self.set_output()
        self.sample_distribution()
        if self.option("is_duplicate"):
            self.group_detail()
        self.end()


def remove_header(file_path, file_name):
    """除去生成文件的header标签，同时将最后一列列名更新为样本名"""
    if os.path.exists(file_path):
        file_path_path = os.path.split(file_path)[0]
        # column_new_name = [ss.split('.bam')[0]for ss in os.path.basename(file_path).split('-')]
        new_file_path = os.path.join(file_path_path, "count")
        os.system('grep \"#\" -v {} > {}'.format(file_path, new_file_path))
        output_file_path = os.path.join(file_path_path, file_name)
        file1 = open(output_file_path, 'w+')
        with open(new_file_path, 'r+') as files:
            data1 = files.readline().strip().split('\t')
            data = data1[:6]
            column_name = data1[6:]
            _column_name = [os.path.basename(ss).split('.bam')[0]for ss in column_name]
            file1.write('\t'.join(data)+'\t'+'\t'.join(_column_name)+'\n')
            i = 0
            for f in files:
                i += 1
                line = f.strip().split("\t")
                for f1 in range(len(line)):
                    if f1 != len(line) - 1:
                        file1.write(line[f1] + "\t")
                    else:
                        #if i == 1:
                        #    file1.write("count" +"\n")
                        #else:
                        file1.write(line[f1] + "\n")
        os.remove(new_file_path)
        file1.close()
        if os.path.exists(output_file_path):
             pass
    else:
        raise Exception("输入文件不存在，无法除掉文件标签！")
