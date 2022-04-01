#!/usr/bin/python
# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna.express.express_distribution import distribution
import os
import re
from mbio.packages.ref_rna.express.single_sample import group_express
from mbio.files.meta.otu.group_table import *
from mbio.packages.ref_rna.express.cmp_ref_cls import *
from mbio.packages.ref_rna.express.single_sample import *
import shutil


class MergeRsemAgent(Agent):
    """
    调用align_and_estimate_abundance.pl脚本，运行rsem，进行表达量计算分析
    version v1.0
    author: qiuping
    last_modify: 2016.06.20
    """

    def __init__(self, parent):
        super(MergeRsemAgent, self).__init__(parent)
        options = [
            {"name": "gtf_ref", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因组gtf文件
            {"name": "gtf_merged", "type": "infile", "format": "gene_structure.gtf"},
            # cat ref.gtf和new_transcript.gtf文件生成的新的gtf文件
            # {"name": "transcript_fa", "type": "infile", "format": "sequence.fasta"},  # 转录本的fa序列信息，用来提取最长的转录本序列(替代基因的序列)
            {"name": "is_class_code", "type": "bool"},  # 是否计算class_code信息(表达量RNA流程分析时设置此参数为False)
            {"name": "rsem_files", "type": "infile", "format": "rna.rsem_dir"},  # SE测序，包含所有样本的fq文件的文件夹
            {"name": "tran_count", "type": "outfile", "format": "rna.express_matrix"},
            {"name": "gene_count", "type": "outfile", "format": "rna.express_matrix"},
            {"name": "tran_fpkm", "type": "outfile", "format": "rna.express_matrix"},
            {"name": "gene_fpkm", "type": "outfile", "format": "rna.express_matrix"},
            {"name": "is_duplicate", "type": "bool"},  # 是否有生物学重复 
            {"name": "class_code", "type": "outfile", "format": "rna.express_matrix"},  # class_code 文件基本信息
            {"name": "edger_group", "type": "infile", "format": "sample.group_table"},
            {"name": "exp_way", "type": "string", "default": "fpkm"},
        ]
        self.add_option(options)
        self.step.add_steps("rsem")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.rsem.start()
        self.step.update()

    def stepfinish(self):
        self.step.rsem.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('rsem_files'):
            raise OptionError('必须设置输入文件：rsem结果文件')
        if self.option("exp_way") not in ['fpkm', 'tpm']:
            raise OptionError("所设表达量的代表指标不在范围内，请检查")
        if self.option("is_duplicate"):
            if not self.option("edger_group").is_set:
                raise Exception("有生物学重复时，请设置样本生物学分组信息！")
        if self.option("is_class_code"):
            if not os.path.exists(self.option("gtf_ref").prop['path']) and not os.path.exists(
                    self.option("gtf_merged").prop['path']):
                raise Exception("请设置gtf_ref和gtf_merged参数！")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [r"matrix$", "xls", "表达量矩阵"]
        ])
        super(MergeRsemAgent, self).end()


class MergeRsemTool(Tool):
    """
    Lefse tool
    """

    def __init__(self, config):
        super(MergeRsemTool, self).__init__(config)
        self._version = '1.0.1'
        self.fpkm = "/bioinfo/rna/scripts/abundance_estimates_to_matrix.pl"
        self.tpm = "/bioinfo/rna/trinityrnaseq-2.2.0/util/abundance_estimates_to_matrix.pl"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:$PATH"
        self._r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.long_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/"
        self.python_path = "program/Python/bin/"
        self.r_path1 = "/program/R-3.3.1/bin/Rscript"

    def merge_rsem(self):
        files = os.listdir(self.option('rsem_files').prop['path'])
        if self.option('exp_way') == 'fpkm':
            merge_gene_cmd = self.fpkm + ' --est_method RSEM --out_prefix genes '
            merge_tran_cmd = self.fpkm + ' --est_method RSEM --out_prefix transcripts '
        else:
            merge_gene_cmd = self.tpm + ' --est_method RSEM --out_prefix genes '
            merge_tran_cmd = self.tpm + ' --est_method RSEM --out_prefix transcripts '
        for f in files:
            if re.search(r'genes\Wresults$', f):
                merge_gene_cmd += '{} '.format(self.option('rsem_files').prop['path'] + '/' + f)
            elif re.search(r'isoforms\Wresults$', f):
                merge_tran_cmd += '{} '.format(self.option('rsem_files').prop['path'] + '/' + f)
        gene_com = self.add_command("merge_gene_cmd", merge_gene_cmd).run()
        self.wait(gene_com)
        if gene_com.return_code == 0:
            self.logger.info("运行merge_gene_cmd成功")
        else:
            self.logger.info("运行merge_gene_cmd出错")
            raise Exception("运行merge_gene_cmd出错")
        tran_com = self.add_command("merge_tran_cmd", merge_tran_cmd).run()
        self.wait(tran_com)
        if tran_com.return_code == 0:
            self.logger.info("运行merge_tran_cmd成功")
        else:
            self.logger.info("运行merge_tran_cmd出错")
            raise Exception("运行merge_tran_cmd出错")

    def run_class_code(self):
        if os.path.exists(self.work_dir + "/changed_id_merged.gtf"):
            os.remove(self.work_dir + "/changed_id_merged.gtf")
        os.link(self.option("gtf_merged").prop['path'], self.work_dir + "/changed_id_merged.gtf")
        ## removed by shicaiping, in order to keep identical with annotation
        #os.system(""" sed -i "s/gene://g" %s """ % (self.work_dir + "/changed_id_merged.gtf"))
        #os.system(""" sed -i "s/transcript://g" %s """ % (self.work_dir + "/changed_id_merged.gtf"))
        if os.path.exists(self.option("gtf_ref").prop['path']) and os.path.exists(
                        self.work_dir + "/changed_id_merged.gtf"):
            class_code_path = self.work_dir + '/class_code'
            write_relation_text_for_merge_gtf(self.option("gtf_ref").prop['path'],
                                              self.work_dir + "/changed_id_merged.gtf", class_code_path)
            self.option("class_code").set_path(self.work_dir + "/class_code")

    def get_distribution(self):
        """获取表达量分布图的作图数据"""
        # gene
        distribution(rfile='./gene_distribution.r', input_matrix=self.option('gene_fpkm').prop['path'],
                     outputfile=self.work_dir, filename="gene")
        # transcript
        distribution(rfile='./tran_distribution.r', input_matrix=self.option('tran_fpkm').prop['path'],
                     outputfile=self.work_dir, filename="transcript")
        gcmd = self.r_path1 + " gene_distribution.r"
        tcmd = self.r_path1 + " tran_distribution.r"
        self.logger.info("开始运行表达量分布图的数据分析")
        cmd1 = self.add_command("gene_cmd", gcmd).run()
        cmd2 = self.add_command("tran_cmd", tcmd).run()
        self.wait()
        if cmd1.return_code == 0 and cmd2.return_code == 0:
            self.logger.info("表达量分布图的数据分析成功")
        else:
            self.set_error("表达量分布图的数据分析出错")

    def get_count_distribution(self):
        """获取count分布图"""
        # gene
        count_distribution = self.work_dir + "/sample_count_distribution"
        if not os.path.exists(count_distribution):
            os.mkdir(count_distribution)
        distribution(rfile='./gene_count_distribution.r', input_matrix=self.option("gene_count").prop['path'],
                     outputfile=count_distribution, filename='gene')
        # transcript
        distribution(rfile='./tran_count_distribution.r', input_matrix=self.option("tran_count").prop['path'],
                     outputfile=count_distribution, filename='transcript')
        genecmd = self.r_path1 + " gene_count_distribution.r"
        trancmd = self.r_path1 + " tran_count_distribution.r"
        self.logger.info("开始运行count的数据分析")
        genecmd1 = self.add_command("gene_count_cmd", genecmd).run()
        trancmd1 = self.add_command("tran_count_cmd", trancmd).run()
        self.wait()
        if genecmd1.return_code == 0 and trancmd1.return_code == 0:
            self.logger.info("count表达量分布图的数据分析成功！")
        else:
            self.set_error("count表达量分布图的数据分析出错！")

    def group_distribution(self, old_fpkm, new_fpkm, old_count, new_count, sample_group_info, outputfile, filename,
                           _type=None):  # add by khl 
        import random
        group_express(old_fpkm=old_fpkm, new_fpkm=new_fpkm, old_count=old_count, new_count=new_count,
                      sample_group_info=sample_group_info, \
                      outputfile=outputfile, filename=filename)
        # new_fpkm
        group_fpkm_distribution = self.work_dir + "/group_{}_fpkm_distribution".format(_type)
        if not os.path.exists(group_fpkm_distribution):
            os.mkdir(group_fpkm_distribution)
        group_count_distribution = self.work_dir + "/group_{}_count_distribution".format(_type)
        if not os.path.exists(group_count_distribution):
            os.mkdir(group_count_distribution)
        distribution(rfile='./group_fpkm_distribution.r', input_matrix=new_fpkm, outputfile=group_fpkm_distribution,
                     filename=filename)
        distribution(rfile='./group_count_distribution.r', input_matrix=new_count, outputfile=group_count_distribution,
                     filename=filename)
        cmd = self.r_path1 + " group_fpkm_distribution.r"
        cmd1 = self.r_path1 + " group_count_distribution.r"

        cmd_run = self.add_command("{}".format((filename + "a").lower()), cmd).run()
        cmd1_run = self.add_command("{}".format((filename + "b").lower()), cmd1).run()

        self.logger.info(cmd)
        self.wait()
        if cmd_run.return_code == 0 and cmd1_run.return_code == 0:
            self.logger.info("计算group密度分布成功")
        else:
            self.logger.info("计算group密度分布失败")

    def group_detail(self):  # add by khl
        g = GroupTableFile()
        if self.option("edger_group").is_set:
            g.set_path(self.option("edger_group").prop['path'])
            sample_group_info = g.get_group_spname()
            self.logger.info("打印sample_group_info信息")
            self.logger.info(sample_group_info)
            outputfile = self.output_dir
            results = os.listdir(outputfile)
            new_group_path = self.work_dir + "/group"
            if not os.path.exists(new_group_path):
                os.mkdir(new_group_path)
            for f in results:
                if re.search(r'^(transcripts\.TMM)(.+)(matrix)$', f):
                    self.logger.info("开始计算trans group分布信息")
                    trans_new_fpkm = self.work_dir + "/Group.trans_" + f
                    trans_new_count = self.work_dir + "/Group.trans_count_" + f
                    self.logger.info("生成trans group 文件")
                    trans_filename = "GroupTrans"
                    old_count = self.work_dir + "/transcripts.counts.matrix"
                    self.group_distribution(old_fpkm=self.work_dir + "/" + f, new_fpkm=trans_new_fpkm,
                                            old_count=old_count, new_count=trans_new_count,
                                            sample_group_info=sample_group_info, outputfile=outputfile,
                                            filename=trans_filename, _type='trans')
                    if os.path.exists(new_group_path + "/Group.trans_" + f):
                        os.remove(new_group_path + "/Group.trans_" + f)
                    os.link(trans_new_fpkm, new_group_path + "/Group.trans_" + f)
                elif re.search(r'^(genes\.TMM)(.+)(matrix)$', f):
                    self.logger.info("开始计算gene group分布信息")
                    genes_new_fpkm = self.work_dir + "/Group.genes_" + f
                    genes_new_count = self.work_dir + "/Group.genes_count_" + f
                    old_count = self.work_dir + "/genes.counts.matrix"
                    self.logger.info("生成gene group 文件")
                    genes_filename = "GroupGenes"
                    self.group_distribution(old_fpkm=self.work_dir + "/" + f, new_fpkm=genes_new_fpkm,
                                            old_count=old_count, new_count=genes_new_count,
                                            sample_group_info=sample_group_info, outputfile=outputfile,
                                            filename=genes_filename, _type='genes')
                    if os.path.exists(new_group_path + "/Group.genes_" + f):
                        os.remove(new_group_path + "/Group.genes_" + f)
                    os.link(genes_new_fpkm, new_group_path + "/Group.genes_" + f)
            self.logger.info("计算基因转录本group成功！")

        else:
            raise Exception("有生物学重复时，请设置样本生物学分组信息！")

    def run_long_transcript(self):
        exon_path = self.option("transcript_fa").prop['path']
        cmd = "{}python {}annotation_longest.py -i {}".format(self.python_path, self.long_path, exon_path)
        self.logger.info("提取最长序列")
        command = self.add_command("the_longest", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运提取最长序列完成")
        else:
            self.set_error("提取最长序列出错")
        output1 = os.path.join(self.work_dir, "exons.fa")
        if os.path.exists(self.output_dir + "/exons.fa"):
            os.remove(self.output_dir + "/exons.fa")
        # shutil.move(output1, self.output_dir + "/exons.fa")
        output2 = os.path.join(self.work_dir, "the_longest_exons.fa")
        if os.path.exists(self.output_dir + "/the_longest_exons.fa"):
            os.remove(self.output_dir + "/the_longest_exons.fa")
        # shutil.move(output2, self.output_dir + "/the_longest_exons.fa")
        # self.option('query', self.work_dir + '/output/exons.fa')

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :params: is_duplicate, 是否计算生物学重复
        :return:
        """
        self.logger.info("设置merge_rsem结果目录")
        results = os.listdir(self.work_dir)
        gene_list_path = self.output_dir + "/gene_list"
        trans_list_path = self.output_dir + "/trans_list"
        if self.option("is_class_code"):
            old_rsem_path = os.path.join(self.work_dir, "oldrsem")
            if not os.path.exists(old_rsem_path):
                os.mkdir(old_rsem_path)
            ref_rsem_path = os.path.join(self.work_dir,'ref_rsem')
            if not os.path.exists(ref_rsem_path):
                os.mkdir(ref_rsem_path)
            else:
                shutil.rmtree(ref_rsem_path)
                os.mkdir(ref_rsem_path)
        ss = 0
        for f in results:
            if re.search(r'^(transcripts\.TMM)(.+)(matrix)$', f):
                ss += 1
                #os.system(""" sed -i "s/gene://g" %s """ % (self.work_dir + '/' + f))
                #os.system(""" sed -i "s/transcript://g" %s """ % (self.work_dir + '/' + f))
                if os.path.exists(self.output_dir + '/' + f):
                    os.remove(self.output_dir + '/' + f)
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                all_gene_list(self.work_dir + '/' + f, trans_list_path)
                input_file = self.output_dir + "/"+f
                if self.option("is_class_code"):
                    filter_ref_gene_or_transcript(input_file, class_code=self.option("class_code").prop['path'],
                                                  output_path=ref_rsem_path, query_type='transcript', gene_list=True)
                    add_gene_name(self.work_dir + '/' + f, os.path.join(old_rsem_path, f),
                                  self.work_dir + "/class_code", type="transcript")
                    self.logger.info("修改第{}个文件".format(str(ss)))
                self.option('tran_fpkm').set_path(self.output_dir + '/' + f)
            elif re.search(r'^(genes\.TMM)(.+)(matrix)$', f):
                ss += 1
                #os.system(""" sed -i "s/gene://g" %s """ % (self.work_dir + '/' + f))
                #os.system(""" sed -i "s/transcript://g" %s """ % (self.work_dir + '/' + f))
                if os.path.exists(self.output_dir + '/' + f):
                    os.remove(self.output_dir + '/' + f)
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                all_gene_list(self.work_dir + '/' + f, gene_list_path)
                input_file = self.output_dir + "/" + f
                if self.option("is_class_code"):
                    filter_ref_gene_or_transcript(input_file, class_code=self.option("class_code").prop['path'],
                                                  output_path=ref_rsem_path, query_type='gene', gene_list=True)
                    add_gene_name(self.work_dir + '/' + f, os.path.join(old_rsem_path, f),
                                  self.work_dir + "/class_code", type="gene")
                self.option('gene_fpkm').set_path(self.output_dir + '/' + f)
            elif re.search(r'^(transcripts\.count)(.+)(matrix)$', f):
                ss += 1
                #os.system(""" sed -i "s/gene://g" %s """ % (self.work_dir + '/' + f))
                #os.system(""" sed -i "s/transcript://g" %s """ % (self.work_dir + '/' + f))
                if os.path.exists(self.output_dir + '/' + f):
                    os.remove(self.output_dir + '/' + f)
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                input_file = self.output_dir + "/" + f
                if self.option("is_class_code"):
                    filter_ref_gene_or_transcript(input_file, class_code=self.option("class_code").prop['path'],
                                                  output_path=ref_rsem_path, query_type='transcript', gene_list=True)
                    add_gene_name(self.work_dir + '/' + f, os.path.join(old_rsem_path, f),
                                  self.work_dir + "/class_code", type="transcript")
                self.option('tran_count').set_path(self.output_dir + '/' + f)
            elif re.search(r'^(genes\.count)(.+)(matrix)$', f):
                ss += 1
                #os.system(""" sed -i "s/gene://g" %s """ % (self.work_dir + '/' + f))
                #os.system(""" sed -i "s/transcript://g" %s """ % (self.work_dir + '/' + f))
                if os.path.exists(self.output_dir + '/' + f):
                    os.remove(self.output_dir + '/' + f)
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                input_file = self.output_dir + "/" + f
                if self.option("is_class_code"):
                    # os.system(""" sed -i "s/gene://g" %s """ % (self.work_dir + '/' + f))
                    # os.system(""" sed -i "s/transcript://g" %s """ % (self.work_dir + '/' + f))
                    filter_ref_gene_or_transcript(input_file, class_code=self.option("class_code").prop['path'],
                                                  output_path=ref_rsem_path, query_type='gene', gene_list=True)
                    add_gene_name(self.work_dir + '/' + f, os.path.join(old_rsem_path, f),
                                  self.work_dir + "/class_code", type="gene")
                self.option('gene_count').set_path(self.output_dir + '/' + f)
        self.logger.info("设置merge_rsem分析结果目录成功")
        if self.option("is_class_code"):
            """表达量RNA分析流程时设置此参数为False"""
            if not os.path.exists(self.work_dir+"/new_gene_location"):
                os.mkdir(self.work_dir+"/new_gene_location")
            new_gene_location(changed_id_gtf=self.option("gtf_merged").prop['path'], output_path=self.work_dir+"/new_gene_location", filename="new_gene_location")
            self.logger.info("提取新基因坐标信息完毕!")

    def run(self):
        super(MergeRsemTool, self).run()
        if self.option("is_class_code"):
            self.run_class_code()
            self.logger.info("开始生成class_code信息！")
        self.merge_rsem()
        self.set_output()
        self.get_distribution()
        self.get_count_distribution()
        if self.option("is_duplicate"):
            self.group_detail()
        # self.run_long_transcript()  #暂时不用再提供此功能
        self.end()

