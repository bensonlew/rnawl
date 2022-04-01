# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.express.express_distribution import *
import os
import re, shutil
from mbio.packages.ref_rna.express.single_sample import group_express
from mbio.files.meta.otu.group_table import *


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
            {"name": "rsem_files", "type": "infile", "format": "rna.rsem_dir"},  # SE测序，包含所有样本的fq文件的文件夹
            {"name": "tran_count", "type": "outfile", "format": "denovo_rna_v2.express_matrix"},
            {"name": "gene_count", "type": "outfile", "format": "denovo_rna_v2.express_matrix"},
            {"name": "tran_fpkm", "type": "outfile", "format": "denovo_rna_v2.express_matrix"},
            {"name": "gene_fpkm", "type": "outfile", "format": "denovo_rna_v2.express_matrix"},
            {"name": "is_duplicate", "type": "bool"}, # 是否有生物学重复 
            {"name": "edger_group", "type":"infile", "format":"meta.otu.group_table"},
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
            raise OptionError('必须设置输入文件：rsem结果文件', code="31900801")
        if self.option("exp_way") not in ['fpkm', 'tpm']:
            raise OptionError("所设表达量的代表指标不在范围内，请检查", code="31900802")
        if self.option("is_duplicate"):
            if not self.option("edger_group").is_set:
                self.set_error("有生物学重复时，请设置样本生物学分组信息！", code="31900801")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = ''

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
        self.r_path1 = "/program/R-3.3.1/bin/Rscript"
        self.distribution_path = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/denovo_rna/express'
        
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
        self.logger.info(merge_tran_cmd)
        self.logger.info(merge_gene_cmd)
        self.logger.info("开始运行merge_gene_cmd")
        self.logger.info("开始运行merge_tran_cmd")
        gene_com = self.add_command("merge_gene_cmd", merge_gene_cmd).run()
        self.wait(gene_com)
        if gene_com.return_code == 0:
            self.logger.info("运行merge_gene_cmd成功")
        else:
            self.logger.info("运行merge_gene_cmd出错")
            self.set_error("运行merge_gene_cmd出错", code="31900802")
        tran_com = self.add_command("merge_tran_cmd", merge_tran_cmd).run()
        self.wait(tran_com)
        if tran_com.return_code == 0:
            self.logger.info("运行merge_tran_cmd成功")
        else:
            self.logger.info("运行merge_tran_cmd出错")
            self.set_error("运行merge_tran_cmd出错", code="31900803")

    def get_distribution(self):
        """获取表达量分布图的作图数据"""
        # gene
        distribution(rfile='./gene_distribution.r', input_matrix=self.option('gene_fpkm').prop['path'], outputfile=self.work_dir, filename="gene")
        # transcript
        distribution(rfile='./tran_distribution.r', input_matrix=self.option('tran_fpkm').prop['path'], outputfile=self.work_dir, filename="transcript")
        gcmd = self.r_path1 + " gene_distribution.r"
        tcmd = self.r_path1 + " tran_distribution.r"
        self.logger.info("开始运行表达量分布图的数据分析")
        cmd1 = self.add_command("gene_cmd", gcmd).run()
        cmd2 = self.add_command("tran_cmd", tcmd).run()
        self.wait()
        if cmd1.return_code == 0 and cmd2.return_code == 0:
            self.logger.info("表达量分布图的数据分析成功")
        else:
            self.set_error("表达量分布图的数据分析出错", code="31900804")
    
    def group_distribution(self, old_fpkm, new_fpkm,sample_group_info, outputfile, filename):
        import random
        shutil.copy2(self.distribution_path+"/express_distribution.py",self.work_dir+"/express_distribution.py")
        shutil.copy2(self.distribution_path+"/express_distribution.r",self.work_dir+"/express_distribution.r")
        self.logger.info("express_distribution py和r文件复制完毕！")
        rfile = self.work_dir+"/express_distribution.r"
        group_express(old_fpkm = old_fpkm, new_fpkm = new_fpkm, sample_group_info = sample_group_info, \
                        rfile=rfile, outputfile = outputfile, filename = filename)
        cmd = self.r_path1 + " {}".format(rfile)
        cmd1 = self.add_command("{}".format("".join(random.sample("abcdeghijk",3))), cmd).run()
        self.logger.info(cmd)
        self.wait()
        if cmd1.return_code == 0:
            self.logger.info("计算group密度分布成功")
        else:
            self.logger.info("计算group密度分布失败")
    
    def group_detail(self):
        g = GroupTableFile()
        if self.option("edger_group").is_set:
            g.set_path(self.option("edger_group").prop['path'])
            sample_group_info = g.get_group_spname()
            self.logger.info("打印sample_group_info信息")
            self.logger.info(sample_group_info)
            outputfile = self.work_dir
            results=os.listdir(outputfile)
            new_group_path = self.work_dir + "/group"
            if not os.path.exists(new_group_path):
                os.mkdir(new_group_path)
            for f in results:
                if re.search(r'^(transcripts\.TMM)(.+)(matrix)$', f):
                    self.logger.info("开始计算trans group分布信息")
                    trans_new_fpkm = self.work_dir + "/Group.trans_"+f
                    self.logger.info("生成trans group 文件")
                    trans_filename = "GroupTrans"
                    self.group_distribution(old_fpkm=self.work_dir + "/"+f, new_fpkm=trans_new_fpkm,\
                            sample_group_info=sample_group_info, outputfile=outputfile, filename = trans_filename)
                    shutil.copy2(trans_new_fpkm, new_group_path+"/Group.trans_"+f)
                elif re.search(r'^(genes\.TMM)(.+)(matrix)$', f):
                    self.logger.info("开始计算gene group分布信息")
                    genes_new_fpkm = self.work_dir + "/Group.genes_"+f
                    self.logger.info("生成gene group 文件")
                    genes_filename = "GroupGenes"
                    self.group_distribution(old_fpkm=self.work_dir + "/"+f, new_fpkm=genes_new_fpkm,\
                            sample_group_info=sample_group_info, outputfile=outputfile, filename = genes_filename)
                    shutil.copy2(genes_new_fpkm, new_group_path+"/Group.genes_"+f)
            self.logger.info("计算基因转录本group成功！")
        else:
            self.set_error("有生物学重复时，请设置样本生物学分组信息！", code="31900805")
        
    def set_output(self, is_duplicate = None):
        """
        将结果文件link到output文件夹下面
        :params: is_duplicate, 是否计算生物学重复
        :return:
        """
        
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置merge_rsem结果目录")
        results = os.listdir(self.work_dir)
        try:
            for f in results:
                if re.search(r'^(transcripts\.TMM)(.+)(matrix)$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                    self.option('tran_fpkm').set_path(self.output_dir + '/' + f)
                elif re.search(r'^(genes\.TMM)(.+)(matrix)$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                    self.option('gene_fpkm').set_path(self.output_dir + '/' + f)
                elif re.search(r'^(transcripts\.count)(.+)(matrix)$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                    self.option('tran_count').set_path(self.output_dir + '/' + f)
                elif re.search(r'^(genes\.count)(.+)(matrix)$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                    self.option('gene_count').set_path(self.output_dir + '/' + f)
            self.logger.info("设置merge_rsem分析结果目录成功")
        except Exception as e:
            self.logger.info("设置merge_rsem分析结果目录失败{}".format(e))
        
            
    def run(self):
        super(MergeRsemTool, self).run()
        self.merge_rsem()
        self.set_output(is_duplicate=False)
        self.get_distribution()
        if self.option("is_duplicate"):
            self.group_detail()
        self.end()
