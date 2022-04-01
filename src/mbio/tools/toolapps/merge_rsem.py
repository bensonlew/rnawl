# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.express.express_distribution import *
import os
import re


class MergeRsemAgent(Agent):
    """
    调用align_and_estimate_abundance.pl脚本，运行rsem，进行表达量计算分析
    宏转录组 合并样本结果，借用Trinity软件统计的脚本
    样本数必须大于2才可以，否则脚本会报错
    """
    def __init__(self, parent):
        super(MergeRsemAgent, self).__init__(parent)
        options = [
            {"name": "rsem_files", "type": "infile", "format": "rna.rsem_dir"},  # SE测序，包含所有样本的fq文件的文件夹
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
        参数二次检查
        :return:
        """
        if not self.option('rsem_files'):
            raise OptionError('必须设置输入文件：rsem结果文件', code="34400901")
        if self.option("exp_way") not in ['fpkm', 'tpm']:
            raise OptionError("所设表达量的代表指标不在范围内，请检查", code="34400902")
        return True

    def set_resource(self):
        """
        设置所需资源
        :return:
        """
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        """
        结束
        """
        super(MergeRsemAgent, self).end()


class MergeRsemTool(Tool):
    """
    调用align_and_estimate_abundance.pl脚本，运行rsem，进行表达量计算分析
    宏转录组 合并样本结果，借用Trinity软件统计的脚本
    """
    def __init__(self, config):
        super(MergeRsemTool, self).__init__(config)
        self._version = '1.0.1'
        self.fpkm = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Trinity-v2.9.1/fpkm/abundance_estimates_to_matrix.pl"
        self.tpm = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Trinity-v2.9.1/tpm/abundance_estimates_to_matrix.pl"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.perl = "/program/perl/perls/perl-5.24.0/bin/perl"
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:$PATH"
        self._r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def merge_rsem(self):
        """
        根据方法计算不同样本的结果
        :return:
        """
        self.logger.info("对RSEM软件的结果进行合并")
        files = os.listdir(self.option('rsem_files').prop['path'])
        if self.option('exp_way') == 'fpkm':
            merge_gene_cmd = self.perl + " "+ self.fpkm + ' --est_method RSEM --out_prefix genes '
            merge_tran_cmd = self.perl + " "+ self.fpkm + ' --est_method RSEM --out_prefix transcripts '
        else:
            merge_gene_cmd = self.perl + " "+ self.tpm + ' --est_method RSEM --out_prefix genes '
            merge_tran_cmd = self.perl + " "+ self.tpm + ' --est_method RSEM --out_prefix transcripts '
        for f in files:
            if re.search(r'genes\Wresults', f):
                merge_gene_cmd += '{} '.format(self.option('rsem_files').prop['path'] + '/' + f)
            elif re.search(r'isoforms\Wresults', f):
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
            self.set_error("运行merge_gene_cmd出错", code="34400901")
        tran_com = self.add_command("merge_tran_cmd", merge_tran_cmd).run()
        self.wait(tran_com)
        if tran_com.return_code == 0:
            self.logger.info("运行merge_tran_cmd成功")
        else:
            self.logger.info("运行merge_tran_cmd出错")
            self.set_error("运行merge_tran_cmd出错", code="34400902")
        
    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置merge_rsem结果目录")
        results = os.listdir(self.work_dir)
        try:
            for f in results:
                if re.search(r'^(transcripts\.TMM)(.+)(matrix)$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                elif re.search(r'^(genes\.TMM)(.+)(matrix)$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                elif re.search(r'^(transcripts\.counts)\.(matrix)$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                elif re.search(r'^(genes\.counts)\.(matrix)$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
            self.logger.info("设置merge_rsem分析结果目录成功")
        except Exception as e:
            self.logger.info("设置merge_rsem分析结果目录失败{}".format(e))
        
            
    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("开始运行tool")
        super(MergeRsemTool, self).run()
        self.merge_rsem()
        self.set_output()
        self.end()
