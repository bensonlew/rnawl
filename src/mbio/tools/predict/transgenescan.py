# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os
import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class TransgenescanAgent(Agent):
    """
    宏转录组 TransGeneScan 进行基因预测
    一般需要用优质的组装拼接结果
    """

    def __init__(self, parent):
        super(TransgenescanAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件/contig文件
            {"name": "gene_faa", "type": "outfile", "format": "sequence.fasta"}, # 输出文件faa文件
            {"name": "gene_ffn", "type": "outfile", "format": "sequence.fasta"}, #输出文件ffn文件
            {"name": "sample", "type": "string"} ##输入样本的名称
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome", code="33301201")
        if not self.option("sample"):
            raise OptionError("必须设置参数sample", code="33301202")

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = 2
        file_size = os.path.getsize(self.option("input_genome").prop['path']) / (1024*1024)
        if file_size <= 1024:
            self._memory = '50G'
        else:
            self._memory = '100G'

    def end(self):
        """
        结束
        :return:
        """
        super(TransgenescanAgent, self).end()


class TransgenescanTool(Tool):
    """
    tool 运行
    """
    def __init__(self, config):
        super(TransgenescanTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.perl_path = "/program/perl/perls/perl-5.24.0/bin/perl"
        self.transgenescan_path = self.config.SOFTWARE_DIR +"/bioinfo/metaGenomic/TransGeneScan-1.2.1/TransGeneScan1.2.1/"
        self.transgenescan_stat = self.config.PACKAGE_DIR + "/metagenomic/scripts/"

    def run_transgenescan(self):
        """
        运行TransGeneScan软件进行预测
        :return:
        """
        self.logger.info("运行TransGeneScan软件进行预测！")
        output_prefix = os.path.join(self.work_dir, 'predict.gene')
        cmd = "{} {}run_TransGeneScan.pl -in {} -out {}".format(self.perl_path, self.transgenescan_path, self.genome_fasta, output_prefix)
        self.logger.info('运行TransGeneScan软件，进行预测')
        self.logger.info(cmd)
        command = self.add_command("transgenescan_cmd", cmd).run()
        self.wait(command)
        if command.return_code in [0]:
            self.logger.info("TransGeneScan运行完成")
        else:
            self.set_error("TransGeneScan运行出错!", code="33301201")

    def run_stat(self):
        """
        对预测结果进行统计
        gene_count total_length average_length max_length min_length GC_content
        :return:
        """
        self.logger.info("开始对预测结果文件进行统计")
        if os.path.exists(os.path.join(self.work_dir, "predict.gene.faa")) and os.path.exists(os.path.join(self.work_dir, "predict.gene.ffn")):
            self.logger.info("对预测结果ffn和faa文件进行整理")
            os.system("sed -i s/_+/_plus/ {}".format(os.path.join(self.work_dir, "predict.gene.faa")))
            os.system("sed -i s/_-/_minus/ {}".format(os.path.join(self.work_dir, "predict.gene.faa")))
            os.system("sed -i s/_+/_plus/ {}".format(os.path.join(self.work_dir, "predict.gene.ffn")))
            os.system("sed -i s/_-/_minus/ {}".format(os.path.join(self.work_dir, "predict.gene.ffn")))
        else:
            self.set_error("TransGeneScan未能生成正确的结果文件", code="33301202")
        ffn_path = os.path.join(self.work_dir, "predict.gene.ffn")
        summary_prefix = os.path.join(self.work_dir, "predict.gene")
        cmd = "{} {}transgenescan_stat.pl {} {}".format(self.perl_path, self.transgenescan_stat, ffn_path, summary_prefix)
        self.logger.info(cmd)
        command = self.add_command("stat_cmd", cmd).run()
        self.wait(command)
        if command.return_code in [0]:
            self.logger.info("TransGeneScan统计运行完成")
        else:
            self.set_error("TransGeneScan统计运行失败!", code="33301203")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        sample_dir = os.path.join(self.output_dir, self.option("sample"))
        if os.path.exists(sample_dir):
            shutil.rmtree(sample_dir)
            os.mkdir(sample_dir)
        else:
            os.mkdir(sample_dir)
        if os.path.exists(os.path.join(self.output_dir, self.option("sample"), self.option("sample") + '.gene.faa')):
            os.remove(os.path.join(self.output_dir,self.option("sample"), self.option("sample") + '.gene.faa'))
        os.link(os.path.join(self.work_dir, 'predict.gene.faa'), os.path.join(self.output_dir, self.option("sample"),self.option("sample") + '.gene.faa'))
        if os.path.exists(os.path.join(self.output_dir,self.option("sample"), self.option("sample") + '.gene.ffn')):
            os.remove(os.path.join(self.output_dir,self.option("sample"), self.option("sample") + '.gene.ffn'))
        os.link(os.path.join(self.work_dir, 'predict.gene.ffn'), os.path.join(self.output_dir, self.option("sample"), self.option("sample") + '.gene.ffn'))
        if os.path.exists(os.path.join(self.output_dir,self.option("sample"), self.option("sample") + '.summary.xls')):
            os.remove(os.path.join(self.output_dir, self.option("sample"), self.option("sample") + '.summary.xls'))
        os.link(os.path.join(self.work_dir, 'predict.gene.summary'), os.path.join(self.output_dir,self.option("sample"), self.option("sample") + '.summary.xls'))
        self.option('gene_ffn',   os.path.join(self.output_dir, self.option("sample"),self.option("sample") + '.gene.ffn'))
        self.option('gene_faa',   os.path.join(self.output_dir, self.option("sample"),self.option("sample") + '.gene.faa'))
        self.logger.info("设置结果文件目录完成")

    def run(self):
        """
        运行
        :return:
        """
        super(TransgenescanTool, self).run()
        self.run_transgenescan()
        self.run_stat()
        self.set_output()
        self.end()
