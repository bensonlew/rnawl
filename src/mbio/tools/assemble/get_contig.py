# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue & guhaidong'

import os
# import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class GetContigAgent(Agent):
    """
    宏基因组装之序列处理，去掉低质量的序列
    version: v1.0
    author: wangzhaoyue & guhaidong
    last_modify: 2017.09.12
    """

    def __init__(self, parent):
        super(GetContigAgent, self).__init__(parent)
        options = [
            {"name": "scafSeq", "type": "infile", "format": "sequence.fasta"},  # 输入文件,sample.scafSeq
            {"name": "min_contig", "type": "string", "default": "500"},  # 输入最短contig长度，默认500
            {"name": "scaftig", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，scaffold去掉N后的序列
            {"name": "cut_more_scaftig", "type": "outfile", "format": "sequence.fasta"},
            # 输出文件，去掉小于最短contig长度的序列
            # {"name": "scaftig_stat", "type": "outfile", "format": "sequence.profile_table"},  # 输出文件，对组装后的序列进行信息统计
        ]
        self.add_option(options)
        self.step.add_steps("GetContig")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.GetContig.start()
        self.step.update()

    def stepfinish(self):
        self.step.GetContig.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('scafSeq'):
            raise OptionError('必须输入样本的scafSeq文件', code="31300801")
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "2G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            [".scagtig", "fasta", "组装处理后的结果文件"]
        ])
        super(GetContigAgent, self).end()


class GetContigTool(Tool):
    def __init__(self, config):
        super(GetContigTool, self).__init__(config)
        # self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.perl_path = '/miniconda2/bin/perl '
        self.logger.info(self.perl_path)
        # self.get_scaftig_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/get_scaftig.pl '
        self.get_scaftig_path = self.config.PACKAGE_DIR + '/assemble/scripts/get_scaftig.pl '
        # self.cut_more_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/cut_more.pl '
        self.cut_more_path = self.config.PACKAGE_DIR + '/sequence/scripts/cut_more.pl '
        # self.scaftig_stat_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/contig_stat.pl '

    def run(self):
        """
        运行
        :return:
        """
        super(GetContigTool, self).run()
        self.run_get_scaftig()
        self.run_cut_more()
        # self.run_scaftig_stat()  #工作流更改，不进行质量统计
        self.set_output()
        self.end()

    def run_get_scaftig(self):
        """
        perl get_scaftig.pl [scafseq.file] [输出文件的名称前缀]
        :return:
        """
        self.logger.info("path is : " + self.option('scafSeq').path)
        # sample_name = os.path.basename(self.option('scafSeq').prop['path']).split('.scafSeq')[0]
        sample_name = os.path.basename(self.option('scafSeq').path).split('.scafSeq')[0]
        cmd = self.perl_path + self.get_scaftig_path + self.option('scafSeq').path + ' ' + sample_name
        # cmd = self.perl_path + self.get_scaftig_path + self.option('scafSeq').prop['path'] + ' ' + sample_name
        command = self.add_command("get_scaftig", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行get_scaftig完成")
        else:
            self.set_error("get_scaftig运行出错!", code="31300801")

    def run_cut_more(self):
        """
        perl cut_more.pl [run_get_scaftig的输出文件] [最短contig长度] [输出文件的名称前缀]
        :return:
        """
        sample_name = os.path.basename(self.option('scafSeq').prop['path']).split('.scafSeq')[0]
        scaftig_file = self.work_dir + "/" + sample_name + ".scaftig"
        outfile_name = os.path.basename(scaftig_file)
        cmd = self.perl_path + self.cut_more_path + scaftig_file + ' ' + self.option('min_contig') + ' ' + outfile_name
        command = self.add_command("cut_more", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行cut_more完成")
        else:
            self.set_error("运行cut_more运行出错!", code="31300802")

    def run_scaftig_stat(self):
        """
        perl contig_stat.pl [run_cut_more的输出文件] [最短contig长度] [输出文件的名称全称]
        :return:
        """
        '''
        sample_name = os.path.basename(self.option('scafSeq').prop['path']).split('.scafSeq')[0]
        kmer = sample_name.strip().split("_K")[-1]
        sample_ID = sample_name.strip().split("_K" + kmer)[0]
        cut_more_scaftig_file = self.work_dir + "/" + sample_name + ".scaftig.more" + self.option('min_contig')
        outfile_name = os.path.basename(cut_more_scaftig_file)
        cmd = self.perl_path + self.scaftig_stat_path + cut_more_scaftig_file + ' ' + self.option('min_contig') + ' ' + \
              sample_ID + ' ' + kmer + ' ' + outfile_name + '.stat'
        command = self.add_command("scaftig_stat", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行scaftig_stat完成")
        else:
            self.set_error("运行scaftig_stat运行出错!")
        '''
        pass

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        sample_name = os.path.basename(self.option('scafSeq').prop['path']).split('.scafSeq')[0]
        # shutil.copy2(self.work_dir + "/" + sample_name + ".scaftig", self.output_dir + "/" + sample_name + ".scaftig")
        if os.path.exists(self.output_dir + "/" + sample_name + ".contig.fa"):
            os.remove(self.output_dir + "/" + sample_name + ".contig.fa")
        os.link(self.work_dir + "/" + sample_name + ".scaftig.more" + self.option('min_contig'), self.output_dir +
                "/" + sample_name + ".contig.fa")
        # shutil.copy2(self.work_dir + "/" + sample_name + ".scaftig.more" + self.option('min_contig') + '.stat', self.output_dir +
        #             "/" + sample_name + ".scaftig.more" + self.option('min_contig') + '.stat')
        self.logger.info("设置get_contig分析结果目录成功")
