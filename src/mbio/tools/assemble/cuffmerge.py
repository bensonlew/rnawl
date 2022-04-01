# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil
import re


class CuffmergeAgent(Agent):
    """
    有参转录组cuffmerge合并
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2016.09.12
    """
    def __init__(self, parent):
        super(CuffmergeAgent, self).__init__(parent)
        options = [
            {"name": "assembly_GTF_list.txt", "type": "infile", "format": "assembly.merge_txt"},
            # 所有样本比对之后的bam文件
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 参考基因文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因的注释文件
            {"name": "cpu", "type": "int", "default": 10},  # cufflinks软件所分配的cpu数
            {"name": "merged_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 输出的合并文件
        ]
        self.add_option(options)
        self.step.add_steps("cuffmerge")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.cuffmerge.start()
        self.step.update()

    def stepfinish(self):
        self.step.cuffmerge.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('assembly_GTF_list.txt'):
            raise OptionError('必须输入所有样本gtf路径文件为txt格式')
        if not self.option('ref_fa'):
            raise OptionError('必须输入参考序列ref.fa')
        if not self.option('ref_gtf'):
            raise OptionError('必须输入参考序列ref.gtf')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["merged.gtf", "gtf", "样本合并之后的gtf文件"]
        ])
        super(CuffmergeAgent, self).end()


class CuffmergeTool(Tool):
    def __init__(self, config):
        super(CuffmergeTool, self).__init__(config)
        self._version = "v1.0.1"
        self.cuffmerge_path = 'bioinfo/rna/cufflinks-2.2.1/'
        tmp = os.path.join(self.config.SOFTWARE_DIR, self.cuffmerge_path)
        tmp_new = tmp + ":$PATH"
        self.logger.debug(tmp_new)
        self.set_environ(PATH=tmp_new)

    def run(self):
        """
        运行
        :return:
        """
        super(CuffmergeTool, self).run()
        self.run_cuffmerge()
        self.run_gtf_to_fa()
        self.set_output()
        self.end()

    def run_cuffmerge(self):
        """
        运行cufflinks软件，进行拼接合并
        """
        cmd = self.cuffmerge_path + ('cuffmerge -p {} -g {} -s {}  -o {}merge_out'.format(
            self.option('cpu'), self.option('ref_gtf').prop['path'], self.option('ref_fa').prop['path'],
            self.work_dir + "/")) + ' %s' % (self.option('assembly_GTF_list.txt').prop['path'])
        self.logger.info('运行cufflinks软件，进行拼接合并')
        command = self.add_command("cuffmerge_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cuffmerge运行完成")
        else:
            self.set_error("cuffmerge运行出错!")

    def run_gtf_to_fa(self):
        """
        运行gtf_to_fasta，转录本gtf文件转fa文件
        """
        cmd = self.cuffmerge_path + "gffread %s -g %s -w merged.fa" % (
        self.work_dir + "/merge_out/"+"merged.gtf", self.option('ref_fa').prop['path'])
        self.logger.info('运行gtf_to_fasta，形成fasta文件')
        command = self.add_command("gtf_to_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            if os.path.exists(self.output_dir + "/merged.gtf"):
                os.remove(self.output_dir + "/merged.gtf")
            os.link(self.work_dir + "/merge_out/merged.gtf", self.output_dir + "/merged.gtf")
            self.option('merged_gtf').set_path(self.work_dir + "/" + "merge_out"+"/"+"merged.gtf")
            if os.path.exists(self.output_dir + "/merged.fa"):
                os.remove(self.output_dir + "/merged.fa")
            os.link(self.work_dir + "/merged.fa", self.output_dir + "/merged.fa")
            self.logger.info("设置拼接合并分析结果目录成功")

        except Exception as e:
            self.logger.info("设置拼接合并分析结果目录失败{}".format(e))
            self.set_error("设置拼接合并分析结果目录失败{}".format(e))
