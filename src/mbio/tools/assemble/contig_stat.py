# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna.trans_step import step_count

class ContigStatAgent(Agent):
    """
    进行单样品contig统计，对SOAPdenovo2拼接结果挑选出最佳kmer结果
    version: v1.0
    author: guhaidong
    last_modify: 2017.08.23
    """
    def __init__(self, parent):
        super(ContigStatAgent, self).__init__(parent)
        options = [
            {"name": "contig_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 输入文件，对组装后的序列路径
            {"name": "choose_kmer", "type": "string", "default": "false"},  # 是否进行kmer筛选，默认不筛选
            {"name": "min_contig", "type": "string", "default": "300"},  #是否只统计序列长度高于此值的结果
            {"name": "assembly_stat", "type": "string"}#, "format": "sequence.profile_table"},  # 输出文件，统计组装后的结果文件
        ]
        self.add_option(options)
        self.step.add_steps("ContigStat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.ContigStat.start()
        self.step.update()

    def stepfinish(self):
        self.step.ContigStat.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('contig_dir'):
            raise OptionError('必须输入拼接结果文件夹', code="31300201")
        if not self.option('assembly_stat'):
            raise OptionError('必须输入结果表名称', code="31300202")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
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
            [".contig.fa", "fasta", "组装处理后的结果文件"]
        ])
        super(ContigStatAgent, self).end()


class ContigStatTool(Tool):
    def __init__(self, config):
        super(ContigStatTool, self).__init__(config)
        self.Python_path = '/program/Python/bin/python '
        # self.stat_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/contig_stat.py '
        self.stat_path = self.config.PACKAGE_DIR + '/assemble/contig_stat.py '

    def run(self):
        """
        运行
        :return:
        """
        super(ContigStatTool, self).run()
        self.get_assemble_stat()
        #self.step_count()
        self.set_output()
        self.end()

    def get_assemble_stat(self):
        """
        将每个组装的统计数据文件按照顺序汇总成一张表格,并挑选出最佳的组装结果
        :return:
        """
        #assemble_stat = self.work_dir + '/all_assembly.stat'
        #select_stat = self.work_dir + '/select.stat'
        final_stat = self.work_dir + '/' + self.option('assembly_stat')#'/assembly.stat'
        if self.option('choose_kmer') == "false":
            self.logger.info("运行assemble_stat,但不进行kmer优选")
            #cmd = self.Python_path + self.stat_path + ' -contig_dir %s -final_stat %s'\
            #                                      % (self.option('contig_dir').prop['path'], final_stat)
            cmd = self.Python_path + self.stat_path + ' -contig_dir %s -final_stat %s -min_contig %s'\
                                                  % (self.option('contig_dir').prop['path'], final_stat,
                                                     self.option("min_contig"))
        else:
            self.logger.info("运行assemble_stat，且kmer优选")
            cmd = self.Python_path + self.stat_path + ' -contig_dir %s -select_kmer %s -final_stat %s'\
                                                  % (self.option('contig_dir').prop['path'], self.option('choose_kmer'), final_stat)
        self.logger.info("cmd logger:" + cmd)
        command = self.add_command("assemble_stat", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行assemble_stat完成")
        else:
            self.set_error("assemble_stat运行出错!", code="31200201")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        all_files = os.listdir(self.work_dir)
        for files in all_files:
            #file_name = os.path.basename(files)
            if files.find("contig.fa") != -1:
                if os.path.exists(self.output_dir + '/' + files):
                    self.logger.info("have files:" + files)
                    os.remove(self.output_dir + '/' + files)
                os.link(self.work_dir + '/' + files, self.output_dir + '/' + files)
        if os.path.exists(self.output_dir + '/' + self.option('assembly_stat')):
            os.remove(self.output_dir + '/' + self.option('assembly_stat'))
        os.link(self.work_dir + '/' + self.option('assembly_stat'), self.output_dir + '/' + self.option('assembly_stat'))
        self.logger.info("设置assemble_stat分析结果目录成功")