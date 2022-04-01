# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'

import os
import re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna.trans_step import step_count


class MetagPredictCombineAgent(Agent):
    def __init__(self, parent):
        super(MetagPredictCombineAgent, self).__init__(parent)
        options = [
            {"name": "fna_path", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "faa_path", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "len_range", "type": "string"},
            # 输入文件，预测后的序列路径
            {"name": "sample_stat", "type": "outfile", "format": "sequence.profile_table"},  # 输出文件，对各基因预测结果进行统计
            {"name": "fasta", "type": "outfile", "format": "sequence.fasta"},  # 输出核酸文件，输出序列用于构建非冗余基因集
            {"name": "faa", "type": "outfile", "format": "sequence.fasta"},  # 输出蛋白文件，输出序列用于构建非冗余基因集
            {"name": "fasta_sample", "type": "outfile", "format": "sequence.fasta"},  # 输出核酸文件，输出序列用于构建非冗余基因集
            {"name": "fasta_mix", "type": "outfile", "format": "sequence.fasta"},  # 输出核酸文件，输出序列用于构建非冗余基因集
            {"name": "faa_sample", "type": "outfile", "format": "sequence.fasta"},  # 输出单拼样本蛋白文件，输出序列用于构建非冗余基因集
            {"name": "faa_mix", "type": "outfile", "format": "sequence.fasta"},  # 输出混拼样本蛋白文件，输出序列用于构建非冗余基因集
        ]
        self.add_option(options)
        self.step.add_steps("MetagPredictCombine")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.MetagPredictCombine.start()
        self.step.update()

    def stepfinish(self):
        self.step.MetagPredictCombine.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 5
        #self._memory = "50G"
        self._memory = str((50 + self._rerun_time * 20)) + 'G' # by xieshichang 20200424

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MetagPredictCombineAgent, self).end()


class MetagPredictCombineTool(Tool):
    def __init__(self, config):
        super(MetagPredictCombineTool, self).__init__(config)
        self._version = "1"
        self.python_path = '/program/Python/bin/python '
        # self.gene_stat_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/gene_stat.py'
        self.gene_stat_path = self.config.PACKAGE_DIR + '/gene_structure/gene_stat.py'
        self.set_environ(LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + '/bioinfo/seq/EMBOSS-6.6.0/lib')
        self.emboss_path = "/bioinfo/seq/EMBOSS-6.6.0/emboss/"

    def run(self):
        """
        运行
        :return:
        """
        super(MetagPredictCombineTool, self).run()
        self.run_stat()
        self.step_count()
        self.combine_faa()
        self.set_output()
        self.end()

    def run_stat(self):
        """
        gene_stat -gene_dir  gene.directory -output_stat  stat_file -output_fa  fasta_file
        :return:
        """
        cmd = self.python_path + '{} -gene_dir {} -output_stat {} -output_fa {} -sp_name {}'
        cmd = cmd.format(self.gene_stat_path, self.option('fna_path').prop['path'],
                         self.work_dir + '/sample.metagene.stat', self.work_dir + '/Total.metagene.fa', 1)
        command = self.add_command("metagenestat", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行metagenestat的cmd完成")
        else:
            self.set_error("运行metagenestat的cmd运行出错!", code="32201201")

    def step_count(self):
        """
        调用函数，对挑选出来的序列,整理到文件夹汇总，并进行步长统计
        :return:
        """
        self.logger.info("开始组装部分进行步长统计")
        step_list = self.option('len_range').split(',')
        all_files = self.option('fna_path').fastas_full
        all_files.append(self.work_dir + '/Total.metagene.fa')
        if not os.path.exists(os.path.join(self.output_dir, "len_distribute")):
            os.mkdir(os.path.join(self.output_dir, "len_distribute"))
        for file_full in all_files:
            file_name = os.path.basename(file_full)
            sample_name = file_name.split('.')[0]
            for step in step_list:
                output1 = self.work_dir + '/' + sample_name + '_step_' + str(step) + '.txt'
                output2 = self.output_dir + '/len_distribute/' + sample_name + '_step_' + str(step) + '.final.txt'
                if os.path.exists(output2):
                    os.remove(output2)
                step_count(file_full, output1, 20, int(step), output2)
        self.logger.info("序列步长统计结束")

    def combine_faa(self):
        faa_mix = self.work_dir + '/Total.metagene_mix.faa'
        faa_sample = self.work_dir + '/Total.metagene_sample.faa'
        faa = self.work_dir + '/Total.metagene.faa'
        if os.path.exists(self.work_dir + '/Total.metagene.fa_mix'):
            mix = True
        else:
            mix = False
        for f in [faa_mix, faa_sample, faa]:
            if os.path.exists(f):
                os.remove(f)
        for f in os.listdir(self.option('faa_path').path):
            f_path = os.path.join(self.option('faa_path').path, f)
            f_name = f.partition('.')[0]
            cmd1 = ''
            if f_name in ["newbler", "Megahit_Mix"]:
                cmd1 = 'cat {} >> {}'.format(f_path, faa_mix)
            elif mix:
                cmd1 = 'cat {} >> {}'.format(f_path, faa_sample)
            cmd2 = 'cat {} >> {}'.format(f_path, faa)
            if os.system(cmd2) != 0:
                self.set_error('蛋白序列合并失败')
            if mix and os.system(cmd1) != 0:
                self.set_error('蛋白序列合并失败')

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + '/sample.metagene.stat'):
            os.remove(self.output_dir + '/sample.metagene.stat')
        if os.path.exists(self.output_dir + '/Total.metagene.fa'):
            os.remove(self.output_dir + '/Total.metagene.fa')
        if os.path.exists(self.output_dir + '/Total.metagene.faa'):
            os.remove(self.output_dir + '/Total.metagene.faa')
        if os.path.exists(self.work_dir + '/Total.metagene.faa'):
            os.link(self.work_dir + '/Total.metagene.faa', self.output_dir + '/Total.metagene.faa')
            self.option('faa').set_path(self.output_dir + '/Total.metagene.faa')
        if os.path.exists(self.work_dir + '/sample.metagene.stat'):
            os.link(self.work_dir + '/sample.metagene.stat', self.output_dir + '/sample.metagene.stat')
            self.option('sample_stat').set_path(self.output_dir + '/sample.metagene.stat')
        if os.path.exists(self.work_dir + '/Total.metagene.fa'):
            os.link(self.work_dir + '/Total.metagene.fa', self.output_dir + '/Total.metagene.fa')
            self.option('fasta').set_path(self.output_dir + '/Total.metagene.fa')
        if os.path.exists(self.work_dir + '/Total.metagene_sample.faa'):
            self.option('faa_sample').set_path(self.work_dir + '/Total.metagene_sample.faa')
            self.option('fasta_sample').set_path(self.work_dir + '/Total.metagene.fa_sample')
        if os.path.exists(self.work_dir + '/Total.metagene_mix.faa'):
            self.option('faa_mix').set_path(self.work_dir + '/Total.metagene_mix.faa')
            self.option('fasta_mix').set_path(self.work_dir + '/Total.metagene.fa_mix')
        self.logger.info("设置Metagene分析结果目录成功")
