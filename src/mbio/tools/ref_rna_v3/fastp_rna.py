# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class FastpRnaAgent(Agent):
    def __init__(self, parent=None):
        super(FastpRnaAgent, self).__init__(parent)
        options = [
            {'name': 'fq1', 'type': "infile", "format": "sequence.fastq"},  # 原始数据
            {'name': 'fq2', 'type': "infile", "format": "sequence.fastq"},  # 原始数据
            {'name': 'fqs', 'type': "infile", "format": "sequence.fastq"},  # 原始数据(单端）
            {'name': 'fq_type', 'type': "string", "default": "PE"},  # 测序方式
            # -q,一个碱基合格的质量值,默认15表示phred质量> = Q15是合格的。
            {'name': 'qualified_quality_phred', 'type': "string", "default": "15"},
            # phred+64 OR phred+33 #add by 20200401 fwy
            {"name": "quality_score_system", "type": "string", "default": "phred+33"},
            # -l,长度过滤参数，比此值短的读取将被丢弃
            {'name': 'length_required', "type": "string", "default": "30"},
            # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {'name': "cut_by_quality5", "type": "string"},
            # -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            {'name': "cut_by_quality3", "type": "string"},
            # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {'name': 'cut_mean_quality', "type": "string", "default": "20"},
            # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {'name': 'n_base_limit', "type": "string", "default": "5"},
            # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {'name': 'compression', "type": "string", "default": "2"},
            # -w,线程数
            {'name': 'thread', "type": "string", "default": '4'},
            # --adapter_sequence,the adapter for read1
            {"name": "adapter_sequence", "type": "string", "default": "AGATCGGAAGAGCACACGTC"},
            # --adapter_sequence_r2,the adapter for read2 (PE data only)
            {"name": "adapter_sequence_r2", "type": "string", "default": "AGATCGGAAGAGCGTCGTGT"},
            {"name": "adapter_sequence_s", "type": "string", "default": "AGATCGGAAGAGCACACGTC"}
        ]
        self.add_option(options)
        self.step.add_steps("fastp")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.fastp.start()
        self.step.update()

    def stepfinish(self):
        self.step.fastp.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option('fq1'):
            raise OptionError('必须输入1端序列')
        if not self.option('fq2'):
            raise OptionError('必须输入2端序列')
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 8
        self._memory = '10G'

    def end(self):
        super(FastpRnaAgent, self).end()


class FastpRnaTool(Tool):
    def __init__(self, config):
        super(FastpRnaTool, self).__init__(config)
        self.fastp_path = 'bioinfo/ref_rna_v3/R_4.1/miniconda3/bin/fastp'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + 'bioinfo/ref_rna_v3/R_4.1/miniconda3/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + 'bioinfo/ref_rna_v3/R_4.1/miniconda3/lib')

    def run(self):
        super(FastpRnaTool, self).run()
        self.run_fastp()
        self.set_output()
        self.end()

    def run_fastp(self):
        """
        运行fastp,进行质控
        """
        if self.option('fq_type') == "PE":
            try:
                sample_name = '_'.join(os.path.basename(self.option('fq1').path).split('_')[0:-4])
            except:
                sample_name = os.path.basename(self.option('fq1').prop['path'])
            # 增加判断，用于数据命名不规范是导致sample_name为空的情形
            if sample_name == "":
                sample_name = os.path.basename(self.option('fq1').prop['path']).split(".")[0]
            cmd = self.fastp_path + ' -i %s -I %s -o %s -O %s -q %s -l %s -M %s -n %s -j %s -h %s -w %s' % (
                self.option('fq1').prop['path'], self.option('fq2').prop['path'],
                self.work_dir + "/" + sample_name + '.clean.1.fastq',
                self.work_dir + "/" + sample_name + '.clean.2.fastq',
                self.option('qualified_quality_phred'), self.option('length_required'),
                self.option('cut_mean_quality'), self.option('n_base_limit'),
                self.work_dir + "/" + sample_name + '.json', self.work_dir + "/" + sample_name + '.html',
                self.option('thread'))
            if self.option('cut_by_quality5'):
                cmd += ' -5 %s' % (self.option('cut_by_quality5'))
            if self.option('cut_by_quality3'):
                cmd += ' -3 %s' % (self.option('cut_by_quality3'))
            cmd += " --adapter_sequence {} --adapter_sequence_r2 {}".format(self.option("adapter_sequence"),
                                                                            self.option("adapter_sequence_r2"))
            if self.option("quality_score_system").endswith("64"):
                cmd += " -6"
            self.logger.info('运行fastp，进行质控')
            command = self.add_command("fastp_cmd", cmd).run()
        else:
            try:
                sample_name = '_'.join(os.path.basename(self.option('fqs').prop['path']).split('_')[0:-4])
            except:
                sample_name = os.path.basename(self.option('fqs').prop['path'])
            # 增加判断，用于数据命名不规范是导致sample_name为空的情形
            if sample_name == "":
                sample_name = os.path.basename(self.option('fqs').prop['path']).split(".")[0]
            cmd = self.fastp_path + ' -i %s -o %s -q %s -l %s -M %s -n %s -j %s -h %s -w %s' % (
                self.option('fqs').prop['path'],
                self.work_dir + "/" + sample_name + '.clean.s.fastq',
                self.option('qualified_quality_phred'), self.option('length_required'),
                self.option('cut_mean_quality'), self.option('n_base_limit'),
                self.work_dir + "/" + sample_name + '.json', self.work_dir + "/" + sample_name + '.html',
                self.option('thread'))
            if self.option('cut_by_quality5'):
                cmd += ' -5 %s' % (self.option('cut_by_quality5'))
            if self.option('cut_by_quality3'):
                cmd += ' -3 %s' % (self.option('cut_by_quality3'))
            cmd += " --adapter_sequence {} ".format(self.option("adapter_sequence_s"))
            if self.option("quality_score_system").endswith("64"):
                cmd += " -6"
            self.logger.info('运行fastp，进行质控')
            command = self.add_command("fastp_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fastp运行完成")
        else:
            self.set_error("fastp运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            all_files = os.listdir(self.work_dir)
            for files in all_files:
                if files.endswith(".fastq"):
                    f_ = self.output_dir + '/' + files
                    if os.path.exists(f_):
                        os.remove(f_)
                    os.link(self.work_dir + "/" + files, f_)
            self.logger.info("设置fastp分析结果目录成功")

        except Exception as e:
            self.logger.info("设置fastp分析结果目录失败{}".format(e))
            self.set_error("设置fastp分析结果目录失败{}".format(e))
