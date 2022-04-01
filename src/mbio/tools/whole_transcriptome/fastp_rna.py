# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""fastp质控工具 """
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import unittest

class FastpRnaAgent(Agent):
    """
    fastp
    """
    def __init__(self, parent=None):
        super(FastpRnaAgent, self).__init__(parent)
        options = [
            {'name': 'fq1', 'type': "infile", "format": "sequence.fastq"},  # 原始数据
            {'name': 'fq2', 'type': "infile", "format": "sequence.fastq"},  # 原始数据
            {'name': 'qualified_quality_phred', 'type': "string", "default": "15"},  # -q,一个碱基合格的质量值,默认15表示phred质量> = Q15是合格的。
            {'name': 'length_required', "type": "string", "default": "30"},  # -l,长度过滤参数，比此值短的读取将被丢弃
            {'name': "cut_by_quality5", "type": "string"},  # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {'name': "cut_by_quality3", "type": "string"},  # -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            {'name': 'cut_mean_quality', "type": "string", "default": "20"},  # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {'name': 'n_base_limit', "type": "string", "default": "5"},  # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {'name': 'compression', "type": "string", "default": "2"},  # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {'name': 'thread', "type": "string", "default": "3"},  # -w,线程数
            {"name": "adapter_sequence", "type": "string", "default": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"},  # --adapter_sequence,the adapter for read1
            {"name": "adapter_sequence_r2", "type": "string", "default": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"},  # --adapter_sequence_r2,the adapter for read2 (PE data only)
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
        self._cpu = 6
        self._memory = '20G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(FastpRnaAgent, self).end()


class FastpRnaTool(Tool):
    """
    """
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
        # sample_name = os.path.basename(self.option('fq1').prop['path']).split('_R1.fastq.gz')[0]
        try:
            sample_name = '_'.join(os.path.basename(self.option('fq1').prop['path']).split('_')[0:-4])
        except:
            sample_name = os.path.basename(self.option('fq1').prop['path'])
        ## 增加判断，用于数据命名不规范是导致sample_name为空的情形
        if sample_name == "":
                sample_name = os.path.basename(self.option('fq1').prop['path']).split(".")[0]
        cmd = self.fastp_path + ' -i %s -I %s -o %s -O %s -q %s -l %s -M %s -n %s -j %s -h %s -w %s' % (
        self.option('fq1').prop['path'], self.option('fq2').prop['path'],
        self.work_dir + "/" + sample_name + '.clean.1.fastq', self.work_dir + "/" + sample_name + '.clean.2.fastq',
        self.option('qualified_quality_phred'),  self.option('length_required'),
        self.option('cut_mean_quality'),self.option('n_base_limit'),
        self.work_dir + "/" + sample_name + '.json', self.work_dir + "/" + sample_name + '.html',
        self.option('thread'))
        if self.option('cut_by_quality5'):
            cmd += ' -5 %s' % (self.option('cut_by_quality5'))
        if self.option('cut_by_quality3'):
            cmd += ' -3 %s' % (self.option('cut_by_quality3'))
        cmd += " --adapter_sequence {} --adapter_sequence_r2 {}".format(self.option("adapter_sequence"), self.option("adapter_sequence_r2"))
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
            all_files = os .listdir(self.work_dir)
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

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "whole_transcriptome_fastp" + str(random.randint(1, 10000)),
            "type": "tool",
            "instant": False,
            "name": "whole_transcriptome.fastp_rna",
            "options": dict(
                fq1 = '/mnt/ilustre/users/sanger-dev/sg-users/zhuwengen/A1_S5_L001_R1_001.fastq',
                fq2 = '/mnt/ilustre/users/sanger-dev/sg-users/zhuwengen/A1_S5_L001_R2_001.fastq'
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()


