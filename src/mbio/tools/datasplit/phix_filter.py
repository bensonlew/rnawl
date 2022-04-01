# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""PhixFilter 去污染 """
import os
import re
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class PhixFilterAgent(Agent):
    """
    PhixFilter
    """
    def __init__(self, parent=None):
        super(PhixFilterAgent, self).__init__(parent)
        options = [
            {'name': 'fq1', 'type': "infile", "format": "sequence.fastq"},  # 原始数据
            {'name': 'fq2', 'type': "infile", "format": "sequence.fastq"},  # 原始数据
            {'name': 'sample_name', "type": "string"},  # 样本名
            {'name': 'insert_size', "type": "string", "default": "500"},  # 插入片段长度
            {'name': 'flag', 'type': "string", "default": "4"},  # 提取没有比对上的reads,此处固定取值4，软件默认0
            {'name': 'readl', "type": "string"},  # 切除序列的阈值
            {'name': 'out_fq1', 'type': "outfile", "format": "sequence.fastq"},
            {'name': 'out_fq2', 'type': "outfile", "format": "sequence.fastq"},  # 去Phix之后的fq序列
        ]
        self.add_option(options)
        self.step.add_steps("phix_filter")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.phix_filter.start()
        self.step.update()

    def stepfinish(self):
        self.step.phix_filter.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option('fq1'):
            raise OptionError('必须输入1端序列')
        if not self.option('fq2'):
            raise OptionError('必须输入2端序列')
        if not self.option('sample_name'):
            raise OptionError('必须输入样本名')
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(PhixFilterAgent, self).end()


class PhixFilterTool(Tool):
    """
    """
    def __init__(self, config):
        super(PhixFilterTool, self).__init__(config)
        self.fq_list = []
        self.sai_list = []
        self.readl = 0
        self.unzip_path = '/bioinfo/seq/scripts/unzip.sh '
        self.bwa = '/bioinfo/align/bwa-0.7.16/'
        self.bwa_path = self.config.SOFTWARE_DIR + '/bioinfo/align/bwa-0.7.16/'
        self.bwa_aln_path = '/bioinfo/align/scripts/bwa_aln.sh '
        self.database = self.config.SOFTWARE_DIR + '/database/datasplit/phix/Phix '
        self.samtools_path = '/bioinfo/align/samtools-1.4/bin/samtools '
        self.bamtofastq_path = '/bioinfo/seq/bedtools-2.25.0/bin/bamToFastq '
        self.stat_path = '/bioinfo/seq/fastx_toolkit_0.0.14/fastx_quality_stats '
        self.fastq_cut = '/bioinfo/seq/scripts/fastq_cut.sh'
        self.fastq_cut_path = self.config.SOFTWARE_DIR + '/bioinfo/seq/scripts/'

    def run(self):
        super(PhixFilterTool, self).run()
        self.run_gunzip()
        self.run_bwa_aln()
        self.run_get_sam()
        self.run_samtools()
        self.run_bamtofastq()
        self.logger.info('对序列进行质量统计，得到参数readI')
        self.run_fq_stat()
        self.set_output()
        self.end()

    def run_gunzip(self):
        """
        将两序列解压
        """
        self.logger.info('对压缩文件进行解压')
        fq_list = [self.option('fq1').prop['path'], self.option('fq2').prop['path']]
        cmd_list = []
        for i in range(len(fq_list)):
            cmd = self.unzip_path + '%s %s' % (fq_list[i], self.work_dir + '/' + self.option('sample_name') + '.' + str(i+1) + '.fq')
            cmd_list.append(cmd)
        for i in range(len(cmd_list)):
            command = self.add_command("unzip_cmd{}".format(i), cmd_list[i]).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("unzip{}运行完成".format(i))
            else:
                for i in range(len(cmd_list)):
                    command = self.add_command("rerun_unzip_cmd{}".format(i), cmd_list[i]).run()
                    self.wait(command)
                    if command.return_code == 0:
                        self.logger.info("unzip{}运行完成".format(i))
                    else:
                        self.set_error("unzip{}运行出错!".format(i))

    def run_bwa_aln(self):
        """
        原始reads mapping，生成*sai文件
        :return:
        """
        self.logger.info('原始reads mapping，生成sai文件')
        cmd_list = []
        all_files = os.listdir(self.work_dir)
        self.fq_list = [self.work_dir + '/' + self.option('sample_name') + '.1.fq', self.work_dir + '/' + self.option('sample_name') + '.2.fq']
        for fq in self.fq_list:
            sample_name = os.path.basename(fq).strip().split('.fq')[0]
            cmd = self.bwa_aln_path + '%s %s %s %s' % (self.bwa_path, self.database, fq, self.work_dir + '/' + sample_name + '.sai')
            cmd_list.append(cmd)
        for i in range(len(cmd_list)):
            command = self.add_command("bwa_aln_cmd{}".format(i), cmd_list[i]).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("bwa_aln{}运行完成".format(i))
            else:
                self.set_error("bwa_aln{}运行出错!".format(i))

    def run_get_sam(self):
        """
        得到sam文件
        """
        self.sai_list = [self.work_dir + '/' + self.option('sample_name') + '.1.sai', self.work_dir + '/' + self.option('sample_name') + '.2.sai']
        cmd = self.bwa + 'bwa sampe -a %s -f %s %s %s %s %s %s' % (
            self.option('insert_size'), self.work_dir + '/' + self.option('sample_name') + '.sam', self.database,
            self.sai_list[0], self.sai_list[1], self.fq_list[0], self.fq_list[1])
        self.logger.info('运行bwa，得到sam文件')
        command = self.add_command("bwa_sampe", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bwa运行完成")
        else:
            command = self.add_command("rerun_bwa_sampe", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("bwa运行完成")
            else:
                self.set_error("bwa运行出错!")

    def run_samtools(self):
        """
        得到bam文件
        """
        cmd = self.samtools_path + 'view -S %s -f %s -b -o %s' % (
            self.work_dir + '/' + self.option('sample_name') + '.sam', self.option('flag'),
            self.work_dir + '/' + self.option('sample_name') + '.unmapping_phix.bam')
        self.logger.info('运行samtools，得到bam文件')
        command = self.add_command("samtools_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("samtools运行完成")
        else:
            self.set_error("samtools运行出错!")

    def run_bamtofastq(self):
        """
        得到去phix后的fastq文件
        """
        cmd = self.bamtofastq_path + '-i %s -fq %s -fq2 %s' % (
            self.work_dir + '/' + self.option('sample_name') + '.unmapping_phix.bam', self.work_dir + '/' + self.option('sample_name') + '.filtphix.1.fq',
            self.work_dir + '/' + self.option('sample_name') + '.filtphix.2.fq')
        self.logger.info('运行BamToFastq，得到fq文件')
        command = self.add_command("bamtofastq_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bamtofastq运行完成")
        else:
            self.set_error("bamtofastq运行出错!")

    def run_fq_stat(self):
        """
        对fq序列进行质量统计，根据统计结果判断是否需要进行剪切
        :return:
        """
        cmd_list = []
        # filt_fq = [self.option('sample_name') + '.filtphix.1.fq', self.option('sample_name') + '.filtphix.2.fq']
        filt_fq = [self.option('sample_name') + '.1.fq']
        for fq in filt_fq:
            sample_name = os.path.basename(fq)
            cmd = self.stat_path + '-i %s -Q 33 -o %s' % (fq, self.work_dir + '/' + sample_name + '.stat')
            cmd_list.append(cmd)
        for i in range(len(cmd_list)):
            command = self.add_command("fq_stat{}".format(i), cmd_list[i]).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("fq_stat{}运行完成".format(i))
            else:
                self.set_error("fq_stat{}运行出错!".format(i))
        # cmd2 = 'wc -l %s|awk -F \" \" \'{print $1-1}\'' % (self.work_dir + '/' + self.option('sample_name') + '.filtphix.1.fq.stat')
        cmd2 = 'wc -l %s|awk -F \" \" \'{print $1-1}\'' % (self.work_dir + '/' + self.option('sample_name') + '.1.fq.stat')
        readl = os.popen(cmd2).readlines()[0].strip()
        if self.option("readl"):
            self.readl = int(self.option("readl"))
        else:
            if int(readl) == 1:
                self.readl = 0
            else:
                self.readl = int(readl)
        if self.readl >= 151:
            self.fastq_cut_run()

    def fastq_cut_run(self):
        cmd_list = []
        filt_fq = [self.option('sample_name') + '.filtphix.1.fq', self.option('sample_name') + '.filtphix.2.fq']
        for fq in filt_fq:
            sample_name = os.path.basename(fq).strip().split('.fq')[0]
            cmd = self.fastq_cut + ' %s %s %s %s' % (self.fastq_cut_path, fq, str(self.readl), self.work_dir + '/' + sample_name + '.cut.fq')
            cmd_list.append(cmd)
        for i in range(len(cmd_list)):
            command = self.add_command("fastq_cut{}".format(i), cmd_list[i]).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("fastq_cut{}运行完成".format(i))
            else:
                self.set_error("fastq_cut{}运行出错!".format(i))

    def linkdir(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            all_files = os.listdir(self.work_dir)
            for files in all_files:
                if self.readl >= 151:
                    if files.endswith(".cut.fq"):
                        base_name = os.path.basename(files).strip().split('.cut.fq')[0]
                        self.linkdir(self.work_dir + "/" + files, self.output_dir + '/' + base_name + '.result.fq')
                else:
                    if files.endswith(".filtphix.1.fq"):
                        self.linkdir(self.work_dir + "/" + files, self.output_dir + '/' + self.option('sample_name') + '.filtphix.1.result.fq')
                    if files.endswith(".filtphix.2.fq"):
                        self.linkdir(self.work_dir + "/" + files, self.output_dir + '/' + self.option('sample_name') + '.filtphix.2.result.fq')
            self.option('out_fq1').set_path(self.output_dir + '/' + self.option('sample_name') + ".filtphix.1.result.fq")
            self.option('out_fq2').set_path(self.output_dir + '/' + self.option('sample_name') + ".filtphix.2.result.fq")
            self.logger.info("设置PhixFilter分析结果目录成功")

        except Exception as e:
            self.logger.info("设置PhixFilter分析结果目录失败{}".format(e))
            self.set_error("设置PhixFilter分析结果目录失败{}".format(e))
