# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""fastp质控工具,专用于meta，避免因版本不同导致的个别样本reads数不同 """
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class FastpMetaAgent(Agent):
    """
    fastp
    """
    def __init__(self, parent=None):
        super(FastpMetaAgent, self).__init__(parent)
        options = [
            {'name': 'fq1', 'type': "infile", "format": "datasplit.fastq"},  # 原始数据
            {'name': 'fq2', 'type': "infile", "format": "datasplit.fastq"},  # 原始数据
            {'name': 'length_required', "type": "string", "default": "50"},  # -l,长度过滤参数，比此值短的读取将被丢弃
            {'name': "cut_by_quality5", "type": "string", "default": "0"},  # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {'name': "cut_by_quality3", "type": "string", "default": "20"},  # -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            {"name": "cut_right_mean_quality", "type": "string", "default": "20"},  # cut_right的平均质量要求,默认20
            {"name": "cut_right_window_size", "type": "string", "default": "50"},  # cut_right的窗口大小，默认4
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option('fq1').is_set:
            raise OptionError('必须输入1端序列')
        if not self.option('fq2').is_set:
            raise OptionError('必须输入2端序列')
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(FastpMetaAgent, self).end()


class FastpMetaTool(Tool):
    def __init__(self, config):
        super(FastpMetaTool, self).__init__(config)
        self.fastp_path = 'bioinfo/seq/fastp-0.19.6 '

    def run(self):
        super(FastpMetaTool, self).run()
        self.run_fastp()
        self.set_output()
        self.end()

    def run_fastp(self):
        """
        运行fastp,进行质控
        """
        sample_name = os.path.basename(self.option('fq1').prop['path']).split('.all.raw.valid.1.fq')[0]
        cmd = self.fastp_path + ' -i %s -I %s -o %s -O %s -5 %s -3 %s --cut_right_window_size %s' % (
              self.option('fq1').prop['path'], self.option('fq2').prop['path'],
              os.path.join(self.work_dir, sample_name + '.trim.1.fq'), os.path.join(self.work_dir, sample_name + '.trim.2.fq'),
              self.option('cut_by_quality5'), self.option('cut_by_quality3'), self.option('cut_right_window_size'))
        cmd += " --cut_right_mean_quality %s --length_required %s --json %s --html %s" % (self.option('cut_right_mean_quality'),
              self.option('length_required'), os.path.join(self.work_dir, sample_name + '.json'), os.path.join(self.work_dir, sample_name + '.html'))
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
        all_files = os .listdir(self.work_dir)
        for files in all_files:
            if files.endswith(".fq") or files.endswith(".json"):
                f_ = self.output_dir + '/' + files
                if os.path.exists(f_):
                    os.remove(f_)
                os.link(self.work_dir + "/" + files, f_)
