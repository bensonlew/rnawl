# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""Flash拼接工具 """
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class FlashAgent(Agent):
    """
    Flash拼接
    """
    def __init__(self, parent=None):
        super(FlashAgent, self).__init__(parent)
        options = [
            {'name': 'fq1', 'type': "infile", "format": "datasplit.fastq"},
            {'name': 'fq2', 'type': "infile", "format": "datasplit.fastq"},
            {'name': 'min_lenth', 'type': "string", "default": "10"},  # -m,两个reads之间所需的最小重叠长度，以提供可靠的重叠
            {'name': 'max_lenth', "type": "string", "default": "70"},  # -M,两个reads之间的最大重叠长度
            {'name': "mismatch_rate", "type": "string", "default": "0.25"},  # -x,错配和重叠长度允许的最大比率
            {'name': "pred", "type": "string", "default": "33"},  # -p,FASTQ文件中的碱基的质量值，Pred33/Pred64.
            {'name': 'thread', "type": "string", "default": "1"},  # -t,线程数
            {'name': 'lib_name', 'type': "string"},  # 文库名称
            {'name': 'out_fq', 'type': "outfile", "format": "datasplit.fastq"},  # 拼接之后的序列
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option('fq1').is_set:
            raise OptionError('必须输入1端序列')
        if not self.option('fq2').is_set:
            raise OptionError('必须输入2端序列')
        if not self.option('lib_name'):
            raise OptionError('必须输入文库名称')

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(FlashAgent, self).end()


class FlashTool(Tool):
    def __init__(self, config):
        super(FlashTool, self).__init__(config)
        self.flash_path = '/bioinfo/seq/flash '

    def run(self):
        super(FlashTool, self).run()
        self.run_flash()
        self.set_output()
        self.end()

    def run_flash(self):
        """
        运行flash,进行拼接
        """
        cmd = self.flash_path + ' -m %s -M %s -x %s -p %s -t %s -d %s  -o %s %s %s' % (
            self.option('min_lenth'), self.option('max_lenth'), self.option('mismatch_rate'), self.option('pred'),
            self.option('thread'), self.work_dir, self.option('lib_name') + ".trim", self.option('fq1').prop['path'],
            self.option('fq2').prop['path'])
        self.logger.info('运行flash，进行质控')
        command = self.add_command("flash_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("flash运行完成")
        else:
            self.set_error("flash运行出错!")

    def set_output(self):
        fq_name = self.option("lib_name") + ".trim.extendedFrags.fastq"
        fq_path = os.path.join(self.output_dir, fq_name)
        if os.path.exists(fq_path):
            os.remove(fq_path)
        os.link(os.path.join(self.work_dir, fq_name), fq_path)
        self.option('out_fq').set_path(fq_path)
        hist_name = self.option("lib_name") + ".trim.hist"
        hist_path = os.path.join(self.output_dir, hist_name)
        if os.path.exists(hist_path):
            os.remove(hist_path)
        os.link(os.path.join(self.work_dir, hist_name), hist_path)
