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
            {'name': 'fq1', 'type': "infile", "format": "sequence.fastq"},
            {'name': 'fq2', 'type': "infile", "format": "sequence.fastq"},
            {'name': 'min_lenth', 'type': "string", "default": "10"},  # -m,两个reads之间所需的最小重叠长度，以提供可靠的重叠
            {'name': 'max_lenth', "type": "string", "default": "70"},  # -M,两个reads之间的最大重叠长度
            {'name': "mismatch_rate", "type": "string", "default": "0.25"},  # -x,错配和重叠长度允许的最大比率
            {'name': "pred", "type": "string", "default": "33"},  # -p,FASTQ文件中的碱基的质量值，Pred33/Pred64.
            {'name': 'thread', "type": "string", "default": "1"},  # -t,线程数
            {'name': 'lib_name', 'type': "string"},  # 文库名称
            {'name': 'out_fq', 'type': "outfile", "format": "sequence.fastq"},  # 拼接之后的序列
        ]
        self.add_option(options)
        self.step.add_steps("flash")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.flash.start()
        self.step.update()

    def stepfinish(self):
        self.step.flash.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option('fq1'):
            raise OptionError('必须输入1端序列')
        if not self.option('fq2'):
            raise OptionError('必须输入2端序列')
        if not self.option('lib_name'):
            raise OptionError('必须输入文库名称')
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
        super(FlashAgent, self).end()


class FlashTool(Tool):
    """
    """
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
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            os.link(self.work_dir + "/" + self.option('lib_name') + ".trim.extendedFrags.fastq",
                    self.output_dir + '/' + self.option('lib_name') + ".trim.extendedFrags.fastq")
            self.option('out_fq').set_path(self.output_dir + '/' + self.option('lib_name') + ".trim.extendedFrags.fastq")
            self.logger.info("设置flash分析结果目录成功")

        except Exception as e:
            self.logger.info("设置flash分析结果目录失败{}".format(e))
            self.set_error("设置flash分析结果目录失败{}".format(e))


