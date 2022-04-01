# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""长度过滤工具 """
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class TrimLengthAgent(Agent):
    """
    TrimFqseq
    """
    def __init__(self, parent=None):
        super(TrimLengthAgent, self).__init__(parent)
        options = [
            {'name': 'fq', 'type': "infile", "format": "datasplit.fastq"},
            {'name': 'start_pos', 'type': "string", "default": "1"},  # -s,开始位点
            {'name': 'valid_len', "type": "string"},  # -l,长度过滤阈值
            {'name': "min_len", "type": "string", "default": "0"},  # -m,最小长度
            # {'name': "max_len", "type": "string"},  # -x,最大长度
            {'name': 'trim_end_len', "type": "string", "default": "0"},  # -e,在末端修剪长度序列
            {'name': 'lib_name', 'type': "string"},  # 文库名称
            {'name': 'out_fq', 'type': "outfile", "format": "datasplit.fastq"},  # 长度过滤之后的序列
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option('fq'):
            raise OptionError('必须输入fq序列')
        if not self.option('lib_name'):
            raise OptionError('必须输入文库名称')
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(TrimLengthAgent, self).end()


class TrimLengthTool(Tool):
    def __init__(self, config):
        super(TrimLengthTool, self).__init__(config)
        self.perl = "program/perl-5.24.0/bin/perl "
        # self.trim_fqseq_path = '/bioinfo/seq/scripts/trim_fqSeq.pl '
        self.trim_fqseq_path = self.config.PACKAGE_DIR + '/datasplit/trim_fqSeq.pl '

    def run_trim_fqseq(self):
        """
        运行TrimFqseq,进行质控
        """
        # out_fq = self.work_dir + "/" + self.option("lib_name") + '.fq'
        out_fq = os.path.join(self.work_dir, os.path.basename(self.option('fq').prop['path']))
        if self.option("valid_len"):
            # cmd = self.trim_fqseq_path + '-i %s -o %s -s %s -l %s -m %s  -e %s' % (
            cmd = self.perl + self.trim_fqseq_path + '-i %s -o %s -s %s -l %s -m %s  -e %s' % (
                self.option('fq').prop['path'], out_fq,
                self.option("start_pos"), self.option("valid_len"), self.option("min_len"),
                self.option("trim_end_len"))
        else:
            # cmd = self.trim_fqseq_path + ' -i %s -o %s -s %s -m %s -e %s' % (
            cmd = self.perl + self.trim_fqseq_path + '-i %s -o %s -s %s -m %s -e %s' % (
                self.option('fq').prop['path'], out_fq,
                self.option("start_pos"), self.option("min_len"), self.option("trim_end_len"))
        command = self.add_command("trim_fqseq_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("trim_fqseq运行完成")
        else:
            self.set_error("trim_fqseq运行出错!")

    def set_output(self):
        self.logger.info("设置结果目录")
        # fq = os.path.join(self.output_dir, self.option("lib_name") + ".fq")
        fq = os.path.join(self.output_dir, os.path.basename(self.option('fq').prop['path']))
        if os.path.exists(fq):
            os.remove(fq)
        # os.link(os.path.join(self.work_dir, self.option("lib_name") + ".fq"), fq)
        fq1 = os.path.join(self.work_dir, os.path.basename(self.option('fq').prop['path']))
        os.link(fq1, fq)
        self.option('out_fq').set_path(fq)

    def run(self):
        super(TrimLengthTool, self).run()
        self.run_trim_fqseq()
        self.set_output()
        self.end()
