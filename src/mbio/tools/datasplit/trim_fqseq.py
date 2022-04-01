# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""长度过滤工具 """
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class TrimFqseqAgent(Agent):
    """
    TrimFqseq
    """
    def __init__(self, parent=None):
        super(TrimFqseqAgent, self).__init__(parent)
        options = [
            {'name': 'fq', 'type': "infile", "format": "sequence.fastq"},
            {'name': 'start_pos', 'type': "string", "default": "1"},  # -s,开始位点
            {'name': 'valid_len', "type": "string"},  # -l,长度过滤阈值
            {'name': "min_len", "type": "string", "default": "0"},  # -m,最小长度
            # {'name': "max_len", "type": "string"},  # -x,最大长度
            {'name': 'trim_end_len', "type": "string", "default": "0"},  # -e,在末端修剪长度序列
            {'name': 'lib_name', 'type': "string"},  # 文库名称
            {'name': 'out_fq', 'type': "outfile", "format": "sequence.fastq"},  # 长度过滤之后的序列
        ]
        self.add_option(options)
        self.step.add_steps("trim_fqseq")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.trim_fqseq.start()
        self.step.update()

    def stepfinish(self):
        self.step.trim_fqseq.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option('fq'):
            raise OptionError('必须输入fq序列')
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
        super(TrimFqseqAgent, self).end()


class TrimFqseqTool(Tool):
    """
    """
    def __init__(self, config):
        super(TrimFqseqTool, self).__init__(config)
        self.trim_fqseq_path = '/bioinfo/seq/scripts/trim_fqSeq.pl '

    def run(self):
        super(TrimFqseqTool, self).run()
        self.run_trim_fqseq()
        self.set_output()
        self.end()

    def run_trim_fqseq(self):
        """
        运行TrimFqseq,进行质控
        """
        if self.option("valid_len"):
            cmd = self.trim_fqseq_path + ' -i %s -o %s -s %s -l %s -m %s  -e %s' % (
                self.option('fq').prop['path'], self.work_dir + "/" + self.option("lib_name") + '.fq',
                self.option("start_pos"), self.option("valid_len"), self.option("min_len"),
                self.option("trim_end_len"))
        else:
            cmd = self.trim_fqseq_path + ' -i %s -o %s -s %s -m %s -e %s' % (
                self.option('fq').prop['path'], self.work_dir + "/" + self.option("lib_name") + '.fq',
                self.option("start_pos"), self.option("min_len"),
                self.option("trim_end_len"))
        self.logger.info('运行trim_fqseq.pl，进行长度过滤')
        command = self.add_command("trim_fqseq_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("trim_fqseq运行完成")
        else:
            self.set_error("trim_fqseq运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            f_ = self.output_dir + '/' + self.option("lib_name") + ".fq"
            if os.path.exists(f_):
                os.remove(f_)
            os.link(self.work_dir + "/" + self.option("lib_name") + ".fq", f_)
            self.option('out_fq').set_path(self.output_dir + '/' + self.option("lib_name") + ".fq")
            self.logger.info("设置trim_fqseq分析结果目录成功")

        except Exception as e:
            self.logger.info("设置trim_fqseq分析结果目录失败{}".format(e))
            self.set_error("设置trim_fqseq分析结果目录失败{}".format(e))
