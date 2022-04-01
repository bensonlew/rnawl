# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""多样性去接头 """
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os


class GetvaildBybarcodeAgent(Agent):
    """
    GetvaildBybarcode
    """
    def __init__(self, parent=None):
        super(GetvaildBybarcodeAgent, self).__init__(parent)
        options = [
            {'name': 'fq1', 'type': "infile", "format": "datasplit.fastq"},  # 一次拆分后的文库序列
            {'name': 'fq2', 'type': "infile", "format": "datasplit.fastq"},
            {'name': "barcode_info", "type": "infile", "format": "datasplit.meta_barcode_info"},  # 样本总信息
            {'name': 'out_fq1', 'type': "outfile", "format": "datasplit.fastq"},
            {'name': 'out_fq2', 'type': "outfile", "format": "datasplit.fastq"},  # 去嵌合体之后的fq序列
            {'name': 'lib_name', 'type': "string"},  # 文库名称
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        """
        参数检测
        """
        if not self.option('fq1'):
            raise OptionError('必须输入1端序列')
        if not self.option('fq2'):
            raise OptionError('必须输入2端序列')
        if not self.option('barcode_info'):
            raise OptionError('必须输入barcode信息表')
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"],
        # ])
        super(GetvaildBybarcodeAgent, self).end()


class GetvaildBybarcodeTool(Tool):
    def __init__(self, config):
        super(GetvaildBybarcodeTool, self).__init__(config)
        self.perl = "program/perl-5.24.0/bin/perl "
        # self.getvaild_bybarcode_path = '/bioinfo/seq/scripts/getValid_rawSeq_byBarcode.2.pl '
        self.getvaild_bybarcode_path = self.config.PACKAGE_DIR + '/datasplit/getValid_rawSeq_byBarcode.2.pl '

    def run(self):
        super(GetvaildBybarcodeTool, self).run()
        self.run_getvaild_bybarcode()
        self.set_output()
        self.end()

    def run_getvaild_bybarcode(self):
        """
        运行GetvaildBybarcode,进行质控
        """
        # cmd = self.getvaild_bybarcode_path + ' %s %s %s %s' % (
        cmd = self.perl + self.getvaild_bybarcode_path + ' %s %s %s %s' % (
            self.option('fq1').prop['path'], self.option('fq2').prop['path'], self.option('barcode_info').prop['path'],
            self.work_dir + "/" + self.option('lib_name') + ".all.raw")
        self.logger.info('运行getvaild_bybarcode，进行去低值')
        command = self.add_command("getvaild_bybarcode_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("getvaild_bybarcode运行完成")
        else:
            self.set_error("getvaild_bybarcode运行出错!")

    def set_output(self):
        self.logger.info("设置结果目录")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        os.link(os.path.join(self.work_dir, self.option('lib_name') + ".all.raw.valid.1.fq"),
                os.path.join(self.output_dir, self.option('lib_name') + ".all.raw.valid.1.fq"))
        os.link(os.path.join(self.work_dir, self.option('lib_name') + ".all.raw.valid.2.fq"),
                os.path.join(self.output_dir, self.option('lib_name') + ".all.raw.valid.2.fq"))
        # os.link(os.path.join(self.work_dir, self.option('lib_name') + ".all.raw.seq2sam.stat"),
        #         os.path.join(self.output_dir, self.option('lib_name') + ".all.raw.seq2sam.stat"))
        self.option('out_fq1').set_path(self.output_dir + '/' + self.option('lib_name') + ".all.raw.valid.1.fq")
        self.option('out_fq2').set_path(self.output_dir + '/' + self.option('lib_name') + ".all.raw.valid.2.fq")
        self.logger.info("设置getvaild_bybarcode分析结果目录成功")
