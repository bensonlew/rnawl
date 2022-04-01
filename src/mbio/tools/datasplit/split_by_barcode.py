# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""多样性根据barcode拆样本 """
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SplitByBarcodeAgent(Agent):
    """
    SplitByBarcode
    """
    def __init__(self, parent=None):
        super(SplitByBarcodeAgent, self).__init__(parent)
        options = [
            {'name': 'fq', 'type': "infile", "format": "sequence.fastq"},  # 需要拆样本的fq文件
            {'name': "barcode_info", "type": "infile", "format": "sequence.barcode_info"},  # barcode信息
            {'name': "mismatch", "type": "string", "default": "2"},  # -l，允许引物错配数
            {'name': 'split_type', 'type': "string", "default": "Pair"},  # 拆分样本序列类型 Pair or Single
            {'name': 'lib_name', 'type': "string"},  # 文库名称
            {'name': 'out_fq', 'type': "outfile", "format": "sequence.fastq"},  # 拆出样本补上样本之后的序列
        ]
        self.add_option(options)
        self.step.add_steps("split_by_barcode")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.split_by_barcode.start()
        self.step.update()

    def stepfinish(self):
        self.step.split_by_barcode.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option('fq'):
            raise OptionError('必须输入需要拆样本的序列')
        if not self.option('barcode_info'):
            raise OptionError('必须输入barcode信息表')
        if self.option('split_type') not in ["Pair", "Single"]:
            raise OptionError('拆分类型只能是Pair或Single')
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
        super(SplitByBarcodeAgent, self).end()


class SplitByBarcodeTool(Tool):
    """
    """
    def __init__(self, config):
        super(SplitByBarcodeTool, self).__init__(config)
        self.split_F_barcode_path = '/bioinfo/seq/scripts/split_byF-Barcode.pl '
        self.split_two_barcode_path = '/bioinfo/seq/scripts/split_byTwoBarcode.pl '

    def run(self):
        super(SplitByBarcodeTool, self).run()
        self.run_split_by_barcode()
        self.set_output()
        self.end()

    def run_split_by_barcode(self):
        """
        运行SplitByBarcode,拆样本
        """
        if self.option('split_type') == 'Pair':
            cmd = self.split_two_barcode_path + ' -i %s -f %s -l %s -o %s' % (
                self.option('fq').prop['path'], self.option('barcode_info').prop['path'], self.option("mismatch"),
                self.work_dir + "/" + self.option('lib_name') + ".trim.merge")
        else:
            cmd = self.split_F_barcode_path + ' -i %s -f %s -l %s -o %s' % (
                self.option('fq').prop['path'], self.option('barcode_info').prop['path'], self.option("mismatch"),
                self.work_dir + "/" + self.option('lib_name') + ".trim")
        self.logger.info('运行split_by_barcode，进行去低值')
        command = self.add_command("split_by_barcode_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("split_by_barcode运行完成")
        else:
            self.set_error("split_by_barcode运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            if self.option('split_type') == 'Pair':
                os.link(self.work_dir + "/" + self.option('lib_name') + ".trim.merge.split.fq",
                        self.output_dir + '/' + self.option('lib_name') + ".trim.merge.split.allLen.fq")
                self.option('out_fq').set_path(self.output_dir + '/' + self.option('lib_name') + ".trim.merge.split.allLen.fq")
            else:
                os.link(self.work_dir + "/" + self.option('lib_name') + ".trim.split.fq",
                        self.output_dir + '/' + self.option('lib_name') + ".trim.split.allLen.fq")
                self.option('out_fq').set_path(
                    self.output_dir + '/' + self.option('lib_name') + ".trim.split.allLen.fq")
            self.logger.info("设置split_by_barcode分析结果目录成功")
        except Exception as e:
            self.logger.info("设置split_by_barcode分析结果目录失败{}".format(e))
            self.set_error("设置split_by_barcode分析结果目录失败{}".format(e))


