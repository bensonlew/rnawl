# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class AssembleAmphoraModule(Module):
    """
    组装结果注释
    """
    def __init__(self, work_id):
        super(AssembleAmphoraModule, self).__init__(work_id)
        options = [
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  # 最佳的scaffold文件
            {"name": "taxon", "type": "string", "default": "Bacteria"}, #输入的细菌还是古菌
            {'name': 'genome_id', "type": "string"},  # 基因组名称
            {"name": "anno", "type": "outfile", "format": "metagbin.file_table"}, #基因组物种注释结果
        ]
        self.amphora_anno = self.add_tool('metagbin.amphora_anno')
        self.amphora = self.add_tool('metagbin.amphora')
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('genome').is_set:
            raise OptionError('必须输入组装序列')
        if not self.option('taxon'):
            raise OptionError('必须输入类型是古菌还是细菌！')

    def run_amphora(self):
        opts = {
            "bin_fa": self.option('genome'),
            "kingdom": self.option('taxon'),
            "bin_name": self.option('genome_id'),
        }
        self.amphora.set_options(opts)
        self.amphora.on('end', self.run_amphora_summ)
        self.amphora.run()

    def run_amphora_summ(self):
        opts = {
            "amphora_dir": self.amphora.output_dir,
        }
        self.amphora_anno.set_options(opts)
        self.amphora_anno.on("end",self.set_output)
        self.amphora_anno.run()

    def run(self):
        """
        运行
        :return:
        """
        super(AssembleAmphoraModule, self).run()
        self.run_amphora()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + "/assembly.anno.xls"):
            os.remove(self.output_dir + "/assembly.anno.xls")
        os.link(self.amphora_anno.output_dir + '/summary.anno.xls', self.output_dir + "/assembly.anno.xls")
        self.option('anno', self.output_dir + "/assembly.anno.xls")
        self.logger.info('物种注释结果生成成功')
        self.end()

    def end(self):
        super(AssembleAmphoraModule, self).end()