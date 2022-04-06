# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os,glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest


class ExpUnitsConvertAgent(Agent):
    """
    Used for convert RNA-seq expression units convert, count/CPM/RPM/TPM/FPKM/RPKM are implemented.

    RPM/CPM: Reads/Counts of exon model per million mapped reads
    RPM/CPM=Total exon reads/ Mapped reads(Millions)

    RPKM/FPKM: Reads/Fragments Per Kilobase of exon model per Million mapped reads
    RPKM/FPKM=Total exon reads/[Mapped reads(Millions)*Exon length(Kb)]

    TPM is like RPKM and FPKM, except the order of operation is switched.
    """
    def __init__(self, parent):
        super(ExpUnitsConvertAgent, self).__init__(parent)
        options = [
            {"name": "exp_matrix", "type": "infile", "format": "ref_rna_v2.common"},  # 表达量矩阵
            {"name": "gene_length", "type": "infile", "format": "ref_rna_v2.common"},  # 基因长度文件
            {"name": "convert_type", "type": "string", "default": ''},  # 转换类型
            {"name": "float_num", "type": "int", "default": '4'},  # 保留小数位数
            {"name": "intersect", "type": "bool", "default": True},  # 是否取交集
        ]
        self.add_option(options)
        self.step.add_steps("exp_units_convert")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.exp_units_convert.start()
        self.step.update()

    def stepfinish(self):
        self.step.exp_units_convert.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('exp_matrix').is_set:
            raise OptionError('必须输入表达定量文件')
        if not self.option('convert_type') in ["count2tpm", "count2cpm", "count2fpkm", "fpkm2tpm", "cpm2tpm", "cpm2fpkm", "fpkm2cpm"]:
            raise OptionError('不支持该定量指标转换类型')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(ExpUnitsConvertAgent, self).end()


class ExpUnitsConvertTool(Tool):
    def __init__(self, config):
        super(ExpUnitsConvertTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'exp_units_convert': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/exp_units_convert.py')
        }

    def run(self):
        """
        运行
        :return:
        """
        super(ExpUnitsConvertTool, self).run()
        self.run_exp_units_convert()
        self.set_output()
        self.end()

    def run_exp_units_convert(self):
        cmd = '{} {}'.format(self.program['python'], self.script['exp_units_convert'])
        cmd += ' -exp_matrix {}'.format(self.option('exp_matrix').prop["path"])
        cmd += ' -convert_type {}'.format(self.option('convert_type'))
        if self.option("convert_type") in ['count2tpm', 'count2fpkm', 'cpm2fpkm', 'cpm2tpm', 'fpkm2cpm']:
            cmd += ' -gene_length {}'.format(self.option('gene_length').prop["path"])
            cmd += ' -intersect {}'.format(self.option('intersect'))
        cmd += ' -num {}'.format(self.option('float_num'))
        cmd_name = 'run_exp_units_convert'
        runcmd(self, cmd_name, cmd)


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            convert_output = glob.glob("*" + self.option("convert_type") + ".xls")[0]
            os.link(os.path.join(self.work_dir, convert_output), os.path.join(self.output_dir, convert_output))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "exp_units_convert_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.exp_units_convert",
            "options": dict(
                exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/known_seqs_count.matrix",
                convert_type="count2tpm",
                gene_length='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/gene_length.txt',
                float_num=4,
                intersect=True,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
