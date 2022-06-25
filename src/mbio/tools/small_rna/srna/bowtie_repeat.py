# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import glob
import ConfigParser
import unittest
import subprocess


class BowtieRepeatAgent(Agent):
    """
    sRNA种类鉴定
    """
    def __init__(self, parent):
        super(BowtieRepeatAgent, self).__init__(parent)
        options = [
            {"name": "file_type", "type": "string", "default": ""}, # 输入序列类型，fastq OR fasta
            {"name": "seq_type", "type": "string", "default": ""}, # 输入测序类型，PE OR SE
            {"name": "pe1", "type": "string", "default": ""}, # PE1的数据路径
            {"name": "pe2", "type": "string", "default": ""}, # PE2的数据路径
            {"name": "se", "type": "string", "default": ""}, # SE的数据路径
            {"name": "repeat", "type": "string", "default": ""}, # 参考基因组重复区域序列
            {"name": "mismatch", "type": "int", "default": 0}, # max # mismatches in seed alignment; can be 0 or 1
            {"name": "software", "type": "string", "default": "bowtie2"}, # 选取比对软件，bowtie or bowtie2
        ]
        self.add_option(options)
        self.step.add_steps("bowtie_repeat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.bowtie_repeat.start()
        self.step.update()

    def stepfinish(self):
        self.step.bowtie_repeat.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if self.option("seq_type").lower() == "pe":
            if self.option("pe1") == "" or self.option("pe2") == "":
                raise OptionError("测序数据为PE时，必须提供双端测序数据路径")
        if self.option("file_type").lower() not in ["fasta", "fastq"]:
            raise OptionError("输入的序列必须为fastq或者fasta")
        if self.option("mismatch") not in [0, 1]:
            raise OptionError("允许的错误匹配数个只能为1或0")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 11
        self._memory = '50G'

    def end(self):
        super(BowtieRepeatAgent, self).end()


class BowtieRepeatTool(Tool):
    def __init__(self, config):
        super(BowtieRepeatTool, self).__init__(config)
        self.bowtie2_path = "bioinfo/align/bowtie2-2.2.9"
        self.bowtie_path = self.config.SOFTWARE_DIR + "/bioinfo/align/bowtie-1.1.2"
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.set_environ(PATH=self.bowtie_path)
        self.set_environ(PATH=self.bowtie2_path)

    def run(self):
        """
        运行
        :return:
        """
        super(BowtieRepeatTool, self).run()
        self.build_index()
        self.bowtie_repeat()
        self.set_output()
        self.end()

    def build_index(self):
        repeat = self.option("repeat")
        if self.option("software").lower() == "bowtie":
            repeat_index = glob.glob(repeat + "*ebwt")
            if len(repeat_index) != 6:
                if os.path.getsize(repeat) >= 40000000000:
                    cmd = "%s/bowtie-build --large-index  %s %s" %("/bioinfo/align/bowtie-1.1.2", repeat, repeat)
                else:
                    cmd = "%s/bowtie-build %s %s" %("/bioinfo/align/bowtie-1.1.2", repeat, repeat)
                command = self.add_command("bowtie_build", cmd).run()
                self.wait()
                if command.return_code == 0:
                    self.logger.info("bowtie索引建立完成")
                else:
                    self.set_error("bowtie索引建立出错")
            else:
                pass
        else:
            repeat_index = glob.glob(repeat + "*bt2")
            if len(repeat_index) != 6:
                cmd = "%s/bowtie2-build --threads 10 %s %s" %(self.bowtie2_path, repeat, repeat)
                command = self.add_command("bowtie2_build", cmd).run()
                self.wait()
                if command.return_code == 0:
                    self.logger.info("bowtie2索引建立完成")
                else:
                    self.set_error("bowtie2索引建立出错")
            else:
                pass

    def bowtie_repeat(self):
        """
        repeat比对
        """
        self.logger.info("开始进行repeat库比对")
        repeat = self.option("repeat")
        output = os.path.join(self.work_dir, "bowtie_repeat.sam")
        if self.option("software").lower() == "bowtie2":
            if self.option("file_type").lower() == "fasta":
                if self.option("seq_type").lower() == "se":
                    if os.path.getsize(repeat) >= 40000000000:
                        cmd = "%s/bowtie2 -p 10 -k 1 -f -N %s -U %s -x %s -S %s --no-unal --no-head --large-index" % (self.bowtie2_path, self.option("mismatch"), self.option("se"), self.option("repeat"), output)
                    else:
                        cmd = "%s/bowtie2 -p 10 -k 1 -f -N %s -U %s -x %s -S %s --no-unal --no-head" % (self.bowtie2_path, self.option("mismatch"), self.option("se"), self.option("repeat"), output)
                else:
                    if os.path.getsize(repeat) >= 40000000000:
                        cmd = "%s/bowtie2 -p 10 -k 1 -f -N %s -1 %s -2 %s -x %s -S %s --no-unal --no-head --large-index" % (self.bowtie2_path, self.option("mismatch"), self.option("pe1"), self.option("pe2"), self.option("repeat"), output)
                    else:
                        cmd = "%s/bowtie2 -p 10 -k 1 -f -N %s -1 %s -2 %s -x %s -S %s --no-unal --no-head" % (self.bowtie2_path, self.option("mismatch"), self.option("pe1"), self.option("pe2"), self.option("repeat"), output)
            else:
                if self.option("seq_type").lower() == "se":
                    if os.path.getsize(repeat) >= 40000000000:
                        cmd = "%s/bowtie2 -p 10 -k 1 -q -N %s -U %s -x %s -S %s --no-unal --no-head --large-index" % (self.bowtie2_path, self.option("mismatch"), self.option("se"), self.option("repeat"), output)
                    else:
                        cmd = "%s/bowtie2 -p 10 -k 1 -q -N %s -U %s -x %s -S %s --no-unal --no-head" % (self.bowtie2_path, self.option("mismatch"), self.option("se"), self.option("repeat"), output)
                else:
                    if os.path.getsize(repeat) >= 40000000000:
                        cmd = "%s/bowtie2 -p 10 -k 1 -q -N %s -1 %s -2 %s -x %s -S %s --no-unal --no-head --large-index " % (self.bowtie2_path, self.option("mismatch"), self.option("pe1"), self.option("pe2"), self.option("repeat"), output)
                    else:
                        cmd = "%s/bowtie2 -p 10 -k 1 -q -N %s -1 %s -2 %s -x %s -S %s --no-unal --no-head" % (self.bowtie2_path, self.option("mismatch"), self.option("pe1"), self.option("pe2"), self.option("repeat"), output)
            command = self.add_command("bowtie2_repeat", cmd).run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("运行repeat库比对完成")
            else:
                self.set_error("运行repeat库比对出错")
        else:
            if self.option("file_type").lower() == "fasta":
                if self.option("seq_type").lower() == "se":
                    cmd = "%s/bowtie %s -p 10 -k 1 -f -v %s %s -S --best --sam-nohead > %s" % (self.bowtie_path, self.option("repeat"), self.option("mismatch"), self.option("se"), output)
                else:
                    cmd = "%s/bowtie %s -p 10 -k 1 -f -v %s -1 %s -2 %s -S --best --sam-nohead > %s" % (self.bowtie_path, self.option("repeat"), self.option("mismatch"), self.option("pe1"), self.option("pe2"), output)
            else:
                if self.option("seq_type").lower() == "se":
                    cmd = "%s/bowtie %s -p 10 -k 1 -q -v %s %s -S --best --sam-nohead > %s" % (self.bowtie_path, self.option("repeat"), self.option("mismatch"), self.option("se"), output)
                else:
                    cmd = "%s/bowtie %s -p 10 -k 1 -q -v %s -1 %s -2 %s -S --best --sam-nohead > %s" % (self.bowtie_path, self.option("repeat"), self.option("mismatch"), self.option("pe1"), self.option("pe2"), output)
            try:
                self.logger.info(cmd)
                subprocess.check_output(cmd, shell=True)
                self.logger.info('运行repeat库比对完成')
            except subprocess.CalledProcessError:
                self.set_error("运行repeat库比对出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        if os.path.exists(os.path.join(self.output_dir, "bowtie_repeat.sam")):
            os.remove(os.path.join(self.output_dir, "bowtie_repeat.sam"))
        os.link(os.path.join(self.work_dir, "bowtie_repeat.sam"), os.path.join(self.output_dir, "bowtie_repeat.sam"))
        self.logger.info("开始设置结果目录")
        self.logger.info("设置结果完成")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "BowtieRepeat" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.bowtie_repeat",
            "instant": False,
            "options": dict(
                file_type="fasta",
                seq_type="se",
                se="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/srna_test/mouse/input.fa",
                repeat="/mnt/ilustre/users/sanger-dev/workspace/20191114/Smallrna_tsg_36158/Srna/Blast/repeatmasker_merge.repeat.fa",
                software="bowtie",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()