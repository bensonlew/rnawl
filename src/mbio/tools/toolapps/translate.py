# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.toolapps.common import link_dir


class TranslateAgent(Agent):
    """
    宏转录组 运行transeq软件进行翻译蛋白
    配置参数和资源
    """
    def __init__(self, parent):
        super(TranslateAgent, self).__init__(parent)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 输入翻译文件夹
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},  # 输入翻译文件
            {"name": "is_dir", "type": "string", "default": "false"}, # 输入是否是文件夹
            {"name": "table", "type": "int", "default": 11 }, #翻译输出表类型
        ]
        self.add_option(options)
        self.step.add_steps('translate')
        self._memory_increase_step = 50

    def check_options(self):
        """
        参数二次检查
        """
        if (not self.option("fasta_dir").is_set) and (not self.option("fasta").is_set):
            raise OptionError("必须输入需要翻译的序列文件或者文件夹", code="34400701")
        if not self.option("is_dir") in ['false', 'true']:
            raise OptionError("必须输入要翻译的序列文件是否是文件夹", code="34400702")
        if self.option("is_dir") in ['true']:
            if not self.option("fasta_dir").is_set:
                raise OptionError("必须输入需要翻译的序列文件夹", code="34400703")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        number = float(os.path.getsize(self.option('fasta').prop['path'])) / 1024*1024*500
        if number < 20:
            new_number = 20
        else:
            new_number = number
        self._memory = '{}G'.format(int(new_number))

    def end(self):
        """
        运行结束
        """
        super(TranslateAgent, self).end()


class TranslateTool(Tool):
    """
    tool运行模块
    """
    def __init__(self, config):
        super(TranslateTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = '/program/perl-5.24.0/bin/perl '
        self.trans_path = 'bioinfo/seq/emboss6.6/bin/transeq'
        if not os.path.isfile(os.path.join(self.config.SOFTWARE_DIR, self.trans_path)):
            self.trans_path = 'bioinfo/seq/EMBOSS-6.6.0/emboss/transeq'
        if self.option("fasta").is_set:
            self.fasta = self.option("fasta").prop['path']
        if self.option("fasta_dir").is_set:
            self.fasta = self.option("fasta_dir").prop['path']

    def run(self):
        """
        运行
        """
        super(TranslateTool, self).run()
        self.logger.info("开始运行 tool！")
        if self.option("is_dir") in ['false']:
            self.run_translate()
        else:
            self.run_translate_dir()
        self.set_output()
        self.end()

    def run_translate(self):
        """
        对核酸序列进行翻译
        注意：这里如果进行翻译，那么必然是通过核酸进行去冗余
        """
        self.logger.info("开始对核酸文件进行翻译")
        real_fastaa = os.path.join(self.work_dir, 'gene.uniGeneset.faa')
        cmd = '%s -sequence %s -table %s -trim -outseq %s' % (
            self.trans_path, self.fasta, self.option("table"), real_fastaa)
        self.logger.info(cmd)
        command = self.add_command('translate', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("translate success!")
        else:
            self.set_error("translate failed", code="34400701")

    def run_translate_dir(self):
        """
        对核酸序列进行翻译
        注意：这里如果进行翻译，那么必然是通过核酸进行去冗余
        """
        self.logger.info("开始对核酸文件进行翻译")
        all_list = os.listdir(self.fasta)
        output_dir = os.path.join(self.work_dir, "translate")
        i = 1
        for file in all_list:
            fasta_path = os.path.join(self.fasta, file)
            newname = file.split(".fa")[0]
            newname_path = os.path.join(output_dir, newname + ".faa")
            cmd = '%s -sequence %s -table %s -trim -outseq %s' % (
                self.trans_path, fasta_path, self.option("table"), newname_path)
            self.logger.info(cmd)
            command = self.add_command('translate%s'%i, cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                i += 1
                self.logger.info("translate success!")
            else:
                self.set_error("translate failed", code="34400702")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
            os.mkdir(self.output_dir)
        if self.option("is_dir") in ['false']:
            if os.path.exists(os.path.join(self.output_dir, "gene.uniGeneset.faa")):
                os.remove(os.path.join(self.output_dir, "gene.uniGeneset.faa"))
            else:
                os.link(os.path.join(self.work_dir, "gene.uniGeneset.faa"), os.path.join(self.output_dir, "gene.uniGeneset.faa"))
        else:
            link_dir(os.path.join(self.work_dir, "translate"), self.output_dir)
        self.logger.info("设置结果文件目录成功")


