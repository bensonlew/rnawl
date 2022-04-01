# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import os, sys
# import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class BacDenovaAssAgent(Agent):
    """
    宏基因SOAPdenovo2组装细菌基因组
    version: SOAPdenovo2-src-r240
    author: hao.gao
    last_modify: 2017.01.26
    """

    def __init__(self, parent):
        super(BacDenovaAssAgent, self).__init__(parent)
        options = [
            {"name": "config", "type": "infile", "format": "bacgenome.config_file"},  # 输入文件,sample.sickle.l.fastq
            {"name": "kmerFreqCutoff", "type": "int"},  # 输入文件
            {"name": "sample_name", "type": "string"},  # 样品名称
            {"name": "kmer", "type": "int"},  # k_mer值，例"39"
        ]
        self.add_option(options)
        self.step.add_steps("SOAPdenovo")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.SOAPdenovo.start()
        self.step.update()

    def stepfinish(self):
        self.step.SOAPdenovo.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('config'):
            raise OptionError('必须输入组装配置config文件', code="32100201")
        if not self.option('kmerFreqCutoff'):
            raise OptionError('必须输入参数-d', code="32100202")
        if not self.option('sample_name'):
            raise OptionError('必须输入样品名称', code="32100203")
        if not self.option('kmer'):
            raise OptionError('必须输入kmer值', code="32100204")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 5
        self._memory =  "80G"


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BacDenovaAssAgent, self).end()


class BacDenovaAssTool(Tool):
    def __init__(self, config):
        super(BacDenovaAssTool, self).__init__(config)
        self.sample_name = self.option('sample_name')
        self._version = "SOAPdenovo2-src-r240"
        self.d= self.option('kmerFreqCutoff')
        self.kmer = self.option('kmer')
        self.config_file = self.option('config').prop['path']
        self.SOAPdenovo_path = '/bioinfo/metaGenomic/SOAPdenovo2/bin/'

    def run(self):
        """
        运行
        :return:
        """
        super(BacDenovaAssTool, self).run()
        have_result = self.run_SOAPdenovo2()  # 如果已有拼接结果，则stat为1，如果没有，则stat为0
        self.set_output(have_result)
        self.end()

    def run_SOAPdenovo2(self):
        self.logger.info(self.SOAPdenovo_path)
        path =self.work_dir + "/" + self.sample_name + '.' + str(self.kmer) + 'kmer_' + str(self.d) + '.scafSeq'
        cmd =self.SOAPdenovo_path +'SOAPdenovo2-127mer all -s %s -o %s -K %s -p 5 -d %s -D 1 -F T -u -R -M 1 -G -L' %(self.config_file,self.work_dir + "/" + self.sample_name + '.' + str(self.kmer) + 'kmer_' + str(self.d),self.kmer,self.d)
        self.logger.info(cmd)
        if os.path.exists(path):
            self.logger.info("%s.%s.%skmer.scafSeq已存在，跳过拼接" % (self.sample_name, self.kmer,self.d))
            result_stat = 1
        else:
            command = self.add_command("soapdenovo", cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行cmd完成")
            else:
                self.set_error("运行cmd运行出错!", code="32100201")
            result_stat = 0
        return result_stat

    def set_output(self, status):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if status == 0:
            self.logger.info("现在复制结果")
            if os.path.exists(self.output_dir + "/" + self.sample_name + '.' + str(self.kmer) + 'kmer_' + str(self.d) + '.scafSeq'):
                os.remove(self.output_dir + "/" + self.sample_name + '.' + str(self.kmer) + 'kmer_' + str(self.d) + '.scafSeq')
            os.link(self.work_dir + "/" + self.sample_name + '.' + str(self.kmer) + 'kmer_' + str(self.d) + '.scafSeq',self.output_dir + "/" + self.sample_name + '.' + str(self.kmer) + 'kmer_' + str(self.d) + '.scafSeq')
        self.logger.info("设置SOAPdenovo2分析结果目录成功")
