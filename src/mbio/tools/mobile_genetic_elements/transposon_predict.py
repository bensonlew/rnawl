# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.06.18

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from Bio import SeqIO


class TransposonPredictAgent(Agent):
    """
    根据输入文件为组装序列，进行Transposon的预测，这里预测复杂转座子
    """

    def __init__(self, parent):
        super(TransposonPredictAgent, self).__init__(parent)
        options = [
            {"name": "input_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gbk", "type": "bool", "default": False},  #true是生成gbk文件，否则不生成
            {"name": "identity", "type": "int", "default": 80},   #blastn与数据库db比对的identity，这里是筛选的阀值
            {"name": "coverage", "type": "int", "default": 80},   #blastn与数据库db比对的coverage，这里是筛选的阀值
            {"name": "distance", "type": "int", "default": 20000},  #Transposon的最大长度是20000bp
            {"name": "extend", "type": "int", "default": 20},  #向外延伸长度
            {"name": "gene_gff", "type": "string"},  # 预测的基因gff文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_fa").is_set:
            raise OptionError("必须设置参数input_fa文件!")

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(TransposonPredictAgent, self).end()


class TransposonPredictTool(Tool):
    def __init__(self, config):
        super(TransposonPredictTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "/program/Python35/bin:"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib:"
        self.set_environ(PATH=self.path,LD_LIBRARY_PATH=self.lib)
        self.python_s = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/tncomp_finder-master/TnComp_finder.py"
        self.fasta = self.option("input_fa").prop['path']
        self.python = "/program/Python35/bin/python3"
        self.python_path = '/miniconda2/bin/python'
        self.script_path = self.config.PACKAGE_DIR + '/mobile_genetic_elements/'
        self.names = os.path.basename(self.option("input_fa").prop['path']).split(".fna")[0]

    def run_transposon(self):
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        os.mkdir(self.work_dir + "/temp")
        if os.path.exists(self.work_dir + "/result"):
            shutil.rmtree(self.work_dir + "/result")
        os.mkdir(self.work_dir + "/result")
        self.split_fasta(self.fasta, self.work_dir + "/temp")
        n =1
        for i in os.listdir(self.work_dir + "/temp"):
            sample = i.split(".fna")[0]
            file = self.work_dir + "/temp/" + i
            cmd = "{} {} -f {} -o {} -p 4 -i {} -c {} -d {} ".format(self.python, self.python_s, file,
                                                                     self.work_dir + "/tmp"+str(n), self.option("identity"),
                                                                     self.option("coverage"), self.option("distance"))
            if self.option("gbk"):
                cmd += "-g "
            if self.option("extend"):
                cmd += "-e {}".format(self.option("extend"))
            command = self.add_command("run_transposon"+str(n), cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("get_transposon_result处理成功！")
            else:
                self.set_error("get_transposon_result文件处理失败!")
            if os.path.exists(self.work_dir + "/tmp" + str(n) + "/" + sample + "_composite.txt"):
                try:
                    cmd1 = '{} {}transposon_chuli.py --g {} --o {}'.format(self.config.SOFTWARE_DIR + self.python_path, self.script_path,
                                                                           self.work_dir + "/tmp" + str(
                                                                               n) + "/" + sample + "_composite.txt",
                                                                           self.work_dir + "/result/" + sample + ".transposon.xls")
                    self.logger.info(cmd1)
                    os.system(cmd1)
                    self.logger.info("transposon_chuli成功")
                except:
                    self.set_error("transposon_chuli失败")
            n +=1

    def run_transposon_stat(self):
        cmd2 = '{} {}transposon_stat.py --d {} --o {}'.format(self.python_path, self.script_path, self.work_dir + "/result",
                                                               self.work_dir + "/" + self.names + ".transposon.result.xls")
        command2 = self.add_command('get_transposon_result', cmd2).run()
        self.wait(command2)
        self.logger.info(command2.return_code)
        if command2.return_code == 0:
            self.logger.info("get_transposon_result处理成功")
        else:
            self.set_error("get_transposon_result文件处理失败")

    def chuli_result(self):
        """
        根据整理的结果与gff编码取比较，去除部分错误Itransposon
        :return:
        """
        cmd = "{} {}gff_mobile.py --g {} --m {} --o {} --s 2 --e 3 --i 1".format(self.python_path, self.script_path, self.option("gene_gff"), self.work_dir + "/" + self.names + ".transposon.result.xls", self.output_dir + '/' + self.names + ".transposon.xls")
        command = self.add_command("chuli_result", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("chuli_result运行完成！")
        else:
            self.set_error("chuli_result运行完成运行出错!")

    def run(self):
        super(TransposonPredictTool, self).run()
        self.run_transposon()
        if len(os.listdir(self.work_dir + "/result")) > 0:
            self.run_transposon_stat()
            self.chuli_result()
        self.end()

    def split_fasta(self, fasta, dir):
        for seq_record in SeqIO.parse(fasta, "fasta"):
            id = seq_record.id
            SeqIO.write(seq_record,dir + "/" +id + ".fna", "fasta")