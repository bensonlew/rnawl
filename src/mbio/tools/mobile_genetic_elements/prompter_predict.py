# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
from Bio import SeqIO

class PrompterPredictAgent(Agent):
    """
    启动子预测 v1.0
    author: gaohao
    last_modify: 20200911
    """
    def __init__(self, parent):
        super(PrompterPredictAgent, self).__init__(parent)
        options = [
            {"name": "assemble", "type": "infile", "format": "sequence.fasta"},  # 拼接序列
            {"name": "window_size", "type": "int", "default": 100},  # 窗口大小
            {"name": "cut", "type": "int", "default": 200},  # orf间间隔长度
            {"name": "seq_type", "type": "int", "default": 0},  # 0是多序列，1是小于10M的基因组，2是大于10M的基因组
            {"name": "sample", "type": "string"},  # 样品名
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "gff", "type": "string"},  # gff文件路径
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gff"):
            raise OptionError("必须设置输入gff文件", code="32201501")
        if self.option('seq_type') not in [0, 1, 2]:
            raise OptionError("错误的序列类型,只能为0,1,2", code="32201502")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(PrompterPredictAgent, self).end()


class PrompterPredictTool(Tool):
    def __init__(self, config):
        super(PrompterPredictTool, self).__init__(config)
        self._version = "1.0"
        self.sh_path = self.config.SOFTWARE_DIR + '/bioinfo/gene-structure/prompredict/prom.sh'
        self.genome1 = self.config.SOFTWARE_DIR + '/bioinfo/gene-structure/prompredict/PromPredict_genome_V1'
        self.genome2 = self.config.SOFTWARE_DIR + '/bioinfo/gene-structure/prompredict/PromPredict_genome_V2'
        self.mulseq = self.config.SOFTWARE_DIR + '/bioinfo/gene-structure/prompredict/PromPredict_mulseq'
        self.python_path = '/miniconda2/bin/python'
        self.pre_path = self.config.PACKAGE_DIR + '/mobile_genetic_elements/pro_prepare.pl'
        self.script_path = self.config.PACKAGE_DIR + '/mobile_genetic_elements/'
        self.bedtools = self.config.SOFTWARE_DIR + "/bioinfo/seq/bedtools-2.25.0/bin/bedtools"
        self.sh = "../../../../../.." + self.config.PACKAGE_DIR  + '/mobile_genetic_elements/bedtools.sh'

    def run(self):
        """
        运行
        :return:
        """
        super(PrompterPredictTool, self).run()
        self.predict_prom()
        self.run_gff_change()
        self.run_bedtools()
        self.get_result()
        self.end()

    def predict_prom(self):
        if self.option('seq_type') == 0:
            software = self.mulseq
        elif self.option('seq_type') == 1:
            software = self.genome1
        elif self.option('seq_type') == 2:
            software = self.genome1
        self.change_fasta(self.option('assemble').prop['path'], self.work_dir + "/pre.fasta")
        cmd1 = '{} {} {} {} {}'.format(self.sh_path, software, self.work_dir + "/pre.fasta", self.option('window_size'), 'default')
        try:
            os.system(cmd1)
            self.logger.info("启动子预测成功")
        except:
            self.set_error("启动子预测失败")
        pro_result = self.work_dir + '/pre_PPde.txt'
        prefix_path = self.work_dir + '/promoter'
        cmd2 = '{} {}prome_chuli.py --g {} --o {}'.format(self.python_path, self.script_path, pro_result, prefix_path)
        command2 = self.add_command('stat', cmd2).run()
        self.wait(command2)
        self.logger.info(command2.return_code)
        if command2.return_code == 0:
            self.logger.info("pre_PPde.txt文件处理成功")
        else:
            self.set_error("pre_PPde.txt文件处理失败")

    def run_gff_change(self):
        cmd2 = '{} {}gff_chuli.py --g {} --o {}'.format(self.python_path, self.script_path, self.option("gff"), self.work_dir + "/gff.bed")
        command2 = self.add_command('gff_change', cmd2).run()
        self.wait(command2)
        self.logger.info(command2.return_code)
        if command2.return_code == 0:
            self.logger.info("gff_change处理成功")
        else:
            self.set_error("gff_change文件处理失败")

    def run_bedtools(self):
        """
        判断预测的启动子和编码区基因不在一个区间，
        :return:
        """
        cmd2 = '{} {} {} {} {}'.format(self.sh, self.bedtools, self.work_dir + '/promoter.bed', self.work_dir + "/gff.bed", self.work_dir + "/result.xls")
        command2 = self.add_command('run_bedtools', cmd2).run()
        self.wait(command2)
        self.logger.info(command2.return_code)
        if command2.return_code == 0:
            self.logger.info("run_bedtools处理成功")
        else:
            self.set_error("run_bedtools文件处理失败")

    def get_result(self):
        cmd2 = '{} {}get_prompter_result.py --p {} --s {} --o {}'.format(self.python_path, self.script_path, self.work_dir + '/promoter.xls',
                                       self.work_dir + "/result.xls", self.output_dir + "/" + self.option("sample")+ ".prompter.xls")
        command2 = self.add_command('get_result', cmd2).run()
        self.wait(command2)
        self.logger.info(command2.return_code)
        if command2.return_code == 0:
            self.logger.info("get_result处理成功")
        else:
            self.set_error("get_result文件处理失败")

    def change_fasta(self, fa, out):
        seq_list = []
        for seq_record in SeqIO.parse(fa, "fasta"):
            seq_record.description = ""
            seq_list.append(seq_record)
        SeqIO.write(seq_list, out, "fasta")