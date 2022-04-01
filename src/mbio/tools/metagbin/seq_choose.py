# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class SeqChooseAgent(Agent):
    """
    对所有样本抽取数据
    """
    def __init__(self, parent):
        super(SeqChooseAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "seed", "type": "int", "default": 100},#种子数
            {"name": "seq_num", "type": "int","default": "28800,14400,9600,7200,5760,4800"},#抽取数据的条数
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fastq'):
            raise OptionError('必须输入fastq文件', code="")
        if not self.option('seq_num'):
            raise OptionError('必须输入seq_num', code="")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '50G'

    def end(self):
        super(SeqChooseAgent, self).end()


class SeqChooseTool(Tool):
    def __init__(self, config):
        super(SeqChooseTool, self).__init__(config)
        self.seqtk = self.config.SOFTWARE_DIR + '/bioinfo/seq/seqtk-master/seqtk'
        self.seq_reads = self.work_dir + "/" + ''
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        wind_arrry = (self.option("seq_num")).split(',')
        #number = os.path.basename(self.option("fastq").prop['path']).split('.')[-2]
        #name = os.path.basename(self.option("fastq").prop['path']).split('.')[0]
        #for wind in wind_arrry:
            #self.out_fastq = self.work_dir + "/" +str(wind) + "/" +number + "/" + name +"_"+ ".fq"

    def run_seqtk(self):
        """
        用于筛选一定种子数一定条数的序列
        :return:
        """
        self.logger.info("开始抽取")
        self.choose_data = self.option("fastq").prop['path']
        for f in os.listdir(self.choose_data):
            file_path = os.path.join(self.choose_data, f)
            number = f.strip().split('.')[-2]
            wind_arrry = (self.option("seq_num")).split(',')
            #self.logger.info('文件名名称和路径: {}'.format(self.option("fastq")))
            #self.logger.info('文件名名称和路径{} {}'.format(number, path))
            name = f.strip().split('.')[0]
            for wind in wind_arrry:
                out_fastq = self.work_dir + "/" + str(wind) + "_" +number + "_"+ name +"_"+ ".fq"
                wind = int(wind)
                cmd = '{}seqtk.sh {} {} {} {}'.format(self.sh_path, self.seqtk, self.option("fastq").prop['path'],wind,out_fastq)
                self.logger.info(cmd)
                seqtk = 'seqtk' + str(wind)
                command = self.add_command(seqtk, cmd).run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("%s运行完成"%cmd)
                else:
                    self.set_error("%s运行失败", variables=(cmd), code="")

    def set_output(self):
        """
        设置结果文件
        :return:
        """
        wind_arrry = (self.option("seq_num")).split(',')
        path = self.option("fastq").prop['path']
        number = os.path.basename(path).split('.')[2]
        name = os.path.basename(self.option("fastq").prop['path']).split('.')[0]
        for wind in wind_arrry:
            os.makedirs(self.output_dir + "/"+str(wind),0755)
            if os.path.exists(self.output_dir + "/" +str(wind) + "/"+number + "_"+ name +"_"+ ".fq"):
                os.remove(self.output_dir + "/" +str(wind) + "/"+number + "_"+ name +"_"+ ".fq")
            self.logger.info('正在生成结果')#
            os.link(self.work_dir +  "/" + str(wind) + "_" +number + "_"+ name +"_"+ ".fq", self.output_dir + "/" +str(wind) + "/"+number + "_"+ name +"_"+ ".fq")
        self.end()

    def run(self):
        super(SeqChooseTool, self).run()
        self.run_seqtk()
        self.set_output()
        self.end()