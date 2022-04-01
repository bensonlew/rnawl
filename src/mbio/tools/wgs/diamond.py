# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180502

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class DiamondAgent(Agent):
    """
    参考基因组用diamond的数据库注释
    比对数据库 nr kegg go uniprot eggnog
    eg：diamond  blastx --db  /mnt/ilustre/users/long.huang/DataBase/NR/2017-8-30/nr --query
    /mnt/ilustre/users/qingmei.cui/newmdt/Genome/Arabidopsis_thaliana/NCBI/TAIR10_2011-05-11_Columbia/GCF_TAIR10/02.split/sub.0.fa
     --evalue 10e-10 --outfmt 5 --threads 8 --out
    /mnt/ilustre/users/qingmei.cui/newmdt/Genome/Arabidopsis_thaliana/NCBI/TAIR10_2011-05-11_Columbia/GCF_TAIR10/03.NR/sub.0.fa.nr.blast
    lasted modified by hongdong@20181011
    """
    def __init__(self, parent):
        super(DiamondAgent, self).__init__(parent)
        options = [
            {"name": "blast", "type": "string", "default": "blastx"},  # 设定diamond程序有blastp，blastx
            {"name": "database", "type": "string"},                     # 输入数据库，目前输入数据库的的dmnd文件；
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入拆开后的fasta,
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，数字默认为5，输出xml
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "db_type", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('diamond')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.diamond.start()
        self.step.update()

    def step_end(self):
        self.step.diamond.finish()
        self.step.update()

    def check_options(self):
        if self.option('blast') not in ['blastp', 'blastx']:
            raise OptionError('比对类型只能是blastp或blastx:%s',variables=(self.option('blast')), code="34502001")
        if not self.option("database"):
            raise OptionError("请设置database", code="34502002") # 必须有数据库
        if not self.option("query"):
            raise OptionError("请设置query", code="34502003") # 必须有序列
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s',variables=(self.option('evalue')), code="34502004")
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = self.option('num_threads')
#        self._memory = self.option('num_mem')+"G"
        self._memory = '30G'

    def end(self):
        super(DiamondAgent,self).end()


class DiamondTool(Tool):
    def __init__(self, config):
        super(DiamondTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/align/diamond-0.9.11/")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/bioinfo/align/diamond-0.9.11/")
        self._version = "v0.9.10.11"
        self.cmd_path = "bioinfo/align/diamond-0.9.11/diamond"

    def Diamond(self):
        """
        要重新写下！！！
        :return:
        """
        if self.option("db_type"):
            name = os.path.basename(self.option("query").prop["path"]) + ".{}.blast".format(self.option("db_type"))
        else:
            name = os.path.basename(self.option("query").prop["path"]) + ".blast"  # 需要加进去数据库的type吗?如nr.blast
        cmd = self.cmd_path
        cmd += " {} --db {} --query {} --evalue {} --outfmt {} --threads {} --out {}" \
            .format(self.option("blast"), self.option("database"), self.option("query").prop['path'],
                    self.option("evalue"), self.option("outfmt"), self.option("num_threads"),
                    self.output_dir + "/"+name)
        self.logger.info(cmd)
        self.logger.info("开始进行Diamond")
        command = self.add_command("diamond", cmd).run()  # nranno必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Diamond完成！")
        else:
            self.set_error("Diamond出错！", code="34502001")
            self.set_error("Diamond出错！", code="34502004")

    def run(self):
        super(DiamondTool, self).run()
        self.Diamond()
        self.end()
