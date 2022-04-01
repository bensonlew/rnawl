# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class DiamondAgent(Agent):
    """
    ncbi blast+ 2.3.0
    """
    def __init__(self, parent):
        super(DiamondAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "query_type", "type": "string"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "NR"},# 比对数据库 nr swissprot customer
            {"name": "blast", "type": "string"},  # 设定blast程序有blastx，blastp此处需要严格警告使用者必须选择正确的比对程序
            {"name": "ref", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "reference_type", "type": "string"},  # 参考序列(库)的类型  prot
            {"name": "top_num", "type": "int", "default": 1},  # top_num
            {"name": "align_len", "type": "int", "default": 50},  # align_len
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "identity", "type": "int", "default": 30},  # Identity
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
        ]
        self.add_option(options)
        self.queue = 'BLAST'  # 投递到指定的队列BLAST
        self._memory_increase_step = 30

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="31100201")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')))
        if self.option("database") == 'Custom':
            if not self.option("ref").is_set:
                raise OptionError("使用自定义数据库模式时必须设置reference")
            if self.option('reference_type') not in ['prot']:
                raise OptionError('reference_type参考序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')))
        elif self.option("database") not in ["NR", 'Swiss-Prot']:
            raise OptionError("数据库%s不被支持", variables=(self.option("database")))
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s', variables=(self.option('evalue')))
        if not 0 < self.option('top_num') < 50:
            raise OptionError('序列比对保留数必须设置在1-50之间:%s', variables=(self.option('top_num')))
        if self.option('blast') not in ['blastp', 'blastx']:
            raise OptionError(
                '程序不试用于提供的查询序列和库的类型，请仔细检查，核酸比对核酸库只能使用，\
                 核酸比对蛋白库只能使用blastx， 蛋白比对蛋白库只能使用blastp, 或者没有提供blast参数')
        return True

    def set_resource(self):
        self._cpu = 4
        self._memory = '40G'

    def end(self):
        super(DiamondAgent, self).end()


class DiamondTool(Tool):
    def __init__(self, config):
        super(DiamondTool, self).__init__(config)
        if self.option("database") in ['Swiss-Prot']:
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/toolapps/swissprot/diamond/swissprot")
        elif self.option("database") == "NR":
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/toolapps/NR_2018.06.06/diamond/nr")
        elif self.option("database") == 'Custom':
            self.db_path = os.path.join(self.work_dir, 'custom_blastdb')
            if not os.path.exists(self.db_path):
                os.mkdir(self.db_path)
        self.cmd_path = "bioinfo/align/diamond-0.8.35/diamond"  # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        self.set_environ(BLASTDB=self.db_path)
        self.python = "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/toolapps/annno_taxon.py"

    def run_makedb_and_blast(self):
        """
        运行diamond makedb
        :return:
        """
        db_name = os.path.splitext(os.path.basename(self.option("ref").prop['path']))[0]
        cmd = "%s makedb --in %s --db %s " % (self.cmd_path, self.option("ref").prop['path'], os.path.join(self.db_path, db_name))
        self.logger.info("开始运行makeblastdb，生成结果库文件放在工作目录的customer_blastdb下")
        makedb_obj = self.add_command("diamond_makedb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("diamond_makedb运行完成")
        else:
            self.set_error("diamond_makedb运行出错!")

    def run_blast(self):
        """
        运行diamond
        :param db_name: blastdb名称
        :return:
        """
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        self.outputfile = os.path.join(self.output_dir, query_name + ".diamond.m8.xls")
        if self.option("database") == 'Custom':
            db_name = os.path.splitext(os.path.basename(self.option("ref").prop['path']))[0]
            db_path = os.path.join(self.db_path, db_name)
        else:
            db_path = self.db_path
        cmd ="{} {} --db {} --query {} --top {} --id {} -o {} -e {} --threads 4".format(self.cmd_path, self.option("blast"), db_path, self.option("query").prop['path'], self.option("top_num"), self.option("identity"), self.outputfile, self.option("evalue"))
        self.logger.info("开始运行blast")
        blast_command = self.add_command("blast", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.logger.info("运行blast完成")
        else:
            self.set_error("blast运行出错!")

    def get_taxid(self):
        """

        """
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        self.out = os.path.join(self.output_dir, query_name + ".taxon.xls")
        cmd = "{} {} -i {} -a {} -o {}".format(self.python, self.python_script, self.outputfile, self.option("database"), self.out)
        taxid_command = self.add_command("get_taxid", cmd, ignore_error=True)
        taxid_command.run()
        self.wait()
        if taxid_command.return_code == 0:
            self.logger.info("运行get_taxid完成")
        elif taxid_command.return_code in [1, -6, -9]:
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("get_taxid运行出错!")


    def run(self):
        """
        运行
        :return:
        """
        super(DiamondTool, self).run()
        if self.option("database") == 'Custom':
            self.run_makedb_and_blast()
            self.run_blast()
            self.end()
        elif self.option("database") == 'Swiss-Prot':
            self.run_blast()
            self.end()
        elif self.option("database") in ["NR"]:
            if int(self.option("top_num")) >1 :
                self.run_blast()
                self.end()
            else:
                self.run_blast()
                self.get_taxid()
                self.end()