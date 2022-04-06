# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class BlastAgent(Agent):
    """
    ncbi blast+ 2.3.0
    """

    def __init__(self, parent):
        super(BlastAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "query_type", "type": "string"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "nr"},# 比对数据库 nt nr swissprot customer
            {"name": "blast", "type": "string"},  # 设定blast程序有blastn，blastx，blastp此处需要严格警告使用者必须选择正确的比对程序
            {"name": "ref", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "reference_type", "type": "string"},  # 参考序列(库)的类型  为nucl或者prot
            {"name": "top_num", "type": "int", "default": 1},  # top_num
            {"name": "align_len", "type": "int", "default": 50},  # align_len
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "identity", "type": "int", "default": 30},  # Identity
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
        ]
        self.add_option(options)
        self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="31100201")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')))
        if self.option("database") == 'Custom':
            if not self.option("ref").is_set:
                raise OptionError("使用自定义数据库模式时必须设置reference")
            if self.option('reference_type') not in ['nucl', 'prot']:
                raise OptionError('reference_type参考序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')))
        elif self.option("database") not in ["NT", "NR", 'Swiss-Prot']:
            raise OptionError("数据库%s不被支持", variables=(self.option("database")))
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s', variables=(self.option('evalue')))
        if not 0 < self.option('top_num') < 50:
            raise OptionError('序列比对保留数必须设置在1-50之间:%s', variables=(self.option('top_num')))
        if self.option('blast') not in ['blastn', 'blastp', 'blastx']:
            raise OptionError(
                '程序不试用于提供的查询序列和库的类型，请仔细检查，核酸比对核酸库只能使用blastn或者tblastn，\
                 核酸比对蛋白库只能使用blastp， 蛋白比对蛋白库只能使用blastp, 或者没有提供blast参数')
        return True

    def set_resource(self):
        self._cpu = 4
        self._memory = '40G'

    def end(self):
        super(BlastAgent, self).end()


class BlastTool(Tool):
    def __init__(self, config):
        super(BlastTool, self).__init__(config)
        if self.option("database") in ['NT']:
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/toolapps/NT_2018.07.26/blast/nt")
        elif self.option("database") in ['Swiss-Prot']:
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/toolapps/swissprot/blast/swissprot")
        elif self.option("database") == "NR":
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/toolapps/NR_2018.06.06/blast/nr")
        elif self.option("database") == 'Custom':
            self.db_path = os.path.join(self.work_dir, 'custom_blastdb')
        self.cmd_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"  # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        self.set_environ(BLASTDB=self.db_path)
        self.python = "/miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/toolapps/annno_taxon.py"

    def run_makedb_and_blast(self):
        """
        运行makeblastdb和blast
        :return:
        """
        db_name = os.path.splitext(os.path.basename(self.option("ref").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, "makeblastdb")
        cmd += " -dbtype %s -in %s -parse_seqids -out %s " % (self.option('reference_type'),
                                                              self.option("ref").prop['path'],
                                                              os.path.join(self.db_path, db_name))
        self.logger.info("开始运行makeblastdb，生成结果库文件放在工作目录的customer_blastdb下")
        makedb_obj = self.add_command("makeblastdb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("makeblastdb运行完成")
        else:
            self.set_error("makeblastdb运行出错!")

    def run_blast(self):
        """
        运行Blast
        :param db_name: blastdb名称
        :return:
        """
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, self.option('blast'))
        self.outputfile = os.path.join(self.output_dir, query_name + ".blast.m8.xls")
        if self.option("database") == 'Custom':
            db_name = os.path.splitext(os.path.basename(self.option("ref").prop['path']))[0]
            db_path = os.path.join(self.db_path, db_name)
        else:
            db_path = self.db_path
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt 6 -max_hsps 10 -num_threads 4 -max_target_seqs %s" % (
            self.option("query").prop['path'], db_path, self.outputfile,
            self.option("evalue"), self.option('top_num'))
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
        cmd = "{} {} -i {} -a {} -o {}".format(self.python, self.python_script, self.outputfile,
                                               self.option("database"), self.out)
        taxid_command = self.add_command("get_taxid", cmd)
        taxid_command.run()
        self.wait()
        if taxid_command.return_code == 0:
            self.logger.info("运行get_taxid完成")
        else:
            self.set_error("get_taxid运行出错!")

    def run(self):
        """
        运行
        :return:
        """
        super(BlastTool, self).run()
        if self.option("database") == 'Custom':
            self.run_makedb_and_blast()
            self.run_blast()
            self.end()
        elif self.option("database") == 'Swiss-Prot':
            self.run_blast()
            self.end()
        elif self.option("database") in ['NT',"NR"]:
            if int(self.option("top_num")) > 1:
                self.run_blast()
                self.end()
            elif int(self.option("top_num")) == 1:
                self.run_blast()
                self.get_taxid()
                self.end()