# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class BlastAgent(Agent):
    """
    ncbi blast+   请详细编写使用说明
    version 1.0  你的程序版本
    author: guoquan  作者
    last_modify: 2015.9.21  最后修改日期
    """

    def __init__(self, parent):
        super(BlastAgent, self).__init__(parent)
        options = [
            {"name": "customer_mode", "type": "bool", "default": False},  # customer 自定义数据库
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "database", "type": "string", "default": "nr"},  # 比对数据库 nt nr string GO swissprot uniprot KEGG
            {"name": "reference", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "output", "type": "outfile", "format": "sequence.fasta"}  # cpu数
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query")
        if self.option("customer_mode") is True and not self.option("reference").is_set:
            raise OptionError("使用自定义数据库模式时必须设置reference")
        if self.option("database") not in ["nt", "nr", "string"]:
            raise OptionError("数据库%s不被支持" % self.option("database"))
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        可以在其中编写复杂逻辑，比如通过判断输入文件大小来动态给出所需资源，使其尽量匹配实际情况
        :return:
        """
        self._cpu = 10
        self._memory = ''

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [r".*\.xml", "align.blast.blastxml", "Blast比对结果，XML格式"]
        ])
        super(BlastAgent, self).end()


class BlastTool(Tool):
    def __init__(self, config):  # 注意 初始化Tool子类时是需要带参数的
        super(BlastTool, self).__init__(config)  # 调用父类初始化
        self._version = "2.2.31"  # 定义程序版本
        self.db_path = os.path.join(self.config.SOFTWARE_DIR, "align/ncbi/db/")
        self.cmd_path = "align/ncbi/blast-2.2.31+/bin"   # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        self.relation = {
            "blastn": ["DNA", "DNA"],
            "blastp": ["Protein", "Protein"],
            "blastx": ["DNA", "Protein"],
            "tblastn": ["Protein", "DNA"]
        }
        self.db_type = {
            "nt": "DNA",
            "nr": "Protein",
            "strings": "Protein",
            "go": "Protein",
            "swissprot": "Protein",
            "uniprot": "Protein",
            "kegg": "Protein"
        }

    def get_blast_type(self):
        """
        根据输入文件及数据库类型获取blast类型名

        :return: string blastn/ blastp/ blastx
        """
        input_type = self.option("query").prop['seq_type']
        if self.option("customer_mode"):
            blast_db_type = self.option("reference").prop['seq_type']
        else:
            blast_db_type = self.db_type[self.option("database").lower()]
        for key, value in self.relation.items():
            if input_type == value[0] and blast_db_type == value[1]:
                return key
        raise Exception("不支持此类型的序列比对: input:%s  reference: %s" % (input_type, blast_db_type))

    def run_makedb_and_blast(self):
        """
        运行makeblastdb和blast

        :return:
        """
        db_name = os.path.basename(self.option("reference").prop['path'])
        cmd = os.path.join(self.cmd_path, "makeblastdb")
        seq_type = "nucl" if self.option("reference").prop['seq_type'] == "DNA" else "prot"
        cmd += " -dbtype %s -in %s -parse_seqids -title %s -out %s " % (seq_type, self.option("reference").prop['path'], db_name, db_name)
        self.logger.info("开始运行makeblastdb")
        makedb_obj = self.add_command("makeblastdb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("makeblastdb运行完成")
            self.run_blast(db_name)
        else:
            self.set_error("makeblastdb运行出错!")

    def run_blast(self, db_name):
        """
        运行Blast

        :param db_name: blastdb名称
        :return:
        """

        cmd = os.path.join(self.cmd_path, self.get_blast_type())
        outputfile = os.path.join(self.output_dir, os.path.basename(self.option("query").prop['path']) + "_vs_"
                                  + db_name + ".xml")
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt 5 -max_hsps 10 -max_target_seqs 10 -num_threads %s" % (
            self.option("query").prop['path'], db_name, outputfile, self.option("evalue"), self.option("num_threads"))
        self.logger.info("开始运行blast")  # 尽量多给出提示，便于调试和阅读程序进度
        blast_command = self.add_command("blast", cmd)   # 添加命令对象
        if self.option("customer_mode"):
            self.db_path = os.getcwd()
        self.set_environ(BLASTDB=self.db_path)   # 设置运行命令所需的环境变量，对整个Tool运行时生效,每个Tool应该设置环境变量保证自身的运行
        blast_command.run()   # 开始运行命令
        self.wait()  # 等待命令结束
        if blast_command.return_code == 0:  # 判断命令是否正常完成，需要根据命令实际情况编写 也可编写_check函数
            self.logger.info("运行blast完成")
            self.end()        # 设置Tool为正常完成状态，并将状态发送远程Agent
        else:
            self.set_error("blast运行出错!")  # 设置Tool为异常错误状态，并将状态发送远程Agent
            # 也可获取错误类型，根据情况调整参数后重新运行命令，使命令正确完成
            # 或者发送自定义State状态给远程Agent 由Module或Workflow中定义运行逻辑

    def run(self):
        """
        运行
        :return:
        """
        super(BlastTool, self).run()
        if self.option("customer_mode"):
            self.run_makedb_and_blast()
        else:
            db_name = self.option("database")
            self.run_blast(db_name)
