# -*- coding: utf-8 -*-
# __author__ = 'Shijin'
#last_modify by: shaohua.yuan

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from biocluster.core.exceptions import OptionError


class MetaDiamondAgent(Agent):
    """
    diamond version: 0.8.35
    version 1.0
    author: shijin
    last_modify by: shaohua.yuan
    last_modify : 20170829
    """
    def __init__(self, parent):
        super(MetaDiamondAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "query_type", "type": "string", "default": "prot"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "plant"},
            # 比对数据库 plant, nr, etc.
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，数字默认为5，输出xml
            {"name": "blast", "type": "string", "default": "blastp"},  # 设定diamond程序有blastp，blastx
            {"name": "reference", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 9},  # cpu数
            {"name": "sensitive", "type": "int", "default": 0},
            {"name": "unalign", "type": "int", "default": 0},         ## report unaligned queries (0=no, 1=yes)
            {"name": "target_num", "type": "int", "default": 1}         ### maximum number of target sequences to report alignments for
            #{"name":"unalign_file","type":"outfile,"format":"sequence.fasta""}         ####未比对上的序列fasta
            ]
        self.add_option(options)
        self.step.add_steps('diamond')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def step_start(self):
        self.step.diamond.start()
        self.step.update()

    def step_end(self):
        self.step.diamond.finish()
        self.step.update()

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="31101801")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')), code="31101802")
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s', variables=(self.option('evalue')), code="31101803")
        if not 0 <= self.option("sensitive") <= 2:
            raise OptionError('敏感度设定必须为[0-2]之间：%s', variables=(self.option('evalue')), code="31101804")
        return True

    def set_resource(self):
        self._cpu = self.option('num_threads')
        self._memory = '30G'  # modified by guhaidong ,20G -> 50G 20171024

    def end(self):
        super(MetaDiamondAgent, self).end()


class MetaDiamondTool(Tool):
    def __init__(self, config):
        super(MetaDiamondTool, self).__init__(config)
        self._version = "0.8.35"
        self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/align/diamond")
        self.cmd_path = "bioinfo/align/diamond-0.8.35"   # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        if self.option("query_type") == "nucl":
            self.blast_type = "blastx"
        else:
            self.blast_type = "blastp"
        self.ori = []
        self.repl = []

    def run_makedb_and_diamond(self):
        """
        创建diamond数据库并运行diamond

        :return:
        """
        db_name = os.path.splitext(os.path.basename(self.option("reference").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, "makedb")
        self.db_path = os.path.join(self.work_dir, 'diamond')
        cmd += " makedb -in {} -d {}".format(self.option("reference").prop['path'], db_name)
        self.logger.info("开始创建diamond数据库，生成结果库文件放在工作目录的customer_blastdb下")
        makedb_obj = self.add_command("makedb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("创建diamond数据库完成")
            self.run_diamond(db_name)
        else:
            self.set_error("创建diamond数据库出错!", code="31101801")

    def run_diamond(self, db_name):
        """
        运行diaomond

        :param db_name: blastdb名称
        :return:
        """
        if self.option("database") in ['kegg_v94.2']:
            db = os.path.join(self.db_path, "kegg_v94.2")
        else:
            db = os.path.join(self.db_path, db_name)
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, "diamond")
        outputfile1 = os.path.join(self.output_dir, query_name + "_vs_" + db_name)
        outfmt = self.option('outfmt')
        if self.option('outfmt') == 5:
            outputfile = outputfile1 + '.xml'  # outfmt默认为5
        elif self.option('outfmt') == 6:       #zouguanqing 20180526
            outputfile = outputfile1 + '.xls'
        else:
            outputfile = outputfile1 + '.xml'
        if os.path.exists(outputfile):
            os.remove(outputfile)
        #outfmt = 5
        cmd += " {} -q {} -d {} -o {} -e {} -f {} -p {} -k {}".format(
            self.blast_type, self.option("query").prop['path'], db, outputfile,
            self.option("evalue"), outfmt, self.option("num_threads"),self.option("target_num"))
        if self.option("sensitive") == 1:
            cmd += " --sensitive"
        elif self.option("sensitive") == 2:
            cmd += " --more-sensitive"
        if self.option("unalign") == 1:
            unalignfile = os.path.join(self.output_dir,query_name + "_unalign.fasta")
            cmd += " --unal 1 --un {}".format(unalignfile)
        self.logger.info("开始运行blast")
        blast_command = self.add_command("diamond", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.logger.info("运行diamond完成")
            self.logger.info(outputfile)
            if self.option("outfmt") == 5:         #zouguanqing 20180526 增加判断
                self.change_version(outputfile)
            else :
                self.end()   #guanqing.zou 20180528
        elif blast_command.return_code == None:
            self.logger.info("重新运行diamond")
            blast_command.rerun()
            self.wait(blast_command)
            if blast_command.return_code == 0:
                self.logger.info("重新运行diamond成功")  #zouguanqing 20180526 增加判断
                # self.end()
                if self.option("outfmt") == 5:
                    self.change_version(outputfile)
                else:
                    self.end()   #guanqing.zou 20180528

        else:
            self.set_error("diamond运行出错!", code="31101802")
            # raise Exception("diamond运行出错!")

    def run(self):
        """
        运行
        :return:
        """
        super(MetaDiamondTool, self).run()
        if self.option("database") == 'customer_mode':
            self.run_makedb_and_diamond()
        else:
            self.run_diamond(self.option("database"))

    def change_version(self, outputfile):
        path = outputfile
        with open(path,"r") as file:
            for line in file:
                line = line.strip()
                if line.lstrip().startswith("<Hit_id>"):
                    m = re.match("<Hit_id>(.+)</Hit_id>", line.lstrip())
                    if m:
                        self.ori.append(m.group(1))
                        line = file.next()
                        n = re.match("<Hit_def>(.+)</Hit_def>", line.lstrip())
                        try:
                            self.repl.append(n.group(1))
                        except:
                            print line
        with open(path,"r") as file, open(path + "_new", "w") as fw:
            i = 0
            for line in file:
                if line.lstrip().startswith("<BlastOutput_db>"):
                    line = line.replace("<BlastOutput_db>", "<BlastOutput_db>" + self.option("database"))
                if line.lstrip().startswith("<BlastOutput_version>"):
                    line = line.replace("diamond 0.8.35", "BLASTX 2.3.0+")
                if line.lstrip().startswith("<Hit_id>"):
                    m = re.match("<Hit_id>(.+)</Hit_id>", line.lstrip())
                    if m:
                        line = line.replace(self.ori[i],self.repl[i])
                if line.lstrip().startswith("<Hit_def>"):
                    m = re.match("<Hit_def>(.+)</Hit_def>", line.lstrip())
                    if m:
                        line = line.replace(self.repl[i],self.ori[i])
                        i += 1
                fw.write(line)
        # os.system("mv {} {}".format(path + "_new", path))
        os.remove(path)
        os.link(path + "_new", path)
        self.end()
