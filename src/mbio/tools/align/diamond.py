# -*- coding: utf-8 -*-
# __author__ = 'Shijin'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
import shutil
import xml.etree.ElementTree as ET
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
import time

class DiamondAgent(Agent):
    """
    diamond version: 0.8.35
    version 1.0
    author: shijin
    last_modify: 20170316
    """
    def __init__(self, parent):
        super(DiamondAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "query_type", "type": "string", "default": "nucl"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "plant"},
            # 比对数据库 plant, nr, etc.
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，数字默认为5，输出xml
            {"name": "blast", "type": "string", "default": "blastp"},  # 设定diamond程序有blastp，blastx
            {"name": "reference", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 9},  # cpu数
            {"name": "sensitive", "type": "int", "default": 2},
            {"name": "version", "type": "string", "default": "1"} # diamond 版本 可选2
            ]
        self.add_option(options)
        self.step.add_steps('diamond')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.queue = 'BLAST'  # 投递到指定的队列BLAST
        self._memory_increase_step = 20  # 每次重运行增加10G内存 add by qingchen.zhang @ 20191204

    def step_start(self):
        self.step.diamond.start()
        self.step.update()

    def step_end(self):
        self.step.diamond.finish()
        self.step.update()

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="31101201")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')), code="31101202")
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：{}',variables=(self.option('evalue')), code="31101203")
        if not 0 <= self.option("sensitive") <= 2:
            raise OptionError('敏感度设定必须为[0-2]之间：%s', variables=(self.option('evalue')), code="31101204")
        return True

    def set_resource(self):
        self._cpu = self.option('num_threads')
        self._memory = '30G'

    def end(self):
        super(DiamondAgent, self).end()


class DiamondTool(Tool):
    def __init__(self, config):
        super(DiamondTool, self).__init__(config)
        if self.option("version") == "1":
            self._version = "0.8.35"
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/align/diamond")
            self.cmd_path = "bioinfo/align/diamond-0.8.35"  # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        else:
            self._version = "2.0.13"
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/align/diamond")
            self.cmd_path = "bioinfo/align/diamond-2.0.13"   # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        if self.option("query_type") == "nucl":
            self.blast_type = "blastx"
        else:
            self.blast_type = "blastp"
        # self.mongodb_nr = Config().biodb_mongo_client.sanger_biodb.NR_sequence
        if self.option("database") in ['nr']:
            self.mongodb_nr = Config().get_mongo_client(mtype="ref_rna", ref=True)[Config().get_mongo_dbname("ref_rna", ref=True)].NR_sequence
        elif self.option("database") in ['nr_v20200604']:
            self.mongodb_nr = Config().get_mongo_client(mtype="ref_rna", ref=True)[Config().get_mongo_dbname("ref_rna", ref=True)]['NR_sequence_20200604']
        elif self.option("database") in ['nr_v20210917']:
            self.mongodb_nr = Config().get_mongo_client(mtype="ref_rna", ref=True)[Config().get_mongo_dbname("ref_rna", ref=True)]['NR_sequence_202110']
        self.ori = []
        self.repl = []
        self.process_rerun = 0

    def run_makedb_and_diamond(self):
        """
        创建diamond数据库并运行diamond

        :return:
        """
        db_name = os.path.splitext(os.path.basename(self.option("reference").prop['path']))[0]
        #cmd = os.path.join(self.cmd_path, "makedb")
        cmd = os.path.join(self.cmd_path, "diamond")  #guanqing.zou 20180525
        self.db_path = os.path.join(self.work_dir, 'diamond')
        #cmd += " makedb -in {} -d {}".format(self.option("reference").prop['path'], db_name)
        cmd += " makedb --in {} -d {}".format(self.option("reference").prop['path'], db_name)  #guanqing 20180525
        self.logger.info("开始创建diamond数据库，生成结果库文件放在工作目录的customer_blastdb下")
        makedb_obj = self.add_command("makedb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("创建diamond数据库完成")
            self.run_diamond(db_name)
        else:
            self.set_error("创建diamond数据库出错!", code="31101201")

    def run_diamond(self, db_name):
        """
        运行diaomond

        :param db_name: blastdb名称
        :return:
        """
        if self.option("database") in ['kegg']:
            db = os.path.join(self.db_path, "kegg_v94.2")
        else:
            db = os.path.join(self.db_path, db_name)
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, "diamond")
        outputfile = os.path.join(self.output_dir, query_name + "_vs_" + db_name)
        outfmt = self.option('outfmt')
        # if self.option('outfmt') == 5:
        outputfile += '.xml'  # outfmt默认为5
        #outfmt = 5  #guanqing 20180525
        cmd += " {} -q {} -d {} -o {} -e {} -f {} -p {} -k 5".format(
            self.blast_type, self.option("query").prop['path'], db, outputfile,
            self.option("evalue"), outfmt, self.option("num_threads"))
        if self.option("sensitive") == 1:
            cmd += " --sensitive"
        elif self.option("sensitive") == 2:
            cmd += " --more-sensitive"
        self.logger.info("开始运行blast")
        blast_command = self.add_command("diamond", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.logger.info("运行diamond完成")
            self.logger.info(outputfile)
            # if db_name in ["nr", "animal", "fungi", "metazoa", "plant", "protist", "vertebrates"]:
            #     self.get_nrxml_gi_description(outputfile)
            self.change_version(outputfile)
        elif blast_command.return_code == None:
            self.logger.info("重新运行diamond")
            blast_command.rerun()
            self.wait(blast_command)
            if blast_command.return_code == 0:
                self.logger.info("重新运行diamond成功")
                # if db_name in ["nr", "animal", "fungi", "metazoa", "plant", "protist", "vertebrates"]:
                    # self.get_nrxml_gi_description(outputfile)
                self.change_version(outputfile)
        else:
            self.set_error("diamond运行出错!", code="31101202")
            self.set_error("diamond运行出错!", code="31101203")

    def run(self):
        """
        运行
        :return:
        """
        super(DiamondTool, self).run()
        if self.option("database") == 'customer_mode':
            self.run_makedb_and_diamond()
        else:
            self.run_diamond(self.option("database"))

    def get_nrxml_gi_description(self, xml_path):
        """
        从参考库NR_sequence中找到gi号对应的description
        """
        try:
            xml = ET.parse(xml_path)
            root = xml.getroot()
            self.logger.info("开始改注释信息内容!")
            BlastOutput_iterations = root.find('BlastOutput_iterations')
            for one_query in BlastOutput_iterations.findall('Iteration'):
                iteration_hits = one_query.findall('Iteration_hits')
                for iteration_hit in iteration_hits:
                    hits = iteration_hit.findall('Hit')
                    for hit in hits:
                        hit_id = hit.find('Hit_id')
                        hit_def = hit.find('Hit_def')
                        # hit_def = hit.find('Hit_id')
                        # hit_id = hit.find('Hit_def')
                        hit_id_split = re.split(r'\s', hit_id.text, maxsplit=1)
                        if self.option("database") in ['nr']:
                            gi_id_list = hit_id_split[0].split("|")
                            if len(gi_id_list) == 1:
                                hit_id_split = re.split(r'\s', hit_def.text, maxsplit=1)  # hit_ref gai hit_def
                                gi_id_list = hit_id_split[0].split("|")
                            gi_id = gi_id_list[1]
                            try:
                                int(gi_id)
                            except:
                                self.logger.info(hit_def)
                            result = self.mongodb_nr.find_one({"_id": int(gi_id)})
                            if result:
                                desc = result["desc"]
                                desc = desc.split("|")
                                description = desc[4]
                                hit_def.text = description
                            else:
                                hits.remove(hit)
                                print "没找到gid:{}".format(gi_id)
                        elif self.option("database") in ['nr_v20200604']:
                            seq_id = hit_id_split[0]
                            result = self.mongodb_nr.find_one({"origin_sequence_id": seq_id})
                            if result:
                                desc = result["desc"]
                                desc = desc.split("|")
                                description = desc[-1]
                                hit_def.text = description.strip()
                            else:
                                hits.remove(hit)
                                print "没找到gid:{}".format(seq_id)
                        elif self.option("database") in ['nr_v20210917']:
                            seq_id = hit_id_split[0]
                            result = self.mongodb_nr.find_one({"_id": seq_id})
                            if result:
                                desc = result["desc"]
                                desc = desc.split("|")
                                description = desc[-1]
                                hit_def.text = description.strip()
                            else:
                                hits.remove(hit)
                                print "没找到gid:{}".format(seq_id)
            xml.write('tmp.txt')
            with open('tmp.txt', 'rb') as f, open('tmp.xml', 'wb') as w:
                lines = f.readlines()
                a = '<?xml version=\"1.0\"?>\n<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n'
                w.write(a)
                w.writelines(lines)
            os.remove('tmp.txt')
            shutil.move('tmp.xml', xml_path)
        except (ServerSelectionTimeoutError, NetworkTimeout):  # 捕获因为mongo服务器问题导致的异常后重运行此方法
            if self.process_rerun < 5:
                self.process_rerun += 1
                self.logger.info("检测到TimeoutError, 第{}次重运行方法".format(self.process_rerun))
                time.sleep(5)
                self.get_nrxml_gi_description(xml_path)
            else:
                self.add_state('memory_limit', '检测到TimeoutError, 重运行tool')


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
                    if self.option("version") == "2":
                        line = line.replace("diamond 2.0.13", "BLASTX 2.3.0+")
                    else:
                        line = line.replace("diamond 0.8.35", "BLASTX 2.3.0+")
                if line.lstrip().startswith("<Hit_id>"):
                    m = re.match("<Hit_id>(.+)</Hit_id>", line.lstrip())
                    if m:
                        if self.repl:
                            line = line.replace(self.ori[i],self.repl[i])
                if line.lstrip().startswith("<Hit_def>"):
                    m = re.match("<Hit_def>(.+)</Hit_def>", line.lstrip())
                    if m:
                        if self.repl:
                            line = line.replace(self.repl[i],self.ori[i])
                        i += 1
                fw.write(line)
        # os.system("mv {} {}".format(path + "_new", path))
        os.remove(path)
        os.link(path + "_new", path)
        db_name = self.option("database")
        if db_name in ["nr", "animal", "fungi", "metazoa", "plant", "protist", "vertebrates", 'nr_v20200604',"nr_v20210917"]:
            if os.path.getsize(path) > 0:
                self.get_nrxml_gi_description(path)
        self.end()
