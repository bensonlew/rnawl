# -*- coding: utf-8 -*-
# __author__ = 'Shijin,zoujiaxun'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
import shutil
import sqlite3
import xml.etree.ElementTree as ET
import cPickle as cpickle
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from mbio.packages.rna.annot_config import AnnotConfig

class DiamondAgent(Agent):
    """
    diamond version: 0.8.35
    version 1.0
    author: shijin
    last_modify: 20170316

    diamond version: 0.9.24
    last_modify: 20200110
    database: "database/Annotation/latest2019"
    """
    def __init__(self, parent):
        super(DiamondAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "prok_rna.common"},  # 输入文件
            {"name": "query_type", "type": "string", "default": "nucl"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "plant"},
            # 比对数据库 plant, nr, etc.
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，数字默认为5，输出xml
            {"name": "outxml", "type": "outfile", "format": "prok_rna.blast_xml"},
            {"name": "blast", "type": "string", "default": "blastp"},  # 设定diamond程序有blastp，blastx
            {"name": "reference", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考序列  选择customer时启用
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 20},  # cpu数
            {'name': 'kegg_version', 'type': 'string', 'default': "2019"},
            {'name': 'cog_version', 'type': 'string', 'default': "2020"},
            {"name": "nr_version", "type": "string", "default": "2019"},
            {"name": "swissprot_version", "type": "string", "default": "2019"},
            {"name": "string_version", "type": "string", "default": "2019"},
            {"name": "eggnog_version", "type": "string", "default": "2019"},
            {"name": "sensitive", "type": "int", "default": 2},
            {"name": "version", "type": "string", "default": "2019"},
            {"name": "diamond_version", "type": "string", "default": "v0.9.24.125"},
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
            raise OptionError("必须设置参数query", code = "35000401")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError("query_type查询序列的类型为nucl（核酸）或者prot（蛋白）：%s", variables = (self.option('query_type')), code = "32000402")
        if not 1 > self.option('evalue') >= 0:
            raise OptionError("E-value值设定必须为[0-1)之间：%s", variables = (self.option('evalue')), code = "35000403")
        if not 0 <= self.option("sensitive") <= 2:
            raise OptionError("敏感度设定必须为[0-2]之间：%s", variables = (self.option('evalue')), code = "35000404")
        return True

    def set_resource(self):
        self._cpu = self.option('num_threads')
        self._memory = '20G'
        self._memory_increase_step = 20
        if self.option("database") == "nr":
            self._memory = '50G'

    def end(self):
        super(DiamondAgent, self).end()


class DiamondTool(Tool):
    def __init__(self, config):
        super(DiamondTool, self).__init__(config)
        self._version = "0.9.24"
        self.get_db_path()
        # self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/Annotation/latest")
        '''
        self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/Annotation/latest2019")

        if self.option("database") == "kegg":
            version_path = "latest2019"
            if self.option("kegg_version") != "2019":
                version_path = 'other{}/kegg{}'.format(self.option("kegg_version")[0:4], self.option("kegg_version"))

        self.db_path = AnnotConfig().get_file_path(
            file=self.option("database"),
            version=self.option("version"),
            soft="diamond",
            
            )

        if self.option("database").startwith("kegg"):
            self.db_path = AnnotConfig().get_file_path(
                file=self.option("database"),
                version=self.option("version"),
                soft="diamond",
                soft_version=self.option("diamond_version"))

        if self.option("diamond_version") == "v0.9.24.125":
            self.cmd_path = "miniconda2/bin"
        else:
            self.cmd_path = "bioinfo/align/diamond-0.8.35" 

            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/Annotation/{}".format(version_path))
        # self.cmd_path = "bioinfo/align/diamond-0.8.35"   # 执行程序路径必须相对于 self.config.SOFTWARE_DIR

        '''
        '''
        if self.option("database").startwith("kegg"):
            self.db_path = AnnotConfig().get_file_path(
                file=self.option("database"),
                version=self.option("version"),
                soft="diamond",
                soft_version=self.option("diamond_version"))
        '''
        if self.option("diamond_version") == "v2.0.4":
            self.cmd_path = "bioinfo/rna/miniconda2_diamond/bin"
        if self.option("diamond_version") == "v2.0.13":
            self.cmd_path = "bioinfo/align/diamond-2.0.13"
        elif self.option("diamond_version") == "v0.9.24.125":
            self.cmd_path = "miniconda2/bin"
        else:
            self.cmd_path = "bioinfo/align/diamond-0.8.35" 

        if self.option("query_type") == "nucl":
            self.blast_type = "blastx"
        else:
            self.blast_type = "blastp"
        # self.mongodb_nr = Config().biodb_mongo_client.sanger_biodb.NR_sequence
        self.mongodb_nr = Config().get_mongo_client(mtype="ref_rna", ref=True)[Config().get_mongo_dbname("ref_rna", ref=True)].NR_sequence
        self.ori = []
        self.repl = []

    def get_db_path(self):
        sub_dict = dict()
        for sub_version in ['kegg_version', 'nr_version', 'swissprot_version', 'eggnog_version', 'string_version', 'cog_version']:
            sub_dict.update({sub_version: self.option(sub_version)})

        self.db_path = AnnotConfig().get_dmnd_path(db_name=self.option("database"),
                                                   version=self.option("version"),
                                                   **sub_dict)
        return self.db_path


    def run_makedb_and_diamond(self):
        """
        创建diamond数据库并运行diamond

        :return:
        """
        db_name = os.path.splitext(os.path.basename(self.option("reference").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, "diamond")
        self.db_path = self.work_dir
        cmd += " makedb --in {} --db {}".format(self.option("reference").prop['path'], db_name)
        self.logger.info("开始创建diamond数据库，生成结果库文件放在工作目录的customer_blastdb下")
        makedb_obj = self.add_command("makedb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("创建diamond数据库完成")
            self.run_diamond(db_name)
        else:
            self.set_error("创建diamond数据库出错", code = "35000405")

    def run_diamond(self, db_name):
        """
        运行diaomond

        :param db_name: blastdb名称
        :return:
        """
        hit_num = 5
        if self.option("database") == "kegg":
            hit_num = 1
        db = self.db_path
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, "diamond")
        outputfile = os.path.join(self.output_dir, query_name + "_vs_" + db_name)
        outfmt = self.option('outfmt')
        # if self.option('outfmt') == 5:
        outputfile += '.xml'  # outfmt默认为5
        outfmt = 5
        cmd += " {} -q {} -d {} -o {} -f {} -p {} -k {}".format(
            self.blast_type, self.option("query").prop['path'], db,
            outputfile, outfmt, self.option("num_threads"), hit_num)
        if self.option("evalue") != None:
            cmd += " -e {}".format(self.option("evalue"))

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
            if blast_command.return_code in [1]:
                self.add_state("memory_limit")
            else:
                self.set_error("diamond运行出错!")
        else:
            self.set_error("diamond运行出错", code = "35000406")

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


        if self.option("database") == "swissprot":
            nracc2des = os.path.dirname(self.db_path) + "/swissprot_acc2des.db"
            if os.path.exists(nracc2des):
                pass
            else:
                nracc2des = AnnotConfig().get_file_path(
                    file="swissprot_acc2des.db",
                    db="swissprot",
                    db_type="file",
                    version=self.option("swissprot_version"))
        else:
            nracc2des = os.path.dirname(self.db_path) + "/nr_acc2des.db"
            if os.path.exists(nracc2des):
                pass
            else:
                nracc2des = AnnotConfig().get_file_path(
                    file="nr_acc2des.db",
                    db="nr",
                    db_type="file",
                    version=self.option("nr_version"))

        conn = sqlite3.connect(nracc2des)
        cursor = conn.cursor()

        xml = ET.parse(xml_path)
        root = xml.getroot()
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
                    
                    acc_id = hit_id_split[0]
                    if len(acc_id.split("|")) > 1 and self.option("database") != "swissprot":
                        gi_id = hit_id_split[0].split("|")[1]
                        acc_id = hit_id_split[0].split("|")[-2]
                    # try:
                    #     int(gi_id)
                    # except:
                    #     self.logger.info(hit_def)
                    # nracc2des = self.db_path + "/nr_acc2des.cpickle"
                    # with open (nracc2des, 'rb') as f:
                    #    acc2des = cpickle.load(f)
                    desc = ""
                    try:
                        cursor.execute('select * from acc2des where acc="{}"'.format(acc_id))
                        desc = cursor.fetchall()[0][1]
                        description = desc
                        hit_def.text = description
                    except:
                        hits.remove(hit)
                        self.logger.info("没找到acc_id:{}".format(acc_id))


        xml.write('tmp.txt')
        with open('tmp.txt', 'rb') as f, open('tmp.xml', 'wb') as w:
            lines = f.readlines()
            a = '<?xml version=\"1.0\"?>\n<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n'
            w.write(a)
            w.writelines(lines)
        os.remove('tmp.txt')
        shutil.move('tmp.xml', xml_path)

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
                    line = line.replace("diamond 0.9.24", "BLASTX 2.9.0+")
                # if line.lstrip().startswith("<Hit_id>"):
                #     m = re.match("<Hit_id>(.+)</Hit_id>", line.lstrip())
                #     if m:
                #         line = line.replace(self.ori[i],self.repl[i])
                # if line.lstrip().startswith("<Hit_def>"):
                #     m = re.match("<Hit_def>(.+)</Hit_def>", line.lstrip())
                #     if m:
                #         line = line.replace(self.repl[i],self.ori[i])
                #         i += 1
                fw.write(line)
        # os.system("mv {} {}".format(path + "_new", path))
        os.remove(path)
        os.link(path + "_new", path)
        self.option('outxml',path)
        db_name = self.option("database")

        if db_name in ["nr", "animal", "fungi", "bacteria", "metazoa", "plant", "protist", "viridiplantae", "vertebrates", "swissprot"]:
            self.get_nrxml_gi_description(path)
            self.option('outxml',path + "_new")
        self.end()
