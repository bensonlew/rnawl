# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import shutil
from biocluster.core.exceptions import OptionError

class BlastAgent(Agent):
    """
    metaasv blast方法 blast+ 2.3.0
    功能不同之处：1.数据库全（不包含nt），2.coverage和identity进行筛选过滤，3.注释
    """
    def __init__(self, parent):
        self.DATABASE = ['unite7.2/its_fungi','unite8.0/its_fungi','fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS','fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA', 'fgr/amoA_archaea', 'fgr/amoA_bacteria','maarjam081/AM', 'Protist_PR2_v4.5','silva132/16s_archaea', 'silva132/16s_bacteria','silva132/18s_eukaryota', 'silva132/16s','silva138/16s_archaea', 'silva138/16s_bacteria','silva138/18s_eukaryota', 'silva138/16s','greengenes135/16s', 'greengenes135/16s_archaea', 'greengenes135/16s_bacteria','rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea', 'nt', 'nt/16s','nt/18s','nt/its','nt_v20200327/16s_archaea', 'nt_v20200327/16s_bacteria','nt_v20200327/16s','nt_v20210917/16s_archaea', 'nt_v20210917/16s_bacteria', 'nt_v20210917/16s',
                         'nt_v20210917/18s_eukaryota', 'nt_v20210917/its_fungi','nt_v20210917','Human_HOMD_v15.2', 'nt_v20200327/18s_eukaryota', 'nt_v20200327/its_fungi',"nt_v20200604",'fgr/amoA_archaea_202012', 'fgr/amoA_bacteria_202012', 'fgr/amoA_AOB_like_202012', 'fgr/amoA_comammox_202012', 'fgr/nosZ_202012', 'fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012', 'fgr/nirK_202012', 'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012', 'fgr/pmoA_202012', 'fgr/mmoX_202012']
        super(BlastAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "query_type", "type": "string"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "nr"},# 比对数据库多种，前面已经做了检查，后面不再做检查 custom_mode
            {"name": "blast", "type": "string", "default": "blastn"},  # 设定blast程序有blastn，blastx，blastp此处需要严格警告使用者必须选择正确的比对程序
            {"name": "ref", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择custom_mode时启用
            {"name": "reference_type", "type": "string"},  # 参考序列(库)的类型  为nucl或者prot
            {"name": "top_num", "type": "int", "default": 2},  # 序列比对最大输出条数，默认1条
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "identity", "type": "int", "default": 80},  # identity
            {"name": "coverage", "type": "int", "default": 80},  # coverage
            {"name": "num_threads", "type": "int", "default": 6},  # cpu数
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
        ]
        self.add_option(options)
        self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')))
        if self.option("database") == 'custom_mode':
            if not self.option("database"):
                raise OptionError("使用自定义数据库模式时必须设置reference")
            if self.option('reference_type') not in ['nucl', 'prot']:
                raise OptionError('reference_type参考序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')))
        elif self.option("database").lower() not in self.DATABASE and (self.option("database") not in self.DATABASE):
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
        self._cpu = self.option("top_num")
        self._memory = '40G'

    def end(self):
        super(BlastAgent, self).end()


class BlastTool(Tool):
    def __init__(self, config):
        self.DATABASE = ['unite7.2/its_fungi','unite8.0/its_fungi','fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS','fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA', 'fgr/amoA_archaea', 'fgr/amoA_bacteria','maarjam081/AM', 'Protist_PR2_v4.5','silva132/16s_archaea', 'silva132/16s_bacteria','silva132/18s_eukaryota', 'silva132/16s','silva138/16s_archaea', 'silva138/16s_bacteria','silva138/18s_eukaryota', 'silva138/16s','greengenes135/16s', 'greengenes135/16s_archaea', 'greengenes135/16s_bacteria','rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea', 'nt', 'nt/16s','nt/18s','nt/its','nt_v20200327/16s_archaea', 'nt_v20200327/16s_bacteria','nt_v20200327/16s','Human_HOMD_v15.2', 'nt_v20200327/18s_eukaryota', 'nt_v20200327/its_fungi',"nt_v20200604",'fgr/amoA_archaea_202012', 'fgr/amoA_bacteria_202012', 'fgr/amoA_AOB_like_202012', 'fgr/amoA_comammox_202012', 'fgr/nosZ_202012', 'fgr/nosZ_atypical_1_202012','fgr/nosZ_atypical_2_202012', 'fgr/nirK_202012', 'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012', 'fgr/pmoA_202012', 'fgr/mmoX_202012''nt_v20210917/16s_archaea', 'nt_v20210917/16s_bacteria', 'nt_v20210917/16s','nt_v20210917/18s_eukaryota', 'nt_v20210917/its_fungi']
        super(BlastTool, self).__init__(config)
        database_name = self.option('database').replace("/", "_") if re.search(r"/", self.option('database')) else self.option('database')
        if self.option("database") in self.DATABASE:
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/taxon_db/blastdb", database_name)
            # self.db_path = os.path.join("/mnt/ilustre/users/sanger-dev/home/zhangqingchen/metaasv/database", database_name)
        elif self.option("database") == "nt_v20210917":
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/align/ncbi/db/nt/nt_v20210917/nt")
        elif self.option("database") == 'custom_mode':
            db_name = os.path.splitext(os.path.basename(self.option("ref").prop['path']))[0]
            self.db_path = os.path.join(self.work_dir, 'custom_blastdb', db_name)
        self.cmd_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"  # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        # self.set_environ(BLASTDB=self.db_path)
        self.python = "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/metaasv/anno_taxonomy.py"
        self.taxon_db = os.path.join(self.config.SOFTWARE_DIR, "database/taxon_db/blastdb", database_name+ ".tax")
        # self.taxon_db = os.path.join("/mnt/ilustre/users/sanger-dev/home/zhangqingchen/metaasv/database", database_name+ ".tax")
        if self.option("blast") == "blastp":
            database_name = self.option('database').split("/")[1]
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/Framebot/fgr", database_name)
            self.taxon_db = os.path.join(self.config.SOFTWARE_DIR, "database/Framebot/tax", database_name + ".tax")

    def run_makedb_and_blast(self):
        """
        运行makeblastdb和blast
        :return:
        """
        if os.path.exists(os.path.join(self.work_dir, 'custom_blastdb')):
            shutil.rmtree(os.path.join(self.work_dir, 'custom_blastdb'))
        os.mkdir(os.path.join(self.work_dir, 'custom_blastdb'))
        cmd = os.path.join(self.cmd_path, "makeblastdb")
        cmd += " -dbtype %s -in %s -parse_seqids -out %s " % (self.option('reference_type'),
                                                              self.option("ref").prop['path'],
                                                              self.db_path)
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
        self.outputfile = os.path.join(self.output_dir, query_name + ".blast.m6.xls")
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt 6 -max_hsps 10 -num_threads %s -max_target_seqs %s -perc_identity %s -qcov_hsp_perc %s" % (
            self.option("query").prop['path'], self.db_path, self.outputfile,
            self.option("evalue"), self.option("num_threads"), self.option('top_num'), self.option('identity'),self.option('coverage'))
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
        获取注释信息
        """
        if os.path.getsize(self.outputfile) == 0:
            self.set_error("blast结果为空，请调整identity和coverage或者更换数据库！")
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        self.out = os.path.join(self.output_dir, "ASV_tax_assignments.txt")
        if self.option("database") in ["custom_mode"]:
            cmd = "{} {} -i {} -a {} -o {}".format(self.python, self.python_script, self.outputfile,
                                                    self.option("ref_taxon").prop['path'],self.out)
        else:
            cmd = "{} {} -i {} -a {} -o {}".format(self.python, self.python_script, self.outputfile,
                                                   self.taxon_db, self.out)
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
        if (len(os.listdir(self.output_dir)) != 0) and os.path.exists(self.output_dir + "/ASV_tax_assignments.txt"):
            self.end()
        else:
            if self.option("database") == 'custom_mode':
                self.run_makedb_and_blast()
                self.run_blast()
                self.get_taxid()
                self.end()
            else:
                self.run_blast()
                self.get_taxid()
                self.end()