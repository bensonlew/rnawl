# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from Bio import SeqIO
import pandas as pd
import unittest


class PpinetworkMapblastAgent(Agent):
    """
    blast比对蛋白互作数据库并且提取相关结果
    ncbi blast+ 2.3.0  注意：在outfmt为6时不按照ncbi格式生成table，而是按照特殊的表头类型，参见packages.align.blast.xml2table
    version 1.0
    author: shenghe
    last_modify: 2016.6.15
    """
    def __init__(self, parent):
        super(PpinetworkMapblastAgent, self).__init__(parent)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        self._blast_type = {'nucl': {'nucl': ['blastn', 'tblastn'],
                                     'prot': ['blastx']},
                            'prot': {'nucl': [],
                                     'prot': ['blastp']}}
        self._database_type = {'nt': 'nucl', 'nr': 'prot', 'kegg': 'prot', 'swissprot': 'prot', 'string': 'prot', 'archaea': 'prot', 'viruses': 'prot', 'fungi': 'prot', 'viridiplantae': 'prot', 'eukaryota_other': 'prot', 'eukaryota': 'prot', 'craniata': 'prot', 'bacteria': 'prot'}
        options = [
            {"name": "fa", "type": "infile", "format": "itraq_and_tmt.fasta"},  # 输入文件
            {"name": "query_type", "type": "string"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "customer_mode"},# 比对数据库 nt nr string swissprot kegg customer_mode 
            {"name": "outfmt", "type": "int", "default": 6},  # 输出格式，数字遵从blast+
            {"name": "mapped_table", "type": "string"},  # 输出格式，数字遵从blast+
            {"name": "unmapped_seq", "type": "string"},
            {"name": "unmapped_db", "type": "string"},
            {"name": "species", "type": "int"},
            {"name": "blast", "type": "string", "default": "blastp"},  # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "reference_type", "type": "string", "default":"prot"},  # 参考序列(库)的类型  为nucl或者prot
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "num_alignment", "type": "int", "default": 1},  # 序列比对最大输出条数，默认500
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            {"name": "outtable", "type": "outfile", "format": "itraq_and_tmt.common"},  # 输出格式为6时输出
            # 当输出格式为非5，6时，只产生文件不作为outfile
            ]
        self.add_option(options)
        self.step.add_steps('blast')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def step_start(self):
        self.step.blast.start()
        self.step.update()
        
    def step_end(self):
        self.step.blast.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fa").is_set:
            raise OptionError("必须设置参数fa", code = "32503801")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables = (self.option('query_type')), code = "32503802")
        else:
            if self._fasta_type[self.option('fa').prop['seq_type']] != self.option('query_type'):
                raise OptionError('文件检查发现查询序列为:%s, 而说明的文件类型为:%s', variables = (self._fasta_type[self.option('fa').prop['seq_type'], self.option('query_type')]), code = "32503803")
        if self.option("database") == 'customer_mode':
            pass
            # if not self.option("reference").is_set:
            #     raise OptionError("使用自定义数据库模式时必须设置reference")
            # if self.option('reference_type') not in ['nucl', 'prot']:
            #    raise OptionError('reference_type参考序列的类型为nucl(核酸)或者prot(蛋白):{}'.format(self.option('query_type')))
            # else:
            #     if self._fasta_type[self.option('reference').prop['seq_type']] != self.option('reference_type'):
            #         raise OptionError(
            #             '文件检查发现参考序列为:{}, 而说明的文件类型为:{}'.format(
            #                 self._fasta_type[self.option('reference').prop['seq_type'], self.option('reference_type')]))
        elif self.option("database").lower() not in ["nt", "nr", "string", 'kegg', 'swissprot', 'archaea', 'viruses', 'fungi', 'viridiplantae', 'eukaryota_other', 'eukaryota', 'craniata', 'bacteria']:
            raise OptionError("数据库%s不被支持", variables = (self.option("database")), code = "32503804")
        else:
            self.option('reference_type', self._database_type[self.option("database").lower()])
        if not 15 > self.option('outfmt') > -1:
            raise OptionError('outfmt遵循blast+输出规则，必须为0-14之间：%s', variables = (str(self.option('outfmt'))), code = "32503805")
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s', variables = (str(self.option('evalue'))), code = "32503806")
        if not 0 < self.option('num_alignment') < 1001:
            raise OptionError('序列比对保留数必须设置在1-1000之间:%s', variables = (str(self.option('num_alignment'))), code = "32503807")
        # if self.option('blast') not in self._blast_type[self.option('query_type')][self.option('reference_type')]:
        #     raise OptionError(
        #'程序不试用于提供的查询序列和库的类型，请仔细检查，核酸比对核酸库只能使用blastn或者tblastn，\
        # blastp， 蛋白比对蛋白库只能使用blastp, 或者没有提供blast参数')
        return True

    def set_resource(self):
        self._cpu = self.option('num_threads')
        self._memory = '20G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ])
        result_dir.add_regexp_rules([
            [r".+_vs_.+\.xml", "xml", "blast比对输出结果，xml格式"],
            [r".+_vs_.+\.xls", "xls", "blast比对输出结果，表格(制表符分隔)格式"],
            [r".+_vs_.+\.txt", "txt", "blast比对输出结果，非xml和表格(制表符分隔)格式"],
            [r".+_vs_.+\.txt_\d+\.xml", "xml", "Blast比对输出多xml结果，输出格式为14的单个比对结果文件,主结果文件在txt文件中"],
            [r".+_vs_.+\.txt_\d+\.json", "json", "Blast比输出对多json结果，输出格式为13的单个比对结果文件,主结果文件在txt文件中"],
            ])
        # print self.get_upload_files()
        super(PpinetworkMapblastAgent, self).end()


class PpinetworkMapblastTool(Tool):
    def __init__(self, config):
        super(PpinetworkMapblastTool, self).__init__(config)
        self._version = "2.3.0"
        if self.option("database") in ['nt', 'string', 'swissprot', 'kegg']:
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/align/ncbi/db")
        elif self.option("database") == "nr":   # add 6 lines by khl 20170217
            if not self.option("nr_species"):
               self.db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/align/ncbi/db/nr_db_20160623/nr')
            else:
               self.db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/align/ncbi/db/nr_db_20160623/{}'.format(self.option('nr_species')))
               self.logger.info(self.db_path)
        self.species_path = self.config.SOFTWARE_DIR + '/database/Annotation/all/String/string11.5/split/' + str(self.option("species")) + ".fa"
        self.cmd_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"   # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        # self.set_environ(BLASTDB=self.db_path)

    def choose_unmapped_seq(self):
        """
        选择没有比对到的序列
        """
        with open(self.option("unmapped_seq") , 'r') as f:
            unmapped_list = [line.strip() for line in  f.readlines()]
        seq_records = SeqIO.parse(self.option('fa').prop['path'], 'fasta')
        with open("unmapped_seq.fa" , 'w') as f:
            for seq_record in seq_records:
                seq_seq = seq_record.seq
                seq_name = seq_record.name
                if seq_name in unmapped_list:
                    f.write('>{}\n{}\n'.format(seq_name, seq_seq))


    def choose_unmapped_db(self):
        """
        选择没有比对到的数据库序列
        """
        with open(self.option("unmapped_db") , 'r') as f:
            unmapped_list = [line.strip() for line in  f.readlines()]
        seq_records = SeqIO.parse(self.species_path, 'fasta')
        with open("unmapped_db.fa" , 'w') as f:
            for seq_record in seq_records:
                seq_seq = seq_record.seq
                seq_name = seq_record.name
                if seq_name in unmapped_list:
                    f.write('>{}\n{}\n'.format(seq_name, seq_seq))

    def run_makedb_and_blast(self):
        """
        运行makeblastdb和blast

        :return:
        """
        db_name = "unmapped_db.fa"
        cmd = os.path.join(self.cmd_path, "makeblastdb")
        self.db_path = os.path.join(self.work_dir, 'customer_blastdb')
        cmd += " -dbtype %s -in %s -parse_seqids -out %s " % (self.option('reference_type'), 'unmapped_db.fa', 'unmapped_db')
        self.logger.info("开始运行makeblastdb，生成结果库文件放在工作目录的customer_blastdb下")
        makedb_obj = self.add_command("makeblastdb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("makeblastdb运行完成")
            self.run_blast()
        else:
            self.set_error("makeblastdb运行出错!", code = "32503808")

    def run_blast(self):
        """
        运行Blast

        :param db_name: blastdb名称
        :return:
        """
        # db = os.path.join(self.db_path, db_name)
        db_name = "unmapped_db"
        query_name = "unmapped_seq.fa"
        cmd = os.path.join(self.cmd_path, self.option('blast'))
        outputfile = os.path.join(self.output_dir, "unmapped_blast")
        outfmt = self.option('outfmt')

        cmd += " -query %s -db %s -out %s -evalue %s -outfmt %s -max_hsps 10 -num_threads %s -max_target_seqs %s" % (
            query_name, db_name, outputfile,
            self.option("evalue"), 6, self.option("num_threads"), self.option('num_alignment'))
        self.logger.info("开始运行blast")
        blast_command = self.add_command("blast", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.option('outtable', outputfile)
        else:
            self.set_error("blast运行出错!", code = "32503809")

    def get_blast_best(self):
        """
        获取blast比对按名称匹配和双向blast匹配最高的基因
        """
        output_xls = os.path.join(self.output_dir, "mapped_all.xls")
        gene_list = []
        gene_list2 = list()
        with open(output_xls , 'w') as f:
            # 写入根据序列名比对的结果
            with open(self.option("mapped_table") , 'r') as f_id:
                head = f_id.readline().strip()
                f.write(head + "\tmethod\n")
                for line in f_id.readlines():
                    if line.split("\t")[1] in gene_list:
                        # 判断根据已知ID对应的关系有没有多个accession对应到一个STRING ID的情况
                        pass
                    else:
                        gene_list.append(line.split("\t")[1])
                        f.write(line.strip() + "\tby_name\n")
            #写入根据blast匹配的结果
            blastfile = os.path.join(self.output_dir, "unmapped_blast")
            if os.path.getsize(blastfile) == 0:
                self.logger.info("比对结果为空，没有序列比对到数据库中")
                pass
            else:
                blast=pd.read_table(blastfile, header = None, dtype=str)
                all_hit=list(set(blast[1]))
                for hit in all_hit:
                    best_hit=blast[blast[1]==hit].sort_values(10, ascending=True).iloc[0]
                    if best_hit[1] in gene_list2:
                        pass
                    else:
                        gene_list2.append(best_hit[1])
                        f.write(str(best_hit[0]) + "\t" + str(best_hit[1]) + "\tby_sequence\n")

    def run(self):
        """
        运行
        :return:
        """
        super(PpinetworkMapblastTool, self).run()
        self.choose_unmapped_db()
        self.choose_unmapped_seq()
        self.run_makedb_and_blast()
        self.get_blast_best()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_itraq_and_tmt/data4'
        data = {
            "id": "mapblast" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "itraq_and_tmt.proteinset.ppinetwork_mapblast",
            "instant": True,
            "options": dict(
                fa = test_dir + "/" + "data4.fa",
                query_type = "prot",
                database = "customer_mode",
                mapped_table = test_dir + "/" + "diff_exp_mapped.txt",
                unmapped_seq =  test_dir + "/" + "unmapped_gene.txt",
                unmapped_db =  test_dir + "/" + "unmapped_db.txt",
                species = 4081,
                blast = "blastp"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
