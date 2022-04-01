# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify: 2018.01.29

from biocluster.agent import Agent
from biocluster.tool import Tool
import os, re
from biocluster.core.exceptions import OptionError


class RefalignDnabacAgent(Agent):
    """
    功能：参考基因组比对注释
    输入：基因序列
    参考库： 整理的参考数据库，放在指定文件夹下（默认为某一物种genebank号），或者提供genebank文件
    比对软件： ncbi blast+ 2.3.0
    """

    def __init__(self, parent):
        super(RefalignDnabacAgent, self).__init__(parent)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "query_type", "type": "string", "default": "prot"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "refalign_database", "type": "string", "default": ""},  # 整理的参考基因组对应Genbank号
            {"name": "outfmt", "type": "int", "default": 6},  # 输出格式，数字遵从blast+
            {"name": "blast", "type": "string", "default": "blastp"},
            # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "reference_gbk", "type": "infile", "format": "gene_structure.gbk"},  # 参考genebank文件
            {"name": "reference_seq", "type": "infile", "format": "sequence.fasta"},  # 参考fasta文件只允许prot格式
            {"name": "reference_type", "type": "string", "default": "prot"},  # 参考序列的格式，为nucl或者prot
            {"name": "evalue", "type": "float", "default": 10},  # evalue值
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "num_alignment", "type": "int", "default": 5},  # 序列比对最大输出条数，默认5
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
            # 当输出格式为非5，6时，只产生文件不作为outfile
            {"name": "coverage","type":"string","default":'F'},  # False : 结果中没有coverage，True 结果中有coverage
            {"name": "analysis_type", "type":"string", "default":""}, #paralogs
        ]
        self.add_option(options)
        self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="31101901")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')), code="31101902")
        else:
            if self._fasta_type[self.option('query').prop['seq_type']] != self.option('query_type'):
                raise OptionError('文件检查发现查询序列为:%s, 而说明的文件类型为:%s', variables=(self._fasta_type[self.option('query').prop['seq_type']], self.option('query_type')),code="31101903")
        if self.option('reference_type') not in ['nucl', 'prot']:
            raise OptionError('reference_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('reference_type')), code="31101904")
        if self.option("refalign_database") == '' and not self.option("reference_gbk").is_set and not self.option(
                "reference_seq").is_set:
            raise OptionError("使用自定义数据库模式时必须设置refalign_database或者提供参考基因组genebank文件", code="31101905")
        if self.option("reference_gbk").is_set:
            self.option("reference_gbk").check_cds_info()
        if not 15 > self.option('outfmt') > -1:
            raise OptionError('outfmt遵循RefalignDnabac+输出规则，必须为0-14之间：%s', variables=(self.option('outfmt')), code="31101906")
        if not 10 >= self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-10]之间：%s', variables=(self.option('evalue')), code="31101907")
        if not 0 < self.option('num_alignment') < 1001:
            raise OptionError('序列比对保留数必须设置在1-1000之间:%s', variables=(self.option('num_alignment')), code="31101908")
        return True

    def set_resource(self):
        self._cpu = self.option('num_threads')
        self._memory = '5G'

    def end(self):
        super(RefalignDnabacAgent, self).end()


class RefalignDnabacTool(Tool):
    def __init__(self, config):
        super(RefalignDnabacTool, self).__init__(config)
        self._version = "2.3.0"
        self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/Genome_prokaryote/")
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.get_seq_by_gbk = self.config.PACKAGE_DIR + "/sequence/scripts/get_seq_by_gbk.pl"  # 从上传的基因组genebank文件中获取基因组参考序列
        self.cmd_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"  # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        self.set_environ(RefalignDnabacDB=self.db_path)
        self.ref_sequence = ''

    def run_makedb_and_blast(self):
        """
        运行makeblastdb和blast
        :return:
        """
        if self.option("reference_gbk").is_set:
            self.ref_sequence = self.option("reference_gbk").prop['genome_name']
            if "," in self.ref_sequence:
                self.ref_sequence = self.ref_sequence.split(",")[0]
            print self.ref_sequence
            cmd_gbk = "%s %s %s %s" % (
                self.perl_path, self.get_seq_by_gbk, self.option("reference_gbk").prop['path'], self.ref_sequence)
            self.logger.info("从genebank文件中提取基因组序列")
            command = self.add_command("get_seq_by_gbk", cmd_gbk).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("get_seq_by_gbk运行完成")
            else:
                self.set_error("get_seq_by_gbk运行出错!", code="31101901")
        elif self.option("reference_seq").is_set:
            self.ref_sequence = self.option("reference_seq").prop['path']
        ref_sequence_name = os.path.splitext(os.path.basename(self.ref_sequence))[0]
        cmd = os.path.join(self.cmd_path, "makeblastdb")
        self.db_path = os.path.join(self.work_dir, 'customer_blastdb')
        cmd += " -dbtype %s -in %s -parse_seqids -out %s " % (self.option("reference_type"),
                                                              self.ref_sequence,
                                                              os.path.join(self.db_path, ref_sequence_name))
        self.logger.info("开始运行makeblastdb，生成结果库文件放在工作目录的customer_blastdb下")
        makedb_obj = self.add_command("makeblastdb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("makeblastdb运行完成")
            self.run_blast(ref_sequence_name)
        else:
            self.set_error("makeblastdb运行出错!", code="31101902")

    def run_blast(self, db_name):
        """
        运行Blast
        :param db_name: blastdb名称
        :return:
        """
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, self.option('blast'))
        outputfile = os.path.join(self.output_dir, query_name + "_vs_" + re.subn("\.faa$", "", db_name)[0])
        outfmt = self.option('outfmt')
        if self.option('outfmt') == 5:
            outputfile += '.xml'
        elif self.option('outfmt') == 6:
            # outputfile += '.xls'
            outputfile += '.xml'
            outfmt = 5  # 为了保证table格式输出表头完全一致，输出为6时，选用5xml为输出结果，后面再通过统一的xml2table转换
        else:
            outputfile += '.txt'
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt %s -max_hsps 10 -num_threads %s -max_target_seqs %s" % (
            self.option("query").prop['path'], os.path.join(self.db_path, db_name), outputfile,
            self.option("evalue"), outfmt, self.option("num_threads"), self.option('num_alignment'))
        self.logger.info("开始运行blast")
        blast_command = self.add_command("blast", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.logger.info("运行blast完成")
            if self.option('outfmt') == 6:
                self.logger.info('程序输出结果为6(table)，实际需要结果为5(xml)，开始调用程序xml2table转换')
                if self.option('coverage')!='F':   #zouguanqing 20190328
                    from mbio.packages.align.blast.xml2table import xml2table_coverage
                    if self.option("analysis_type") == "paralog":
                        xml2table_coverage(outputfile, outputfile[:-3] + 'xls',hit_ref=False, anno_head=True)
                    else:
                        xml2table_coverage(outputfile, outputfile[:-3] + 'xls')
                else:
                    from mbio.packages.align.blast.xml2table import xml2table
                    xml2table(outputfile, outputfile[:-3] + 'xls')
                    self.option('outtable', outputfile[:-3] + 'xls')
                blast_xml = self.work_dir + "/" + os.path.basename(outputfile)
                os.system("mv {} {}".format(outputfile, blast_xml))
                self.option("outxml", blast_xml)
                self.logger.info('程序table格式转换完成，旧xml文件已移除')

            elif self.option('outfmt') == 5:
                self.logger.info(outputfile)
                self.option('outxml', outputfile)
            self.end()
        else:
            self.set_error("blast运行出错!", code="31101903")

    def run(self):
        """
        运行
        :return:
        """
        super(RefalignDnabacTool, self).run()
        if self.option("reference_gbk").is_set or self.option("reference_seq").is_set:
            self.run_makedb_and_blast()
        elif self.option("refalign_database") != "":
            self.run_blast(self.option("refalign_database"))
