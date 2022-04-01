# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class BlastAgent(Agent):
    """
    ncbi blast+ 2.3.0  注意：在outfmt为6时不按照ncbi格式生成table，而是按照特殊的表头类型，参见packages.align.blast.xml2table
    version 1.0
    author: shenghe
    last_modify: 2016.6.15
    """
    def __init__(self, parent):
        super(BlastAgent, self).__init__(parent)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        self._blast_type = {'nucl': {'nucl': ['blastn', 'tblastn'],
                                     'prot': ['blastx']},
                            'prot': {'nucl': [],
                                     'prot': ['blastp']}}
        self._database_type = {'nt': 'nucl', 'nr': 'prot', 'kegg': 'prot', 'swissprot': 'prot', 'string': 'prot', 'archaea': 'prot', 'viruses': 'prot', 'fungi': 'prot', 'viridiplantae': 'prot', 'eukaryota_other': 'prot', 'eukaryota': 'prot', 'craniata': 'prot', 'bacteria': 'prot'}
        options = [
            {"name": "query", "type": "infile", "format": "denovo_rna_v2.trinity_fasta"},  # 输入文件
            {"name": "query_type", "type": "string"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "nr"},
            # 比对数据库 nt nr string swissprot kegg customer_mode 
            {"name": "nr_species", "type": "string"}, # nr分类物种模式：Archaea Viruses Fungi Viridiplantae Eukaryota_other Eukaryota Craniata Bacteria
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，数字遵从blast+
            {"name": "blast", "type": "string"},  # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "reference", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "reference_type", "type": "string"},  # 参考序列(库)的类型  为nucl或者prot
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "num_alignment", "type": "int", "default": 5},  # 序列比对最大输出条数，默认500
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
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
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code = "32000101")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables = (self.option('query_type')), code = "32000102")
        else:
            if self._fasta_type[self.option('query').prop['seq_type']] != self.option('query_type'):
                raise OptionError(
                    '文件检查发现查询序列为:%s, 而说明的文件类型为:%s', variables = (self._fasta_type[self.option('query').prop['seq_type']
                    , self.option('query_type')]), code = "32000103")
        if self.option("database") == 'customer_mode':
            if not self.option("reference").is_set:
                raise OptionError("使用自定义数据库模式时必须设置reference", code = "32000104")
            if self.option('reference_type') not in ['nucl', 'prot']:
                raise OptionError('reference_type参考序列的类型为nucl(核酸)或者prot(蛋白):%s', variables = (self.option('query_type')), code = "32000105")
            else:
                if self._fasta_type[self.option('reference').prop['seq_type']] != self.option('reference_type'):
                    raise OptionError(
                        '文件检查发现参考序列为:%s, 而说明的文件类型为:%s', variables = (
                            self._fasta_type[self.option('reference').prop['seq_type'], self.option('reference_type')]), code = "32000106")
        elif self.option("database").lower() not in ["nt", "nr", "string", 'kegg', 'swissprot', 'archaea', 'viruses', 'fungi', 'viridiplantae', 'eukaryota_other', 'eukaryota', 'craniata', 'bacteria']:
            raise OptionError("数据库%s不被支持", variables = (self.option("database")), code = "32000107")
        else:
            self.option('reference_type', self._database_type[self.option("database").lower()])
        if not 15 > self.option('outfmt') > -1:
            raise OptionError('outfmt遵循blast+输出规则，必须为0-14之间：%s',variables = (self.option('outfmt')), code = "32000108")
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s', variables = (self.option('evalue')), code = "32000109")
        if not 0 < self.option('num_alignment') < 1001:
            raise OptionError('序列比对保留数必须设置在1-1000之间:%s', variables = (self.option('num_alignment')), code = "32000110")
        if self.option('blast') not in self._blast_type[self.option('query_type')][self.option('reference_type')]:
            raise OptionError(
                '程序不试用于提供的查询序列和库的类型，请仔细检查，核酸比对核酸库只能使用blastn或者tblastn，\
                 核酸比对蛋白库只能使用blastp， 蛋白比对蛋白库只能使用blastp, 或者没有提供blast参数', code = "32000111")
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
        super(BlastAgent, self).end()


class BlastTool(Tool):
    def __init__(self, config):
        super(BlastTool, self).__init__(config)
        self._version = "2.3.0"
        if self.option("database") in ['nt', 'string', 'swissprot', 'kegg']:
            self.db_path = os.path.join(self.config.SOFTWARE_DIR, "database/align/ncbi/db")
        elif self.option("database") == "nr":   # add 6 lines by khl 20170217
            if not self.option("nr_species"):
               self.db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/align/ncbi/db/nr_db_20160623/nr')
            else:
               self.db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/align/ncbi/db/nr_db_20160623/{}'.format(self.option('nr_species')))
               self.logger.info(self.db_path)
        self.logger.info(self.option("nr_species"))
        self.cmd_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"   # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        self.set_environ(BLASTDB=self.db_path)

    def run_makedb_and_blast(self):
        """
        运行makeblastdb和blast

        :return:
        """
        db_name = os.path.splitext(os.path.basename(self.option("reference").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, "makeblastdb")
        self.db_path = os.path.join(self.work_dir, 'customer_blastdb')
        cmd += " -dbtype %s -in %s -parse_seqids -out %s " % (self.option('reference_type'),
                                                              self.option("reference").prop['path'],
                                                              os.path.join(self.db_path, db_name))
        self.logger.info("开始运行makeblastdb，生成结果库文件放在工作目录的customer_blastdb下")
        makedb_obj = self.add_command("makeblastdb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("makeblastdb运行完成")
            self.run_blast(db_name)
        else:
            self.set_error("makeblastdb运行出错!", code = "32000112")

    def run_blast(self, db_name):
        """
        运行Blast

        :param db_name: blastdb名称
        :return:
        """
        # db = os.path.join(self.db_path, db_name)
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, self.option('blast'))
        outputfile = os.path.join(self.output_dir, query_name + "_vs_" + db_name)
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
            self.option("query").prop['path'], db_name, outputfile,
            self.option("evalue"), outfmt, self.option("num_threads"), self.option('num_alignment'))
        self.logger.info("开始运行blast")
        blast_command = self.add_command("blast", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.logger.info("运行blast完成")
            if self.option('outfmt') == 6:
                self.logger.info('程序输出结果为6(table)，实际需要结果为5(xml)，开始调用程序xml2table转换')
                from mbio.packages.align.blast.xml2table import xml2table
                xml2table(outputfile, outputfile[:-3] + 'xls')
                blast_xml = self.work_dir + "/" + os.path.basename(outputfile)
                os.system("mv {} {}".format(outputfile, blast_xml))
                self.option("outxml", blast_xml)
                self.logger.info('程序table格式转换完成，旧xml文件已移除')
                self.option('outtable', outputfile[:-3] + 'xls')
            elif self.option('outfmt') == 5:
                self.logger.info(outputfile)
                self.option('outxml', outputfile)
            self.end()
        else:
            self.set_error("blast运行出错!", code = "32000113")

    def run(self):
        """
        运行
        :return:
        """
        super(BlastTool, self).run()
        if self.option("database") == 'customer_mode':
            self.run_makedb_and_blast()
        else:
            if self.option("database")=="nr":
                 if not self.option("nr_species"):
                      self.run_blast(self.option("database"))
                 else:
                      self.run_blast(self.option("nr_species"))
            else:
                   self.run_blast(self.option("database"))

