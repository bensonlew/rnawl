# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.lnc_rna import xml2lncannot
import unittest

class LncrnaAnnotAgent(Agent):
    """
    ncbi blast+ 2.3.0  注意：在outfmt为6时不按照ncbi格式生成table，而是按照特殊的表头类型，参见packages.align.blast.xml2table
    version 1.0
    author: shenghe
    last_modify: 2016.6.15
    """
    def __init__(self, parent):
        super(LncrnaAnnotAgent, self).__init__(parent)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        self._blast_type = {'nucl': {'nucl': ['blastn', 'tblastn'],
                                     'prot': ['blastx']},
                            'prot': {'nucl': [],
                                     'prot': ['blastp']}}
        self._database_type = {'nt': 'nucl', 'nr': 'prot', 'kegg': 'prot', 'swissprot': 'prot', 'string': 'prot', 'archaea': 'prot', 'viruses': 'prot', 'fungi': 'prot', 'viridiplantae': 'prot', 'eukaryota_other': 'prot', 'eukaryota': 'prot', 'craniata': 'prot', 'bacteria': 'prot'}
        options = [
            {"name": "query", "type": "infile", "format": "lnc_rna.fasta"},
            {"name": "query_type", "type": "string", "default": "nucl"},
            {"name": "database", "type": "string", "default": "customer_mode"},

            {"name": "outfmt", "type": "int", "default": 6},
            {"name": "blast", "type": "string", "default": "blastn"},
            {"name": "reference", "type": "infile", "format": "lnc_rna.fasta"},
            {"name": "reference_type", "type": "string", "default": "nucl"},
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "num_threads", "type": "int", "default": 10},
            {"name": "num_alignment", "type": "int", "default": 5},
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},
            
            {"name": "outtable", "type": "outfile", "format": "lnc_rna.common"},
            {"name": "qcov", "type": "float", "default": 80},
            {"name": "scov", "type": "float", "default": 80},
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
        if self.option("database") == 'customer_mode':
            if not self.option("reference").is_set:
                raise OptionError("使用自定义数据库模式时必须设置reference", code = "32000104")
            if self.option('reference_type') not in ['nucl']:
                raise OptionError('reference_type参考序列的类型为nucl(核酸)或者prot(蛋白):%s', variables = (self.option('query_type')), code = "32000105")
        else:
            raise OptionError("必须设置为custom mode")


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
            [r".+_vs_.+\.txt_\d+\.xml", "xml", "LncrnaAnnot比对输出多xml结果，输出格式为14的单个比对结果文件,主结果文件在txt文件中"],
            [r".+_vs_.+\.txt_\d+\.json", "json", "LncrnaAnnot比输出对多json结果，输出格式为13的单个比对结果文件,主结果文件在txt文件中"],
            ])
        # print self.get_upload_files()
        super(LncrnaAnnotAgent, self).end()


class LncrnaAnnotTool(Tool):
    def __init__(self, config):
        super(LncrnaAnnotTool, self).__init__(config)
        self._version = "2.3.0"
        self.cmd_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"   # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        self.python = 'miniconda2/bin/python'
        self.xml2table = self.config.PACKAGE_DIR + '/lnc_rna/xml2lncannot.py'

    def run_makedb_and_blast(self):
        """
        运行makeblastdb和blast

        :return:
        """
        db_name = os.path.splitext(os.path.basename(self.option("reference").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, "makeblastdb")
        self.db_path = os.path.join(self.work_dir, 'customer_blastdb')
        self.db = os.path.join(self.db_path, db_name)
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
        运行LncrnaAnnot

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
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt %s -max_hsps 5 -num_threads %s -max_target_seqs %s" % (
            self.option("query").prop['path'], self.db, outputfile,
            self.option("evalue"), outfmt, self.option("num_threads"), self.option('num_alignment'))
        self.logger.info("开始运行blast")
        blast_command = self.add_command("blast", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.logger.info("运行blast完成")
            if self.option('outfmt') == 6:
                self.logger.info('程序输出结果为6(table)，实际需要结果为5(xml)，开始调用程序xml2table转换')
                cmd = "{} {} -i {} -o {} -q {} -s {}".format(self.python, self.xml2table, outputfile, outputfile[:-3] + 'xls', self.option("qcov"), self.option("scov"))
                self.logger.info("开始运行blast2table")
                blast_command = self.add_command("blast2xls", cmd)
                blast_command.run()
                self.wait()

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
        super(LncrnaAnnotTool, self).run()
        if self.option("database") == 'customer_mode':
            self.run_makedb_and_blast()
        else:
            pass


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/lnc_rna/diff_exp'

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "lnc_annot" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "whole_transcriptome.lncrna_identification.lncrna_annot",
            "instant": False,
            "options": dict(
                reference="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/lncrna.fa",
                query="/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/LargeGush/output/new_lncrna_predict/novel_lncrna.fa",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()




