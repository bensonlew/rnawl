# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import unittest


class BlastForrelaAgent(Agent):
    """
    ncbi blast+ 2.3.0  注意：在outfmt为6时不按照ncbi格式生成table，而是按照特殊的表头类型，参见packages.align.blast.xml2table
    version 1.0
    author: shenghe
    last_modify: 2016.6.15
    """
    def __init__(self, parent):
        super(BlastForrelaAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "itraq_and_tmt.fasta"},  # 输入文件
            {"name": "query_type", "type": "string", "default": "prot"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "outfmt", "type": "int", "default": 6},  # 输出格式，数字遵从blast+
            {"name": "blast", "type": "string", "default": "blastp"},  # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "reference", "type": "infile", "format": "itraq_and_tmt.fasta"},  # 参考序列  选择customer时启用
            {"name": "reference_type", "type": "string", "default": "prot"},  # 参考序列(库)的类型  为nucl或者prot
            {"name": "evalue", "type": "float", "default": 1e-2},  # evalue值
            {"name": "num_threads", "type": "int", "default": 20},  # cpu数
            {"name": "num_alignment", "type": "int", "default": 50},  # 序列比对最大输出条数，默认500
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
            # 当输出格式为非5，6时，只产生文件不作为outfile
            ]
        self.add_option(options)
        self.step.add_steps('blast')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        # self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def step_start(self):
        self.step.blast.start()
        self.step.update()

    def step_end(self):
        self.step.blast.finish()
        self.step.update()

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = self.option('num_threads')
        self._memory = '50G'

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
        super(BlastForrelaAgent, self).end()


class BlastForrelaTool(Tool):
    def __init__(self, config):
        super(BlastForrelaTool, self).__init__(config)
        self._version = "2.3.0"
        self.cmd_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"   # 执行程序路径必须相对于 self.config.SOFTWARE_DIR
        self.filter_pep(self.option('query').prop['path'])
        self.filter_pep(self.option('reference').prop['path'])

    def run_makedb_and_blast(self):
        """
        运行makeblastdb和blast

        :return:
        """
        db_name = os.path.splitext(os.path.basename(self.option("reference").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, "makeblastdb")
        self.db_path = os.path.dirname(self.option('reference').prop['path'])
        cmd += " -dbtype %s -in %s -parse_seqids -out %s " % (self.option('reference_type'),
                                                              self.option("reference").prop['path'],
                                                              os.path.join(self.db_path, db_name))
        self.logger.info("开始运行makeblastdb，生成结果库文件放在工作目录的customer_blastdb下")
        makedb_obj = self.add_command("makeblastdb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("makeblastdb运行完成")
            self.run_blast(os.path.join(self.db_path, db_name))
        else:
            self.set_error("makeblastdb运行出错!")

    def run_blast(self, db_name):
        """
        运行Blast

        :param db_name: blastdb名称
        :return:
        """
        # db = os.path.join(self.db_path, db_name)
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        cmd = os.path.join(self.cmd_path, self.option('blast'))
        outputfile = os.path.join(self.output_dir, query_name + "_vs_" + os.path.basename(db_name))
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
            self.set_error("blast运行出错!")

    def run(self):
        """
        运行
        :return:
        """
        super(BlastForrelaTool, self).run()
        self.run_makedb_and_blast()

    def filter_pep(self, file):
        with open(file, 'r') as file_r:
            pep_info = file_r.read().split('\n>')
        with open(file, 'w') as file_w:
            for block in pep_info:
                block = block.strip().lstrip('>').split('\n')
                id = block[0]
                seq = '\n'.join(block[1:])
                if not u'.' in seq and id and seq:
                    file_w.write('>' + id + '\n' + seq + '\n')
                else:
                    seq = seq.replace('.', '*')
                    file_w.write('>' + id + '\n' + seq + '\n')

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript'
        data = {
            "id": "Blast_rel_" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "protein_transcript.blast_forrela",
            "instant": False,
            "options": dict(
                # query = test_dir + "/" + "exp.fasta",
                reference = test_dir + "/" + "exp.fasta",
                query = test_dir + "/" + "GCF_000826755.1_ZizJuj_1.1_protein.faa",
                # reference = test_dir + "/" + "GCF_000826755.1_ZizJuj_1.1_protein.faa",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
