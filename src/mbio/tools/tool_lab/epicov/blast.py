# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import ConfigParser
import unittest
from Bio.Blast import NCBIXML


class BlastAgent(Agent):
    """
    已知miRNA鉴定
    """
    def __init__(self, parent):
        super(BlastAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "small_rna.fasta"},  # 输入序列文件
            {"name": "query_type", "type": "string", "default": "nucl"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，数字遵从blast+
            {"name": "blast", "type": "string", "default": "blastn"},  # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "evalue", "type": "float", "default": 30000},
            {"name": "num_threads", "type": "int", "default": 10},
            {"name": "database", "type": "infile", "format": "small_rna.fasta"}, # 序列比对文件
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            {"name": "outtable", "type": "outfile", "format": "small_rna.common"},  # 输出格式为5时输出
            {"name": "word_size", "type": "int", "default": 7},
            {"name": "max_hsps", "type": "int", "default": 10},  # Max targets per sequence
            {"name": "max_target_seqs", "type": "int", "default": 1000},  # Max targets to show
        ]
        self.add_option(options)
        self.step.add_steps("blast")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.blast.start()
        self.step.update()

    def stepfinish(self):
        self.step.blast.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("query").is_set:
            raise OptionError("必须提供输入文件")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 11
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "比对结果目录"]
            ])
        result_dir.add_regexp_rules([
            ])
        super(BlastAgent, self).end()


class BlastTool(Tool):
    def __init__(self, config):
        super(BlastTool, self).__init__(config)
        self.blast_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"
        self.set_environ(PATH=self.blast_path)

    def run(self):
        """
        运行
        :return:
        """
        super(BlastTool, self).run()
        self.blast()
        blast_table = self.output_dir + "/" + os.path.basename(self.xml_fp).replace("xml", 'txt')
        self.xml2table(self.xml_fp, blast_table)
        self.set_output()
        self.end()

    def blast(self):
        """
        nt库比对
        """
        self.logger.info("开始进行nt库比对")
        input_file = self.option("query").prop["path"]
        cmd = os.path.join(self.blast_path, self.option('blast'))
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        self.xml_fp = os.path.join(self.output_dir, query_name + "_vs_db.xml")
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt %s -max_hsps %s -num_threads %s -max_target_seqs %s -word_size %s" % (
            input_file, self.option('database').prop['path'], self.xml_fp,
            self.option('evalue'), self.option('outfmt'), self.option('max_hsps'), self.option("num_threads"), self.option('max_target_seqs'), self.option('word_size'))
        command = self.add_command("blast_db", cmd).run()
        self.wait()
        blast_xml = self.output_dir + "/" + os.path.basename(self.xml_fp)
        self.option("outxml", blast_xml)
        if command.return_code == 0:
            self.logger.info("运行比对完成")
        else:
            self.set_error("运行比对出错")

    def xml2table(self, xml_fp, table_out, header=None, anno_head=False, hit_id=None):
        """

        :param xml_fp: 输入xml文件路径
        :param table_out: 输出文件路径
        :param header: 选择列，列的选择必须是上面定义的 all_values中的值组成的列表
        :param anno_head: 是否写入表头
        :param hit_id: 这里为兼容新版本的nt库的代码用，表示含义为是hit_id的取的方式改为从des中取
        """
        default_header = ['Query-Name', 'Hit-Name', 'Identity-%', 'HSP-Len', 'Mismatch', 'Gapopen_num', 'Q-Begin',
                          'Q-End', 'Hsp-Begin', 'Hsp-End', 'Expect', 'Score', 'Hsp_qseq', 'Hsp_hseq', 'Match']

        all_values = ['Score', 'E-Value', 'HSP-Len', 'Identity-%', 'Similarity-%', 'Identity', 'Positives',
                      'Query-Name', 'Q-Len', 'Q-Begin', 'Q-End', 'Q-Frame', 'Hit-Name', 'Hit-Len', 'Hsp-Begin',
                      'Hsp-End', 'Hsp-Frame', 'Hit-Description', 'Q-Strand', 'Hsp-Strand', 'Mismatch', 'Gapopen_num',
                      'Hsp_qseq', 'Hsp_hseq', 'Match', 'Expect']
        if header:
            for i in header:
                if i not in all_values:
                    raise Exception('无法获取的值:{}\n可用的值:{}'.format(i, '\t'.join(all_values)))
        else:
            header = default_header
        if not os.path.isfile(xml_fp):
            raise Exception('输入xml文件不存在:{}'.format(xml_fp))
        with open(xml_fp) as f, open(table_out, 'w') as w:
            if anno_head:
                w.write('\t'.join(header) + '\n')
            records = NCBIXML.parse(f)
            values = {i: 'N/A' for i in all_values}
            for rec in records:
                query = re.split(' ', rec.query, maxsplit=1)[0]
                for align in rec.alignments:
                    for hsp in align.hsps:
                        one_hsp = values.copy()
                        one_hsp['Query-Name'] = query
                        if hit_id:
                            one_hsp['Hit-Name'] = (align.hit_def).split(" ")[0]
                            one_hsp['Hit-Description'] = " ".join(align.hit_def.split(" ")[1:])
                        else:
                            one_hsp['Hit-Name'] = align.hit_id
                            one_hsp['Hit-Description'] = align.hit_def
                        one_hsp['Score'] = str(hsp.score)
                        one_hsp['E-Value'] = str(hsp.expect)
                        one_hsp['HSP-Len'] = str(hsp.align_length)
                        one_hsp['Identity'] = str(hsp.identities)
                        one_hsp['Positives'] = str(hsp.positives)
                        one_hsp['Q-Len'] = str(rec.query_length)
                        one_hsp['Q-Begin'] = str(hsp.query_start)
                        one_hsp['Q-End'] = str(hsp.query_end)
                        one_hsp['Q-Frame'] = str(hsp.frame[0])
                        one_hsp['Hit-Len'] = str(align.length)
                        one_hsp['Hsp-Begin'] = str(hsp.sbjct_start)
                        one_hsp['Hsp-End'] = str(hsp.sbjct_end)
                        one_hsp['Hsp-Frame'] = str(hsp.frame[1])
                        one_hsp['Q-Strand'] = str(hsp.strand[0])
                        one_hsp['Hsp-Strand'] = str(hsp.strand[1])
                        one_hsp['Hsp_qseq'] = str(hsp.query)
                        one_hsp['Hsp_hseq'] = str(hsp.sbjct)
                        one_hsp['Match'] = str(hsp.match)
                        one_hsp['Expect'] = str(hsp.expect)
                        one_hsp['Identity-%'] = str(round(float(hsp.identities) / hsp.align_length, 3) * 100)
                        one_hsp['Similarity-%'] = str(round(float(hsp.positives) / hsp.align_length, 3) * 100)
                        one_hsp['Mismatch'] = str(
                            int(hsp.align_length) - int(len(hsp.match)))  # 此处使用len(hsp.match) 还是hsp.indentities没有具体证据
                        one_hsp['Gapopen_num'] = str(hsp.gaps)  # 此处的gaps不是 gapopen 而是 gaps总数，因为xml获取不到
                        line = list()
                        for i in header:
                            line.append(one_hsp[i])
                        w.write('\t'.join(line) + '\n')
                        # line = ""
                        # for i in header:
                        #     line += one_hsp[i] + '\t'
                        # w.write(line.strip() + '\n')
        return table_out

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("开始设置结果目录")
        self.logger.info("设置结果完成")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "BlastRfam" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna_v2.srna.blast",
            "instant": False,
            "options": dict(
                query="/mnt/lustre/users/sanger/sg-users/shicaiping/primer_blast/primer.fa",
                database="nt"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()