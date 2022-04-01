# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import unittest
import re
from collections import defaultdict


class ExtractGffAgent(Agent):
    def __init__(self, parent):
        super(ExtractGffAgent, self).__init__(parent)
        options = [
            {"name": "type", "type": "string", "default": "gff"},  # 输入文件的类型
            {"name": "gff_gtf", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
            {"name": "outlist", "type": "outfile", "format": "itraq_and_tmt.common"},  # 输出的有关联关系的list文件
            ]
        self.add_option(options)
        self.step.add_steps('gff')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.gff.start()
        self.step.update()

    def step_end(self):
        self.step.gff.finish()
        self.step.update()

    def check_options(self):
        if not os.path.exists(self.option('gff_gtf').prop['path']):
            raise OptionError("没有正确传入gff或者gtf文件")
        if not self.option('type') in ['gff', 'gtf']:
            raise OptionError("只有gff和gtf这两种情况")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(ExtractGffAgent, self).end()


class ExtractGffTool(Tool):
    def __init__(self, config):
        super(ExtractGffTool, self).__init__(config)

    def extract_gff(self):
        g2t = defaultdict(set)
        g2p = defaultdict(set)
        with open(self.option('gff_gtf').prop['path'], 'r') as gff:
            for line in gff:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    if len(line) > 4:
                        if line[2] == 'mRNA' or line[2] == 'exon':
                            gene = trans = ''
                            for att in line[8].split(';'):
                                if re.match('Dbxref=GeneID:(.*?),', att):
                                    gene = re.match('Dbxref=GeneID:(.*?),', att).group(1)
                                if re.match('^transcript_id=(.*?)$', att):
                                    trans = re.match('^transcript_id=(.*?)$', att).group(1)
                            if gene and trans:
                                g2t[gene].add(trans)
                        if line[2] == 'CDS':
                            gene = pro = ''
                            for att in line[8].split(';'):
                                if re.match('Dbxref=GeneID:(.*?),', att):
                                    gene = re.match('Dbxref=GeneID:(.*?),', att).group(1)
                                if re.match('^protein_id=(.*?)$', att):
                                    pro = re.match('^protein_id=(.*?)$', att).group(1)
                            if gene and pro:
                                g2p[gene].add(pro)
        return g2t, g2p

    def extract_gtf(self):
        g2t = dict()
        g2p = dict()
        return g2t, g2p

    def run(self):
        super(ExtractGffTool, self).run()
        if self.option('type') == 'gff':
            g2t, g2p = self.extract_gff()
        else:
            g2t, g2p = self.extract_gtf()
        with open(self.output_dir + '/g2t2p.list', 'w') as list_w:
            genes = sorted(list(set(g2t.keys() + g2p.keys())))
            for gene in genes:
                str_i = gene + '\t'
                # if gene in g2t:
                #     str_i += ';'.join(g2t[gene]) + '\t'
                # else:
                #     str_i += '_' + '\t'
                if gene in g2p:
                    str_i += ';'.join(g2p[gene]) + '\t'
                else:
                    str_i += '_' + '\t'
                list_w.write(str_i.strip() + '\n')
        self.option('outlist').set_path(self.output_dir + '/g2t2p.list')
        self.end()

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
            "id": "Extract_gff_" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "protein_transcript.extract_gff",
            "instant": False,
            "options": dict(
                gff_gtf = test_dir + "/" + "GCF_001605985.1_ASM160598v1_genomic.gff",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
