# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import unittest
import re
from collections import defaultdict


class ExtractDenovoAgent(Agent):
    def __init__(self, parent):
        super(ExtractDenovoAgent, self).__init__(parent)
        options = [
            {"name": "pep", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
            {"name": "protein_faa", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
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
        return True

    def set_resource(self):
        self._cpu = 15
        self._memory = '40G'

    def end(self):
        super(ExtractDenovoAgent, self).end()


class ExtractDenovoTool(Tool):
    def __init__(self, config):
        super(ExtractDenovoTool, self).__init__(config)

    def extract_pep(self):
        # 由于线下会把id给处理一下，没办法，只能匹配序列
        with open(self.option('protein_faa').prop['path'], 'r') as fa_r:
            faa_info = fa_r.read().split('\n>')
            faa_tuples = [(fa.split('\n')[0].strip().lstrip('>'), ''.join(fa.strip().split('\n')[1:])) for fa in faa_info]
        g2p = defaultdict(set)
        with open(self.option('pep').prop['path'], 'r') as pep:
            for block in pep.read().split('\n>'):
                pro = block.split('\n')[0].strip().lstrip(">").split(" ")[0]
                gene = pro.split('::')[0]
                if '_orf' in gene:
                    if 'i' in gene.split('_')[-2]:
                        gene = '_'.join(gene.split('_')[0:-2])
                    else:
                        gene = '_'.join(gene.split('_')[0:-1])
                    g2p[gene].add(pro)
                # 有的时候转录拼接成的fasta序列文件竟然都有断行之类的，会导致根据序列匹配不到什么东西，用id try一下吧还是
                else:
                    try:
                        converted = pro.split('::')
                        converted = '%s_%s' % (converted[1], converted[-1])
                        g2p[gene].add(converted)
                    except:
                        pass
                seq = ''.join(block.strip().split('\n')[1:])
                for id, s in faa_tuples:
                    if s == seq or seq.endswith(s):
                        if gene in g2p:
                            if gene in list(g2p[gene])[0]:
                                continue
                            else:
                                g2p[gene] = set()
                        g2p[gene].add(id)
                        g2p[gene].add(pro)
                        break
        return g2p

    def run(self):
        super(ExtractDenovoTool, self).run()
        g2p = self.extract_pep()
        with open(self.output_dir + '/g2t2p.list', 'w') as list_w:
            genes = sorted(g2p.keys())
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
            "name": "protein_transcript.extract_denovo",
            "instant": False,
            "options": dict(
                pep = test_dir + "/" + "Trinity.filter.fasta.transdecoder.pep",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
