# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import unittest
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import copy
import textwrap
import random


class ExtractNcbiAgent(Agent):
    def __init__(self, parent):
        super(ExtractNcbiAgent, self).__init__(parent)
        options = [
            {"name": "pep", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
            {"name": "relation_file", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
            {"name": "protein_faa", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
            {"name": "outlist", "type": "outfile", "format": "itraq_and_tmt.common"},  # 输出的有关联关系的list文件
            ]
        self.add_option(options)
        self.step.add_steps('relation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.relation.start()
        self.step.update()

    def step_end(self):
        self.step.relation.finish()
        self.step.update()

    def check_options(self):
        if not os.path.exists(self.option('relation_file').prop['path']):
            raise OptionError("没有正确传入biomart文件")
        if not os.path.exists(self.option('pep').prop['path']):
            raise OptionError("没有正确传入蛋白序列文件")
        return True

    def set_resource(self):
        self._cpu = 15
        self._memory = '40G'

    def end(self):
        super(ExtractNcbiAgent, self).end()


class ExtractNcbiTool(Tool):
    def __init__(self, config):
        super(ExtractNcbiTool, self).__init__(config)
        with open(self.option('pep').prop['path'], 'r') as pep_r:
            pep = pep_r.read().split('\n>')
            self.protein_list = [x.lstrip('>').split('\n')[0].strip().split(' ')[0].strip() for x in pep]
        self.gff_path = self.search_gff()

    def search_gff(self):
        # biomart_path = os.path.abspath(self.option('relation_file').prop['path'])
        biomart_path = self.option('relation_file').prop['path']
        self.logger.info(biomart_path)
        gene_path = '/'.join(biomart_path.split('/')[:-2])
        self.logger.info(gene_path)
        for root, dirs, files in os.walk(gene_path):
            for file in files:
                if file.endswith('gff'):
                    self.logger.info('找到了gff文件')
                    return os.path.join(root, file)
        return

    def extract_gff(self):
        rp2p = defaultdict(set)
        with open(self.gff_path) as gff_r:
            for line in gff_r:
                if not line.strip() or line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                if len(line) < 8:
                    continue
                if line[2] == 'gene':
                    for att in line[8].split(';'):
                        if 'ID=' in att:
                            gene = att.split('ID=')[1]
                            break
                if line[2] == 'CDS':
                    for att in line[8].split(';'):
                        if 'Parent=' in att:
                            parent = att.split('Parent=')[1]
                        if 'protein_id=' in att:
                            protein = att.split('protein_id=')[1]
                            rp2p[parent].add(protein)
                            rp2p[gene].add(protein)
        return rp2p

    def extract_pep(self):
        def judge_in_or_not(a, b):
            a = textwrap.wrap(a, width=10)
            e = 0
            for n in range(int(len(a) / 3)):
                tmp = random.choice(a)
                a.remove(tmp)
                if tmp not in b:
                    # return False
                    e += 1
                    if e > 1:
                        return False
            return True

        # 由于线下会把id给处理一下，没办法，只能匹配序列
        with open(self.option('protein_faa').prop['path'], 'r') as fa_r:
            faa_info = fa_r.read().split('\n>')
            faa_tuples = [(fa.split('\n')[0].strip().lstrip('>'), ''.join(fa.strip().split('\n')[1:])) for fa in faa_info]
        rp2p = defaultdict(set)
        with open(self.option('pep').prop['path'], 'r') as pep:
            for block in pep.read().split('\n>'):
                pro = block.split('\n')[0].strip().lstrip(">").split(" ")[0]
                seq = ''.join(block.strip().split('\n')[1:]).strip('.')
                # seq = ''.join(block.strip().split('\n')[1:]).replace('.', '')
                tmp = False
                for id, s in faa_tuples:
                    if s in seq or seq in s:
                        minus = len(s) - len(seq)
                        if minus > 10 or minus < -10:
                            continue
                    # if judge_in_or_not(s, seq):
                        rp2p[pro].add(id)
                        self.logger.info('%s_%s'%(pro,id))
                        tmp = True
                        break
                else:
                    if r'.' in seq or tmp:
                        continue
                    self.logger.info(pro)
                    for id, s in faa_tuples:
                        minus = len(s) - len(seq)
                        if minus > 10 or minus < -10:
                            continue
                        if judge_in_or_not(s, seq):
                            rp2p[pro].add(id)
                            self.logger.info('%s_%s' % (pro, id))
                            break
        return rp2p

    def extract_relation(self, line):
        tmp_list = copy.copy(self.protein_list)
        for pro in tmp_list:
            if re.search('\t%s\t'%pro, line):
                self.g2p[line.split('\t')[0].strip()].add(pro)
                # self.protein_list.remove(pro)
                break

    def deal_biomart(self):
        if self.gff_path:
            rp2p = self.extract_gff()
        else:
            rp2p = self.extract_pep()
        self.logger.info(rp2p)
        self.g2p = dict()
        p_cols = list()
        with open(self.option('relation_file').prop['path'], 'r') as rel:
            for line in rel:
                if line.strip():
                    line = line.strip().split('\t')
                    for n,i in enumerate(line):
                        if i in self.protein_list:
                            p_cols.append(n)
                    if p_cols:
                        break
            if not p_cols:
                self.set_error("无法从biomart文件中得到基因蛋白对应关系%s"%self.option('relation_file').prop['path'], variables = (self.option('relation_file').prop['path']),code="35004101")
            rel.seek(0,0)
            for line in rel:
                if line.strip():
                    line = line.strip().split('\t')
                    if line[0] not in self.g2p:
                        self.g2p[line[0]] = set()
                    for n in p_cols:
                        if line[n] in self.protein_list:
                            self.g2p[line[0]].add(line[n])
                            if line[n] in rp2p:
                                for i in list(rp2p[line[n]]):
                                    self.g2p[line[0]].add(i)
        for rp in rp2p:
            if rp not in self.g2p:
                self.g2p[rp] = rp2p[rp]

    def run(self):
        super(ExtractNcbiTool, self).run()
        # self.g2p = defaultdict(set)
        # with open(self.option('relation_file').prop['path'], 'r') as rel:
        # #     rel_info = rel.readlines()
        # # with ThreadPoolExecutor(6) as pool:
        # #     pool.map(self.extract_relation, rel_info)
        #     for line in rel:
        #         self.extract_relation(line)
        self.deal_biomart()
        with open(self.output_dir + '/g2p.list', 'w') as list_w:
            genes = sorted(list(set(self.g2p.keys())))
            for gene in genes:
                if not self.g2p[gene]:
                    continue
                str_i = gene + '\t'
                if gene in self.g2p:
                    str_i += ';'.join(self.g2p[gene]) + '\t'
                else:
                    str_i += '_' + '\t'
                list_w.write(str_i.strip('\t') + '\n')
        self.option('outlist').set_path(self.output_dir + '/g2p.list')
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87'
        data = {
            "id": "Extract_relation_" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "protein_transcript.extract_relation",
            "instant": False,
            "options": dict(
                pep = test_dir + "/" + "cds/Homo_sapiens.GRCh37.pep.fa",
                relation_file = test_dir + "/" + "biomart/Homo_sapiens.GRCh37.biomart",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
