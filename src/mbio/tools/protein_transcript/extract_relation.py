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


class ExtractRelationAgent(Agent):
    def __init__(self, parent):
        super(ExtractRelationAgent, self).__init__(parent)
        options = [
            {"name": "pep", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
            {"name": "relation_file", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
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
        super(ExtractRelationAgent, self).end()


class ExtractRelationTool(Tool):
    def __init__(self, config):
        super(ExtractRelationTool, self).__init__(config)
        self.g2p = dict()
        with open(self.option('pep').prop['path'], 'r') as pep_r:
            pep = pep_r.read().split('\n>')
            self.protein_list = [x.lstrip('>').split('\n')[0].strip().split(' ')[0].strip() for x in pep]
            if '\t' in self.protein_list[0]:
                self.protein_list = [pro.split('\t')[0] for pro in self.protein_list]
            for n in range(len(self.protein_list)):
                if u'.' in self.protein_list[n] and 'EN' in self.protein_list[n]:
                    self.protein_list[n] = re.sub('\.\d+$', '', self.protein_list[n])
            #     我真是服了，有的对应关系biomart文件里没有，竟然在pep文件里有。。。
            pg_re = re.compile('(ENS\w*P\d+\.\d+).*(ENS\w*G\d+).*?')
            for p in pep:
                p = p.split('\n')[0]
                self.logger.info(p)
                pg = pg_re.search(p)
                if not pg:
                    self.logger.info('not match')
                    continue
                pro, gene = pg.groups(1)
                if gene not in self.g2p:
                    self.g2p[gene] = set()
                self.g2p[gene].add(pro)
                self.g2p[gene].add(re.sub('\.\d+$', '', pro))

    def extract_relation(self, line):
        tmp_list = copy.copy(self.protein_list)
        for pro in tmp_list:
            if re.search('\t%s\t'%pro, line):
                self.g2p[line.split('\t')[0].strip()].add(pro)
                # self.protein_list.remove(pro)
                break

    def deal_biomart(self):
        self.logger.info('已经在处理biomart文件')
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
                # self.set_error("无法从biomart文件中得到基因蛋白对应关系%s"%self.option('relation_file').prop['path'], variables = (self.option('relation_file').prop['path']),code="35004101")
                self.logger.info("无法从biomart文件中得到基因蛋白对应关系%s"%self.option('relation_file').prop['path'])
                tg = self.search_tran2gene()
                if not tg:
                    self.logger.info("无法找到tran2gene，只能blast")
                else:
                    with open(tg) as tgw:
                        for line in tgw:
                            line = line.strip().split('\t')
                            if not line:
                                continue
                            if len(line) >= 5:
                                gene = line[1]
                                protein = line[4]
                                if gene not in self.g2p:
                                    self.g2p[gene] = set()
                                if protein in self.protein_list:
                                    self.g2p[gene].add(protein)

            else:
                rel.seek(0,0)
                for line in rel:
                    if line.strip():
                        line = line.strip().split('\t')
                        if line[0] not in self.g2p:
                            self.g2p[line[0]] = set()
                        for n in p_cols:
                            if line[n] in self.protein_list:
                                self.g2p[line[0]].add(line[n])

    def search_tran2gene(self):
        # biomart_path = os.path.abspath(self.option('relation_file').prop['path'])
        biomart_path = self.option('relation_file').prop['path']
        self.logger.info(biomart_path)
        gene_path = '/'.join(biomart_path.split('/')[:-2])
        self.logger.info(gene_path)
        for root, dirs, files in os.walk(gene_path):
            for file in files:
                if file.endswith('tran2gene.txt'):
                    self.logger.info('tran2gene.txt')
                    return os.path.join(root, file)
        return

    def run(self):
        super(ExtractRelationTool, self).run()
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
