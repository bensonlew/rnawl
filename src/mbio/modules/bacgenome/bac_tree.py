# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify: 20180320
import os
import re
import shutil
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
class BacTreeModule(Module):
    """
    16s,kgene等进化树
    """

    def __init__(self, work_id):
        super(BacTreeModule, self).__init__(work_id)
        options = [
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "seq_faa", "type": "infile", "format": "sequence.fasta"},  #gene蛋白文件
            {"name": "seq_fa", "type": "infile", "format": "sequence.fasta"},  # gene核酸文件
            {"name": "sample_name", "type": "string"},
            {"name": "analysis", "type": "string", "default": "uncomplete"}  ###流程分析模式complete，uncomplete
        ]
        self.add_option(options)
        self.step.add_steps("sss", "hgene","tree_16s","tree_hgene",'ssu_align')
        self.tree_16s = self.add_tool("graph.phy_tree")
        self.tree_hgene = self.add_tool("graph.phy_tree")
        self.sss = self.add_tool("bacgenome.get_tree")
        self.hgene = self.add_tool("bacgenome.get_hgene_tree")
        self.ssu_align = self.add_tool("bacgenome.ssu_align")
        self.list = [self.tree_16s,self.tree_hgene]

    def set_step(self, event):
        if "start" in event["data"].keys():
            event["data"]["start"].start()
        if "end" in event["data"].keys():
            event["data"]["end"].finish()
        self.step.update()

    def check_options(self):
        if self.option('analysis') in ['complete']:
            if not self.option("rrna_gff").is_set:
                raise OptionError("必须设置rRNA的gff3文件！", code="21400801")
            if not self.option("genome_fa").is_set:
                raise OptionError("必须设置参考基因组文件！", code="21400802")
            if not self.option("seq_faa").is_set:
                raise OptionError("必须设置蛋白序列文件！", code="21400803")
        if self.option('analysis') in ['uncomplete']:
            if not self.option("seq_faa").is_set:
                raise OptionError("必须设置蛋白序列文件！", code="21400804")
        if not self.option("sample_name"):
            raise OptionError("必须设置样品名称！", code="21400805")

    def run_get_16s(self):
        options = {
            "rrna_gff": self.option("rrna_gff"),
            "genome_fa": self.option("genome_fa"),
            "sample_name": self.option("sample_name"),
        }
        self.sss.set_options(options)
        self.sss.on("start", self.set_step, {"start": self.step.sss})
        self.sss.on("end", self.set_step, {"end": self.step.sss})
        self.sss.on('end', self.set_output, 'sss')
        self.sss.run()

    def run_get_gene(self):
        options = {
            "seq_faa": self.option("seq_faa"),
            "gene_gff": self.option("gene_gff"),
            "seq_fa": self.option("seq_fa"),
            "sample_name": self.option("sample_name"),
        }
        self.hgene.set_options(options)
        self.hgene.on("start", self.set_step, {"start": self.step.hgene})
        self.hgene.on("end", self.set_step, {"end": self.step.hgene})
        self.hgene.on('end', self.set_output, 'hgene')
        self.hgene.run()

    def run_tree_16s(self):
        opts = {
            "fasta": self.ssu_align.option('out'),
            "align": "false"
        }
        self.tree_16s.set_options(opts)
        self.tree_16s.on("start", self.set_step, {"start": self.step.tree_16s})
        self.tree_16s.on("end", self.set_step, {"end": self.step.tree_16s})
        self.tree_16s.on('end', self.set_output, 'tree_16s')
        self.tree_16s.run()

    def run_ssu_align(self):
        opts = {
            "16s_fa": self.sss.option('out')
        }
        self.ssu_align.set_options(opts)
        self.ssu_align.on("start", self.set_step, {"start": self.step.ssu_align})
        self.ssu_align.on("end", self.set_step, {"end": self.step.ssu_align})
        self.ssu_align.on('end', self.set_output, 'ssu_align')
        self.ssu_align.run()

    def run_tree_hgene(self):
        opts = {
            "fasta": self.hgene.option('out'),
            "sequence_type":'amino_acid',
        }
        self.tree_hgene.set_options(opts)
        self.tree_hgene.on("start", self.set_step, {"start": self.step.tree_hgene})
        self.tree_hgene.on("end", self.set_step, {"end": self.step.tree_hgene})
        self.tree_hgene.on('end', self.set_output, 'tree_hgene')
        self.tree_hgene.run()

    def set_output(self, event):
        self.logger.info("开始set output")
        obj = event["bind_object"]
        if self.option('analysis') in ['uncomplete']:
            if event['data'] == 'hgene':
                if os.path.exists(self.output_dir + '/hgene'):
                    shutil.rmtree(self.output_dir + '/hgene')
                os.mkdir(self.output_dir + '/hgene')
                if os.path.exists(self.hgene.work_dir + '/all.hgene.fa'):  ## add by qingchen.zhang @20200820
                    os.link(self.hgene.work_dir + '/all.hgene.fa',
                            self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_gene.fa')
                if os.path.exists(self.hgene.work_dir + '/' + self.option('sample_name') +".cor_gene.fa"):
                    os.link(self.hgene.work_dir + '/' + self.option('sample_name') +".cor_gene.fa",
                            self.output_dir + '/hgene/' + self.option('sample_name') + '.corgene.fa')
                if os.path.exists(self.hgene.work_dir + '/' + self.option(
                        'sample_name') + '.cor_gene.fnn'):  ## add by gao.hao @20210621
                    os.link(self.hgene.work_dir + '/' + self.option('sample_name') + '.cor_gene.fnn',
                            self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_gene.fnn')
                if os.path.exists(self.hgene.work_dir + '/' + self.option(
                        'sample_name') + '.cor_gene.faa'):  ## add by gao.hao @20210621
                    os.link(self.hgene.work_dir + '/' + self.option('sample_name') + '.cor_gene.faa',
                            self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_gene.faa')
                if os.path.exists(self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_blast.xls'):  #zouguanqing 2019.328
                    os.remove(self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_blast.xls')
                if os.path.exists(self.hgene.work_dir + '/' + self.option('sample_name') + '.coregene.m8'):## add by qingchen.zhang @20200820
                    os.link(self.hgene.work_dir + '/' + self.option('sample_name') + '.coregene.m8' , self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_blast.xls')
                if os.path.exists(self.output_dir + '/hgene/' + self.option(
                        'sample_name') + '.house_blast.xls'):  # gaohao 2021.05.18
                    os.remove(self.output_dir + '/hgene/' + self.option('sample_name') + '.house_blast.xls')
                if os.path.exists(self.hgene.work_dir + '/' + self.option('sample_name') + '.last.m8'):  # gaohao 2021.05.18
                    os.link(self.hgene.work_dir + '/' + self.option('sample_name') + '.last.m8',
                            self.output_dir + '/hgene/' + self.option('sample_name') + '.house_blast.xls')
            if event['data'] == 'tree_hgene':
                if os.path.exists(self.output_dir + '/hgene/' + self.option('sample_name') + '.phylo_tree.nwk'):
                    os.remove(self.output_dir + '/hgene/' + self.option('sample_name') + '.phylo_tree.nwk')
                os.link(self.tree_hgene.work_dir + '/' + 'phylo_tree.nwk',self.output_dir + '/hgene/' + self.option('sample_name') + '.phylo_tree.nwk')
        elif self.option('analysis') in ['complete']:
            if self.lines >= 1:
                if event['data'] == 'hgene':
                    if os.path.exists(self.output_dir + '/hgene'):
                        shutil.rmtree(self.output_dir + '/hgene')
                    os.mkdir(self.output_dir + '/hgene')
                    if os.path.exists(self.hgene.work_dir + '/all.hgene.fa'):  ## add by qingchen.zhang @20200820
                        os.link(self.hgene.work_dir + '/all.hgene.fa',
                                self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_gene.fa')
                    if os.path.exists(self.hgene.work_dir + '/' + self.option('sample_name') + ".cor_gene.fa"):
                        os.link(self.hgene.work_dir + '/' + self.option('sample_name') + ".cor_gene.fa",
                                self.output_dir + '/hgene/' + self.option('sample_name') + '.corgene.fa')
                    if os.path.exists(self.hgene.work_dir + '/' + self.option(
                            'sample_name') + '.cor_gene.fnn'):  ## add by gao.hao @20210621
                        os.link(self.hgene.work_dir + '/' + self.option('sample_name') + '.cor_gene.fnn',
                                self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_gene.fnn')
                    if os.path.exists(self.hgene.work_dir + '/' + self.option(
                            'sample_name') + '.cor_gene.faa'):  ## add by gao.hao @20210621
                        os.link(self.hgene.work_dir + '/' + self.option('sample_name') + '.cor_gene.faa',
                                self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_gene.faa')
                    if os.path.exists(self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_blast.xls'):
                        os.remove(self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_blast.xls')
                    if os.path.exists(self.hgene.work_dir + '/' + self.option('sample_name') + '.coregene.m8'):## add by qingchen.zhang @20200820
                        os.link(self.hgene.work_dir + '/' + self.option('sample_name') + '.coregene.m8',
                            self.output_dir + '/hgene/' + self.option(
                                'sample_name') + '.cor_blast.xls')  # zouguanqing 2019.328
                    if os.path.exists(self.output_dir + '/hgene/' + self.option(
                            'sample_name') + '.house_blast.xls'):  # gaohao 2021.05.18
                        os.remove(self.output_dir + '/hgene/' + self.option('sample_name') + '.house_blast.xls')
                    if os.path.exists(self.hgene.work_dir + '/'+self.option('sample_name') + '.last.m8'):  # gaohao 2021.05.18
                        os.link(self.hgene.work_dir + '/'+self.option('sample_name') + '.last.m8',
                                self.output_dir + '/hgene/' + self.option('sample_name') + '.house_blast.xls')
                if event['data'] == 'tree_hgene':
                    if os.path.exists(self.output_dir + '/hgene/' + self.option('sample_name') + '.phylo_tree.nwk'):
                        os.remove(self.output_dir + '/hgene/' + self.option('sample_name') + '.phylo_tree.nwk')
                    os.link(self.tree_hgene.work_dir + '/' + 'phylo_tree.nwk',
                            self.output_dir + '/hgene/' + self.option('sample_name') + '.phylo_tree.nwk')
                if event['data'] == 'sss':
                    if not os.path.exists(self.output_dir + '/16s'):
                        os.mkdir(self.output_dir + '/16s')
                    self.logger.info(self.sss.work_dir + '/16s.fasta')
                    self.logger.info(self.output_dir + '/16s/' + self.option('sample_name') + '.16s.fa')
                    if os.path.exists(self.output_dir + '/16s/' + self.option('sample_name') + '.16s.fa'):
                        os.remove(self.output_dir + '/16s/' + self.option('sample_name') + '.16s.fa')
                    os.link(self.sss.work_dir + '/16s.fasta',
                            self.output_dir + '/16s/' + self.option('sample_name') + '.16s.fa')
                    if os.path.exists(self.output_dir + '/16s/' + self.option(
                            'sample_name') + '.16s_blast.xls'):  # gaohao 2021.05.18
                        os.remove(self.output_dir + '/16s/' + self.option('sample_name') + '.16s_blast.xls')
                    if os.path.exists(self.sss.work_dir + '/blast_out.xls'):  # gaohao 2021.05.18
                        os.link(self.sss.work_dir + '/blast_out.xls',
                                self.output_dir + '/16s/' + self.option('sample_name') + '.16s_blast.xls')
                if event['data'] == 'tree_16s':
                    if os.path.exists(self.output_dir + '/16s/' + self.option('sample_name') + '.phylo_tree.nwk'):
                        os.remove(self.output_dir + '/16s/' + self.option('sample_name') + '.phylo_tree.nwk')
                    os.link(self.tree_16s.work_dir + '/' + 'phylo_tree.nwk',
                            self.output_dir + '/16s/' + self.option('sample_name') + '.phylo_tree.nwk')
            else:
                if event['data'] == 'hgene':
                    if os.path.exists(self.output_dir + '/hgene'):
                        shutil.rmtree(self.output_dir + '/hgene')
                    os.mkdir(self.output_dir + '/hgene')
                    if os.path.exists(self.hgene.work_dir + '/all.hgene.fa'):## add by qingchen.zhang @20200820
                        os.link(self.hgene.work_dir + '/all.hgene.fa',
                            self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_gene.fa')
                    if os.path.exists(self.hgene.work_dir + '/' + self.option('sample_name') + ".cor_gene.fa"):
                        os.link(self.hgene.work_dir + '/' + self.option('sample_name') + ".cor_gene.fa",
                                self.output_dir + '/hgene/' + self.option('sample_name') + '.corgene.fa')
                    if os.path.exists(self.hgene.work_dir + '/'+ self.option('sample_name') + '.cor_gene.fnn'):## add by gao.hao @20210621
                        os.link(self.hgene.work_dir + '/'+ self.option('sample_name') + '.cor_gene.fnn',
                            self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_gene.fnn')
                    if os.path.exists(self.hgene.work_dir + '/'+ self.option('sample_name') + '.cor_gene.faa'):## add by gao.hao @20210621
                        os.link(self.hgene.work_dir + '/'+ self.option('sample_name') + '.cor_gene.faa',
                            self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_gene.faa')
                    if os.path.exists(self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_blast.xls'):
                        os.remove(self.output_dir + '/hgene/' + self.option('sample_name') + '.cor_blast.xls')
                    if os.path.exists(self.hgene.work_dir + '/' + self.option('sample_name') + '.coregene.m8'): ## add by qingchen.zhang @20200820
                        os.link(self.hgene.work_dir + '/' + self.option('sample_name') + '.coregene.m8',
                            self.output_dir + '/hgene/' + self.option(
                                'sample_name') + '.cor_blast.xls')  # zouguanqing 2019.328
                    if os.path.exists(self.output_dir + '/hgene/' + self.option(
                            'sample_name') + '.house_blast.xls'):  # gaohao 2021.05.18
                        os.remove(self.output_dir + '/hgene/' + self.option('sample_name') + '.house_blast.xls')
                    if os.path.exists(
                            self.hgene.work_dir + '/' + self.option('sample_name') + '.last.m8'):  # gaohao 2021.05.18
                        os.link(self.hgene.work_dir + '/' + self.option('sample_name') + '.last.m8',
                                self.output_dir + '/hgene/' + self.option('sample_name') + '.house_blast.xls')
                if event['data'] == 'tree_hgene':
                    if os.path.exists(self.output_dir + '/hgene/' + self.option('sample_name') + '.phylo_tree.nwk'):
                        os.remove(self.output_dir + '/hgene/' + self.option('sample_name') + '.phylo_tree.nwk')
                    os.link(self.tree_hgene.work_dir + '/' + 'phylo_tree.nwk',
                            self.output_dir + '/hgene/' + self.option('sample_name') + '.phylo_tree.nwk')

    def run(self):
        super(BacTreeModule, self).run()
        if self.option('analysis') in ['uncomplete']:
            self.tree_hgene.on('end',self.end)
            self.hgene.on('end',self.run_tree_hgene)
            self.run_get_gene()
        elif self.option('analysis') in ['complete']:
            self.lines = self.get_line(self.option('rrna_gff').prop['path'])
            if self.lines >= 1:
                self.on_rely(self.list, self.end)
                self.hgene.on('end', self.run_tree_hgene)
                self.ssu_align.on("end", self.run_tree_16s)
                self.sss.on('end',self.run_ssu_align)
                self.run_get_16s()
                self.run_get_gene()
            else:
                self.tree_hgene.on('end',self.end)
                self.hgene.on('end',self.run_tree_hgene)
                self.run_get_gene()

    def get_line(self,file):
        list1 = []
        with open (file,'r') as f:
            lines =f.readlines()
            for line in lines[1:]:
                if re.search(r"16S_rRNA", line.strip()):
                    list1.append(line)
        return len(list1)
        
    def end(self):
        super(BacTreeModule, self).end()