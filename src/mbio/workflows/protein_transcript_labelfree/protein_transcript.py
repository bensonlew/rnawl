# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2018-11-04"

from biocluster.workflow import Workflow
import re
from collections import defaultdict
import unittest
from biocluster.wpm.client import *
import datetime
import copy
import os


class ProteinTranscriptWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinTranscriptWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='genome_info', type='infile', format="labelfree.common"),
            dict(name='pep', type='infile', format="labelfree.common"),
            dict(name="transcript_list", type='infile', format="labelfree.common"),
            dict(name="protein_faa", type='infile', format="labelfree.common"),
            # dict(name="transcript_faa", type='infile', format="labelfree.common"),
            # dict(name="p2t_list", type='outfile', format="labelfree.common"),
            dict(name="rna_type", type='string', default="ref_rna_v2"),
            dict(name="main_id", type='string'),
            dict(name="rna_matrix", type='string', default=""),
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # ncbi的会特殊处理一下
        self.ncbi = False
        with open(self.option('pep').prop['path']) as gr:
            first = gr.readline()
            if re.match('^>rna\d+ ',first) or re.match('^>gene\d+ ',first):
                self.ncbi = True
        if self.option('rna_type') == 'prok_rna':
            self.extract_gff = self.add_tool("protein_transcript_labelfree.extract_ptt")
        elif self.option('rna_type') == 'denovo_rna_v2':
            self.extract_gff = self.add_tool("protein_transcript_labelfree.extract_denovo")
        else:
            if not self.ncbi:
                self.extract_gff = self.add_tool("protein_transcript_labelfree.extract_relation")
            else:
                self.extract_gff = self.add_tool("protein_transcript_labelfree.extract_ncbi")
            # self.extract_gff = self.add_tool("protein_transcript_labelfree.extract_relation_mulity")
        with open(self.option('protein_faa').prop['path']) as pro, open(self.option('genome_info').prop['path']) as trans:
            self.protein_list = [x.split('\n')[0].strip().lstrip('>').split(' ')[0] for x in pro.read().split('\n>')]
            if self.option('rna_type') == 'prok_rna':
                self.transcript_list = list(set([x.strip().split('\t')[6].strip() for x in trans.readlines()]))
            elif self.option('rna_type') == 'denovo_rna_v2':
                self.transcript_list = list(set([x.strip().split('\t')[1].strip() for x in trans.readlines()]))
            elif self.option('rna_type') == 'whole_transcriptome':
                self.filter_exp_matrix(self.option('rna_matrix'), cutoff = "0")
                with open(self.option('rna_matrix'), 'r') as trans2:
                    _ = trans2.readline()
                    self.transcript_list = [x.split('\t')[0].strip() for x in trans2.readlines()]
            else:
                self.transcript_list = list(set([x.strip().split('\t')[0].strip() for x in trans.readlines()]))
        self.copy_protein = copy.copy(self.protein_list)
        print(len(self.protein_list))
        print(len(self.transcript_list))
        self.remove_trans = list()
        self.relation_list = list()

    def run(self):

        self.extract_gff.on("end", self.try_gff)
        self.run_extract_gff()
        super(ProteinTranscriptWorkflow, self).run()

    def set_db(self):
        def get_blast_dict(file):
            q2r = defaultdict(set)
            with open(file, 'r') as blast:
                for line in blast:
                    line = line.strip().split('\t')
                    q2r[line[5]].add(line[10])
            return q2r

        with open(self.output_dir + '/relationship_list', 'w') as rel_w:
            rel_w.write('related' + '\t' + ';'.join(self.relation_list) + '\n')
            rel_w.write('failed_proteins' + '\t' + ';'.join(self.protein_list) + '\n')
            rel_w.write('failed_transcripts' + '\t' + ';'.join(self.transcript_list) + '\n')
        relation = self.api.api("protein_transcript_labelfree.protein_transcript")
        params = dict(
            software = 'blast',
            rna_type = self.option('rna_type'),
        )
        print(len(self.protein_list))
        print(len(self.transcript_list))
        relation.add_relation(main_id = self.option('main_id'), params = params, relation_file= self.output_dir + '/relationship_list')
        relation.update_task_relation(main_id=self.option('main_id'))
        """
        保存结果表到mongo数据库中
        """
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "联合分析蛋白组转录组关联结果"],
            ["relationship_list", "txt", "联合分析蛋白组转录组关联结果"]
        ])
        super(ProteinTranscriptWorkflow, self).end()

    def filter_exp_matrix(self, matrix, cutoff):
        filter_list = list()
        with open(matrix, 'r') as mat_r:
            header = mat_r.readline()
            for line in mat_r:
                tmp = line.strip().split('\t')
                cut = 0
                for i in tmp[1:]:
                    if float(i) < float(cutoff):
                        cut += 1
                if not cut > float(len(tmp))/2:
                    filter_list.append(line)
        with open(matrix, 'w') as mat_w:
            mat_w.write(header)
            for line in filter_list:
                mat_w.write(line)

    def run_extract_gff(self):
        if self.option('rna_type') == 'prok_rna':
            options = dict(
                ptt=self.option('genome_info'),
            )
        elif self.option('rna_type') == 'denovo_rna_v2':
            options = dict(
                pep=self.option('pep'),
                protein_faa=self.option('protein_faa'),
            )
        else:
            if not self.ncbi:
                options = dict(
                    pep = self.option('pep'),
                    relation_file = self.option('genome_info')
                )
            else:
                options = dict(
                    pep=self.option('pep'),
                    protein_faa=self.option('protein_faa'),
                    relation_file=self.option('genome_info')
                )
        self.extract_gff.set_options(options)
        self.extract_gff.run()

    def try_gff(self):
        gff_g2p = defaultdict(set)
        self.gff_p2g = dict()
        with open(self.extract_gff.option('outlist').prop['path'], 'r') as gff_info:
            for line in gff_info:
                line = line.strip().split('\t')
                try:
                    gff_g2p[line[0]] = set(line[1].split(';'))
                    for i in set(line[1].split(';')):
                        self.gff_p2g[i] = line[0]
                except:
                    gff_g2p[line[0]] = set(line[2].split(';'))
                    for i in set(line[2].split(';')):
                        self.gff_p2g[i] = line[0]
        for pro in self.copy_protein:
            for gene in gff_g2p:
                if pro in gff_g2p[gene]:
                    if gene in self.transcript_list or gene in self.remove_trans:
                        self.relation_list.append('%s|%s' % (pro, gene))
                        self.protein_list.remove(pro)
                        try:
                            self.transcript_list.remove(gene)
                            self.remove_trans.append(gene)
                        except:
                            self.remove_trans.append(self.transcript_list.pop())
                        break
                    # else:
                    #     self.relation_list.append('%s|%s' % (pro, gene))
                    #     self.protein_list.remove(pro)
                    #     self.remove_trans.append(gene)
        # 偶尔会有ncbi这种项目gff文件有某个基因而biomart文件没有的现象
        if self.protein_list:
            for pro in self.protein_list[::]:
                for gene in gff_g2p:
                    if pro in gff_g2p[gene]:
                        self.relation_list.append('%s|%s' % (pro, gene))
                        self.protein_list.remove(pro)
                        self.remove_trans.append(gene)
                        break
        if not self.protein_list:
            self.set_db()
        else:
            self.logger.info('对应文件不全，需要blast')
            self.logger.info(self.protein_list)
            self.two_side_blast()

    def two_side_blast(self):
        ref = os.path.join(self.work_dir, 'protein_pep.fasta')
        if os.path.exists(ref):
            os.remove(ref)
        # os.link(self.option('pep').prop['path'], ref)
        cmd = 'cp %s %s'%(self.option('pep').prop['path'], ref)
        os.system(cmd)
        self.query_blast = self.add_tool('protein_transcript_labelfree.blast_forrela')
        self.ref_blast = self.add_tool('protein_transcript_labelfree.blast_forrela')
        options_query = dict(
                query=self.option('protein_faa').prop['path'],
                # reference=self.option('pep').prop['path'],
                reference=ref,
            )
        options_ref = dict(
            reference=self.option('protein_faa').prop['path'],
            # query=self.option('pep').prop['path'],
            query=ref,
        )
        self.query_blast.set_options(options_query)
        self.ref_blast.set_options(options_ref)
        self.on_rely([self.query_blast, self.ref_blast], self.blast_out)
        for blast in [self.query_blast, self.ref_blast]:
            blast.run()

    def blast_out(self):
        query2ref = defaultdict(set)
        with open(self.query_blast.option('outtable').prop['path'], 'r') as query:
            _ = query.readline()
            for line in query.readlines():
                line = line.strip().split('\t')
                if len(line) > 4:
                    query = line[5]
                    # if u'.' in query:
                    #     query = re.sub('\.\d+$', '', query)
                    ref = line[10]
                    if u'.' in ref and self.option('rna_type').startswith('ref'):
                        ref = re.sub('\.\d+$', '', ref)
                    query2ref[query].add(ref)
        # for i in query2ref:
        #     tmp = copy.copy(query2ref[i])
        #     # query2ref[i] = set()
        #     for j in tmp:
        #         if j in self.gff_p2g:
        #             query2ref[i].add(self.gff_p2g[j])
        #         else:
        #             query2ref[i].add(i)
        ref2query = defaultdict(set)
        with open(self.ref_blast.option('outtable').prop['path'], 'r') as ref:
            _ = ref.readline()
            for line in ref.readlines():
                line = line.strip().split('\t')
                if len(line) > 4:
                    query = line[10]
                    # if u'.' in query:
                    #     query = re.sub('\.\d+$', '', query)
                    ref = line[5]
                    if u'.' in ref and self.option('rna_type').startswith('ref'):
                        ref = re.sub('\.\d+$', '', ref)
                    ref2query[ref].add(query)
        # tmp_dict = dict()
        # for i in ref2query:
        #     if i in self.gff_p2g:
        #         tmp_dict[self.gff_p2g[i]] = ref2query[i]
        # ref2query.update(tmp_dict)
        self.copy_protein = copy.copy(self.protein_list)
        # self.logger.info(query2ref)
        for pro in self.copy_protein:
            judge = 0
            if pro in query2ref:
                for ref in query2ref[pro]:
                    try:
                        gene = self.gff_p2g[ref]
                    except:
                        self.logger.info(ref)
                        continue
                    if pro in ref2query[ref] and (gene in self.transcript_list or gene in self.remove_trans):
                        judge = 1
                        self.relation_list.append('%s|%s' % (pro, gene))
                        self.protein_list.remove(pro)
                        try:
                            self.transcript_list.remove(gene)
                            self.remove_trans.append(gene)
                        except:
                            self.remove_trans.append(self.transcript_list.pop())
                        break
                    else:
                        # print(gene,pro)
                        if gene in self.transcript_list or gene in self.remove_trans:
                            self.relation_list.append('%s|%s' % (pro, gene))
                            self.protein_list.remove(pro)
                            try:
                                self.transcript_list.remove(gene)
                                self.remove_trans.append(gene)
                            except:
                                self.remove_trans.append(self.transcript_list.pop())
                            break
                        else:
                            print(gene)
                else:
                    print(judge)
                    if not judge:
                        for ref in query2ref[pro]:
                            try:
                                gene = self.gff_p2g[ref]
                            except:
                                self.logger.info(ref)
                                continue
                            if gene in self.transcript_list or gene in self.remove_trans:
                                self.relation_list.append('%s|%s' % (pro, gene))
                                self.protein_list.remove(pro)
                                try:
                                    self.transcript_list.remove(gene)
                                    self.remove_trans.append(gene)
                                except:
                                    self.remove_trans.append(self.transcript_list.pop())
                                break
            else:
                print(pro)
        self.set_db()


class TestFunction(unittest.TestCase):
    def test(self):
        worker = worker_client()
        id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        data = {
            "id": "protein_transcript_labelfree_" + id,
            "type": "workflow",
            "name": "protein_transcript_labelfree.protein_transcript_labelfree",
            "options": dict(
                rna_type="ref_rna",
                gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/biomart/Homo_sapiens.GRCh37.biomart",
                pep="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/cds/Homo_sapiens.GRCh37.pep.fa",
                protein_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/protein.list",
                transcript_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/gene.list",
                protein_faa="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/DB.fasta",
                )
        }

        info = worker.add_task(data)
        print(info)


    # def test(self):
    #     import random
    #     from mbio.workflows.single import SingleWorkflow
    #     from biocluster.wsheet import Sheet
    #     data = {
    #         "id": "protein_transcript_labelfree",
    #         "type": "workflow",
    #         "name": "protein_transcript_labelfree.protein_transcript_labelfree",
    #         "options": dict(
    #             rna_type="ref_rna",
    #             gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/biomart/Homo_sapiens.GRCh37.biomart",
    #             pep="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/cds/Homo_sapiens.GRCh37.pep.fa",
    #             protein_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/protein.list",
    #             transcript_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/gene.list",
    #             protein_faa="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/DB.fasta",
    #             )
    #     }
    #     # data['options']['method'] = 'rsem'
    #     # wsheet = Sheet(data=data)
    #     # wf = SingleWorkflow(wsheet)
    #     # wf.run()
    #     # #
    #     # data['id'] += '1'
    #     # data['options']['method'] = 'salmon'
    #     # wsheet = Sheet(data=data)
    #     # wf = SingleWorkflow(wsheet)
    #     # wf.run()
    #     #
    #     data['id'] += '_fyt'
    #     wsheet = Sheet(data=data)
    #     wf = SingleWorkflow(wsheet)
    #     wf.run()


if __name__ == '__main__':
    unittest.main()
