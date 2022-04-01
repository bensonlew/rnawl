# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2018-11-04"
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
import re
from collections import defaultdict
import unittest
from biocluster.wpm.client import *
import datetime
import copy
import os
from biocluster.config import Config
import json
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.file import download



class ToolProteinTranscriptWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ToolProteinTranscriptWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'rna_type', 'type': 'string'},
            {'name': 'rna_task_id', 'type': 'string'},
            {'name': 'protein_type', 'type': 'string'},
            {'name': 'protein_task_id', 'type': 'string'},
            dict(name='genome_info', type='infile', format="itraq_and_tmt.common"),
            dict(name='pep', type='infile', format="itraq_and_tmt.common"),
            dict(name="transcript_list", type='infile', format="itraq_and_tmt.common"),
            dict(name="protein_faa", type='infile', format="itraq_and_tmt.common"),
            dict(name="main_id", type='string'),
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        # 获取原本在接口中获取的文件
        # genome_info, pep, transcript_list, protein_faa
        self.raw_workflow_db_version = self.config.DBVersion
        self.logger.info("工作流db_version={}\n".format(self.config.DBVersion))


        if self.option('rna_type') == 'denovo_rna_v2':
            task_info, db_version = get_task_info_new('denovo_rna_v2', self.option('rna_task_id'))
            try:
                self.genome_info_s3 = self.gene = task_info['assemble_t2g']
                self.genome_info = download_from_s3(self.genome_info_s3, inter_dir=self.work_dir)
            except:
                raise OptionError("没能找到此无参转录组组装信息文件")
            try:
                self.pep_s3 = task_info['bedpath'].replace('.bed', '.pep')
                self.pep = download_from_s3(self.pep_s3, inter_dir=self.work_dir)
            except:
                raise OptionError("没能找到此无参转录组的预测蛋白文件")
            try:
                # self.gene = task_info['unigene_fa'].replace('.fasta', '.unigene.fasta')
                self.gene_s3 = task_info['unigene_fa']
                self.gene = download_from_s3(self.gene_s3, inter_dir=self.work_dir)
            except:
                raise OptionError("没能找到此无参转录组的组装后得到的unigene的fasta文件")

        if self.option('rna_type') == 'ref_rna_v2':
            task_info, db_version = get_task_info_new('ref_rna_v2', self.option('rna_task_id'))
            try:
                # biomart, type = self.ref_rna_v2.get_des_type(self.rna_task_id)
                biomart, type = get_des_type(task_info, self.option('rna_task_id'), db_version)
                self.genome_info = self.gene = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish/',
                                                        biomart)
                if not os.path.exists(self.genome_info):
                    raise OptionError("没能找到此有参转录组的物种biomart文件")
            except:
                raise OptionError("没能找到此有参转录组的物种信息")
            try:
                pep = get_pep(task_info, self.option('rna_task_id'), db_version)
                self.pep = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish/', pep)
                if not os.path.exists(self.genome_info):
                    raise OptionError("没能找到此有参转录组的物种pep文件")
            except:
                raise OptionError("没能找到此有参转录组的物种信息")

        if self.option('rna_type') == 'whole_transcriptome':
            task_info, db_version = get_task_info_new_whole_transcript('whole_transcriptome', self.option('rna_task_id'))
            print 'task_info'
            print task_info
            try:
                # biomart, type = self.ref_rna_v2.get_des_type(self.rna_task_id)
                biomart, type = get_des_type(task_info, self.option('rna_task_id'), db_version)
                print 'biomart'
                print biomart
                self.genome_info = self.gene = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish/',
                                                        biomart)
                if not os.path.exists(self.genome_info):
                    raise OptionError("没能找到此有参转录组的物种biomart文件")
            except:
                raise OptionError("没能找到此有参转录组的物种信息")
            try:
                pep = get_pep(task_info, self.option('rna_task_id'), db_version)
                self.pep = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish/', pep)
                if not os.path.exists(self.genome_info):
                    raise OptionError("没能找到此有参转录组的物种pep文件")
            except:
                raise OptionError("没能找到此有参转录组的物种信息")



        if self.option('rna_type') == 'prok_rna':
            task_info, db_version = get_task_info_new('prok_rna', self.option('rna_task_id'))
            if 'rock_index' in task_info:

                self.genome_info_s3 = self.gene_s3 = os.path.join(task_info['rock_index'], 'ptt.bed')
                self.genome_info = download_from_s3(self.genome_info_s3, inter_dir=self.work_dir)
                self.gene = download_from_s3(self.gene_s3, inter_dir=self.work_dir)
                self.pep_s3 = os.path.join(task_info['rock_index'], 'cds.faa')
                self.pep = download_from_s3(self.pep_s3, inter_dir=self.work_dir)
            else:
                raise OptionError("没能找到此原核转录组的基因组信息")
        protein_task, db_version_protein = get_task_info_new(self.option('protein_type'), self.option('protein_task_id'))
        self.protein_fa_s3 = protein_task['protein_fa']
        self.protein_fa = download_from_s3(self.protein_fa_s3, inter_dir=self.work_dir)

        # ncbi的会特殊处理一下
        self.ncbi = False
        with open(self.pep) as gr:
            first = gr.readline()
            if re.match('^>rna\d+ ',first) or re.match('^>gene\d+ ',first):
                self.ncbi = True
        if self.option('rna_type') == 'prok_rna':
            self.extract_gff = self.add_tool("protein_transcript.extract_ptt")
        elif self.option('rna_type') == 'denovo_rna_v2':
            self.extract_gff = self.add_tool("protein_transcript.extract_denovo")
        else:
            if not self.ncbi:
                self.extract_gff = self.add_tool("protein_transcript.extract_relation")
            else:
                self.extract_gff = self.add_tool("protein_transcript.extract_ncbi")
            # self.extract_gff = self.add_tool("protein_transcript.extract_relation_mulity")
        with open(self.protein_fa) as pro, open(self.genome_info) as trans:
            self.protein_list = [x.split('\n')[0].strip().lstrip('>').split(' ')[0] for x in pro.read().split('\n>')]
            if self.option('rna_type') == 'prok_rna':
                self.transcript_list = list(set([x.strip().split('\t')[6].strip() for x in trans.readlines()]))
            elif self.option('rna_type') == 'denovo_rna_v2':
                self.transcript_list = list(set([x.strip().split('\t')[1].strip() for x in trans.readlines()]))
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
        super(ToolProteinTranscriptWorkflow, self).run()

    def set_db(self):
        def get_blast_dict(file):
            q2r = defaultdict(set)
            with open(file, 'r') as blast:
                for line in blast:
                    line = line.strip().split('\t')
                    q2r[line[5]].add(line[10])
            return q2r

        self.config.DBVersion = self.raw_workflow_db_version
        with open(self.work_dir + '/relationship_list', 'w') as rel_w:
            rel_w.write('related' + '\t' + ','.join(self.relation_list) + '\n')
            rel_w.write('failed_proteins' + '\t' + ','.join(self.protein_list) + '\n')
            rel_w.write('failed_transcripts' + '\t' + ','.join(self.transcript_list) + '\n')
        relation = self.api.api("tool_lab.protein_transcript")
        params = dict(
            software = 'blast',
            rna_type = self.option('rna_type'),
        )

        relate_file = os.path.join(self.output_dir, 'gene_protein_relation.txt')
        only_gene = os.path.join(self.output_dir, 'only_gene.txt')
        only_protein = os.path.join(self.output_dir, 'only_protein.txt')
        with open(relate_file, 'w') as r, open(only_gene, 'w') as og, open(only_protein, 'w') as op:
            for i in self.relation_list:
                pro, gene = i.strip().split('|')
                r.write(pro + '\t' + gene + '\n')
            for g in self.transcript_list:
                og.write(g + '\n')
            for p in self.protein_list:
                op.write(p + '\n')

        print(len(self.protein_list))
        print(len(self.transcript_list))
        relation.add_relation(main_id = self.option('main_id'), params = params, relation_file=self.work_dir + '/relationship_list')
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
        super(ToolProteinTranscriptWorkflow, self).end()

    def run_extract_gff(self):
        if self.option('rna_type') == 'prok_rna':
            options = dict(
                ptt=self.genome_info,
            )
        elif self.option('rna_type') == 'denovo_rna_v2':
            options = dict(
                pep=self.pep,
                protein_faa=self.protein_fa,
            )
        else:
            if not self.ncbi:
                options = dict(
                    pep=self.pep,
                    relation_file=self.genome_info
                )
            else:
                options = dict(
                    pep=self.pep,
                    protein_faa=self.protein_fa,
                    relation_file=self.genome_info
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
        cmd = 'cp %s %s'%(self.pep, ref)
        os.system(cmd)
        self.query_blast = self.add_tool('protein_transcript.blast_forrela')
        self.ref_blast = self.add_tool('protein_transcript.blast_forrela')
        options_query = dict(
                query=self.protein_fa,
                # reference=self.option('pep').prop['path'],
                reference=ref,
            )
        options_ref = dict(
            reference=self.protein_fa,
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


def download_from_s3(from_file, to_path="download/", inter_dir="", cover=True):
    base_name = os.path.basename(from_file)
    to_file = os.path.join(inter_dir, base_name)
    print('from {} to {}'.format(from_file, to_file))
    if os.path.exists(to_file) and os.path.getsize(to_file) != 0:
        pass
    else:
        ## 接口不允许多线程下载，20181217
        # transfer = MultiFileTransfer()
        # print('from {} to {}'.format(from_file, to_file))
        # transfer.add_download(os.path.dirname(from_file) + "/" , os.path.dirname(to_file) + "/")
        # transfer.perform()
        download(from_file, to_file)
    return to_file


def get_task_info_new(project_type, task_id):
    try:
        Config().DBVersion = 0
        db = Config().get_mongo_client(mtype=project_type, db_version=0)[Config().get_mongo_dbname(project_type, db_version=0)]
        task_info = db['sg_task'].find_one({"task_id": task_id})
        task_info['_id']
        return task_info, 0
    except:
        Config().DBVersion = 1
        db = Config().get_mongo_client(mtype=project_type, db_version=1)[Config().get_mongo_dbname(project_type, db_version=1)]
        task_info = db['sg_task'].find_one({"task_id": task_id})
        task_info['_id']
        return task_info, 1

def get_task_info_new_whole_transcript(project_type, task_id):
    try:
        Config().DBVersion = 0
        db = Config().get_mongo_client(mtype=project_type, db_version=0)[Config().get_mongo_dbname(project_type, db_version=0)]
        task_info = db['task'].find_one({"task_id": task_id})
        task_info['_id']
        task_info_long_task_id = task_info['options']['long_task_id']
        task_info = db['task'].find_one({"task_id": task_info_long_task_id})['options']

        return task_info, 0
    except:
        Config().DBVersion = 1
        db = Config().get_mongo_client(mtype=project_type, db_version=1)[Config().get_mongo_dbname(project_type, db_version=1)]
        task_info = db['task'].find_one({"task_id": task_id})
        task_info['_id']
        task_info_long_task_id = task_info['options']['long_task_id']
        task_info = db['task'].find_one({"task_id": task_info_long_task_id})['options']
        return task_info, 1


def get_des_type(result_task, task_id, db_version, table_name="sg_annotation_stat"):
    if "genome_id" in result_task.keys():
        genome_id = result_task["genome_id"]
        db = Config().get_mongo_client(mtype='ref_rna_v2', db_version=db_version)[
            Config().get_mongo_dbname('ref_rna_v2', db_version=db_version)]
        col = db["sg_genome_db"]
        genome_info = col.find_one({"genome_id": genome_id})
        db_path = Config().SOFTWARE_DIR + "/database/Genome_DB_finish"
        des = os.path.join(db_path, genome_info["bio_mart_annot"])
        des_type = genome_info["biomart_gene_annotype"]
        return des, des_type
    else:
        Config().DBVersion = db_version
        db = Config().get_mongo_client(mtype='ref_rna_v2', db_version=db_version)[Config().get_mongo_dbname('ref_rna_v2', db_version=db_version)]
        collection = db[table_name]
        result = collection.find_one({"task_id": task_id, "type": "origin"})
        species_name = result['species_name']
        json_path = Config().SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        json_dict = get_json()
        des = os.path.join(os.path.split(json_path)[0], json_dict[species_name]["bio_mart_annot"])
        des_type = json_dict[species_name]["biomart_gene_annotype"]
        return des, des_type


def get_pep(result_task, task_id, db_version, table_name="sg_annotation_stat"):
    if "genome_id" in result_task.keys():
        genome_id = result_task["genome_id"]
        db = Config().get_mongo_client(mtype='ref_rna_v2', db_version=db_version)[
            Config().get_mongo_dbname('ref_rna_v2', db_version=db_version)]
        col = db["sg_genome_db"]
        genome_info = col.find_one({"genome_id": genome_id})
        db_path = Config().SOFTWARE_DIR + "/database/Genome_DB_finish"
        pep = os.path.join(db_path, genome_info["pep"])
        return pep
    else:
        Config().DBVersion = db_version
        db = Config().get_mongo_client(mtype='ref_rna_v2', db_version=db_version)[Config().get_mongo_dbname('ref_rna_v2', db_version=db_version)]
        collection = db[table_name]
        result = collection.find_one({"task_id": task_id, "type": "origin"})
        species_name = result['species_name']
        json_path = Config().SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        json_dict = get_json(json_path)
        pep = os.path.join(os.path.split(json_path)[0], json_dict[species_name]["pep"])
        return pep

def get_json(json_path):
    f = open(json_path, "r")
    json_dict = json.loads(f.read())
    return json_dict


class TestFunction(unittest.TestCase):
    def test(self):
        from mbio.workflows.tool_lab.tool_protein_transcript import ToolProteinTranscriptWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'Tool_protein_transcript_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.tool_protein_transcript',
            'options': {
                'protein_type': 'itraq_tmt',
                'protein_task_id': 'tsg_219017',
                'rna_type': 'ref_rna_v2',
                'rna_task_id': 'tsg_38223'
            }
        }
        wsheet = Sheet(data=data)
        wf = ToolProteinTranscriptWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)