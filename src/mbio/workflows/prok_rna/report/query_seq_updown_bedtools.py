# -*- coding: utf-8 -*-
# __author__ = 'fengyitong, 20180903'


from bson import ObjectId
import datetime
import os
from biocluster.workflow import Workflow
import sqlite3

class QuerySeqUpdownBedtoolsWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(QuerySeqUpdownBedtoolsWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "seq_id", "type": "string"},
            {"name": "ref_genome", "type": "string"},
            {"name": "ptt_path", "type": "string"},
            {"name": "genome_path", "type": "string"},
            {"name": "upstream", "type": "int"},
            {"name": "downstream", "type": "int"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.mongo_db = self.config.get_mongo_client(mtype="prok_rna")[self.config.get_mongo_dbname("prok_rna")]
        self.bedtool_path = self.config.SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'

    def get_fasta(self):
        conn = sqlite3.connect(self.option('genome_path'))
        cursor = conn.cursor()
        genome2len = dict()
        for i in cursor.execute("SELECT genome_id text, length text from seq_genome"):
            genome2len[i[0]]= int(i[1])

        fa_path = os.path.join(self.output_dir, 'seq_updown.fa')
        bed_path = os.path.join(self.output_dir, 'seq_updown.bed')
        with open(self.option('ptt_path'), 'r') as ptt_r, \
            open(bed_path, 'w') as bed_w:
            for line in ptt_r.readlines():
                if line.strip():
                    line = line.strip().split('\t')
                    if line[6] == self.option('seq_id'):
                        if line[2] == '+':
                            start = int(line[1].split('..')[0]) - 1 - self.option('upstream') if int(
                                line[1].split('..')[0]) - 1 - self.option('upstream') > 0 else 0
                            end = int(line[1].split('..')[1]) + self.option('downstream') if int(
                                line[1].split('..')[1]) + self.option('downstream') < genome2len[line[0]] else genome2len[line[0]]
                        if line[2] == '-':
                            start = int(line[1].split('..')[0]) - 1 - self.option('downstream') if int(
                                line[1].split('..')[0]) - 1 - self.option('downstream') > 0 else 0
                            end = int(line[1].split('..')[1]) + self.option('upstream') if int(
                                line[1].split('..')[1]) + self.option('upstream') < genome2len[line[0]] else genome2len[line[0]]
                        bed_w.write(line[0] + '\t' + str(start) + '\t' + str(end) + '\t' + line[6] + '\t0' + '\t' + line[2] + '\n')
                        break

        cmd = "{bedtool_path} getfasta -fi {fna} -bed {bed} -s -name -fo {fapath}".format(
            bedtool_path=self.bedtool_path, fna=self.option('ref_genome'), bed=bed_path, fapath= fa_path)
        print(cmd)
        os.system(cmd)
        return fa_path
    def run(self):
        self.start_listener()
        self.fire("start")
        fa_path = self.get_fasta()
        self.set_db(fa_path)
        self.end()

    def end(self):
        super(QuerySeqUpdownBedtoolsWorkflow, self).end()

    def set_db(self, fa_path):
        with open(fa_path, 'r') as fa_r:
            fa_info = fa_r.read()
            gene_seq = ''.join(fa_info.split('\n')[1:])
        conn = self.mongo_db['sg_query_seq_down']
        conn.update({'_id': ObjectId(self.option('main_id'))}, {'$set': {'gene_seq': gene_seq}})
        all_exp = self.api.api("prok_rna.all_exp")
        all_exp.add_venn_tt(
            self.option('main_id'),
        )