# -*- coding: utf-8 -*-
# __author__ = 'fengyitong, 20180903'

from biocluster.config import Config
from bson import SON
from bson import ObjectId
import datetime
import os
from biocluster.workflow import Workflow
import sqlite3
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import time

class ToolPrimerWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ToolPrimerWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "relate_name", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "seq_id", "type": "string"},
            {"name": "ref_genome", "type": "string"},
            {"name": "ptt_path", "type": "string"},
            {"name": "genome_path", "type": "string"},
            {"name": "upstream", "type": "int"},
            {"name": "downstream", "type": "int"},
            {'name': 'params', 'type': 'string'},
            {'name': 'tool_type', 'type':'string'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.mongo_db = self.config.get_mongo_client(mtype="prok_rna")[self.config.get_mongo_dbname("prok_rna")]
        # self.bedtool_path = self.config.SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'

    def get_fasta(self):
        conn = sqlite3.connect(self.option('genome_path'))
        cursor = conn.cursor()
        genome2len = dict()
        genome2seq = dict()
        for i in cursor.execute("SELECT genome_id text, genome_seq text, length text from seq_genome"):
            genome2len[i[0]]= int(i[2])
            genome2seq[i[0]]= i[1]

        seq = ''
        # fa_path = os.path.join(self.output_dir, 'seq_updown.fa')
        # bed_path = os.path.join(self.output_dir, 'seq_updown.bed')
        with open(self.option('ptt_path'), 'r') as ptt_r, open(os.path.join(self.output_dir, 'seq.txt'), 'w') as s:
            for line in ptt_r.readlines():
                if line.strip():
                    line = line.strip().split('\t')
                    if line[6] == self.option('seq_id'):
                        if line[2] == '+':
                            start = int(line[1].split('..')[0]) - 1 - self.option('upstream') if int(
                                line[1].split('..')[0]) - 1 - self.option('upstream') > 0 else 0
                            end = int(line[1].split('..')[1]) + self.option('downstream') if int(
                                line[1].split('..')[1]) + self.option('downstream') < genome2len[line[0]] else genome2len[line[0]]
                            seq = genome2seq[line[0]][start:end]
                        if line[2] == '-':
                            start = int(line[1].split('..')[0]) - 1 - self.option('downstream') if int(
                                line[1].split('..')[0]) - 1 - self.option('downstream') > 0 else 0
                            end = int(line[1].split('..')[1]) + self.option('upstream') if int(
                                line[1].split('..')[1]) + self.option('upstream') < genome2len[line[0]] else genome2len[line[0]]
                            seq = genome2seq[line[0]][start:end]
                            seq = Seq(seq).reverse_complement()
                            seq = str(seq)
                        # bed_w.write(line[0] + '\t' + str(start) + '\t' + str(end) + '\t' + line[6] + '\t0' + '\t' + line[2] + '\n')
                        s.write(seq)
                        break

        # cmd = "{bedtool_path} getfasta -fi {fna} -bed {bed} -s -name -fo {fapath}".format(
        #     bedtool_path=self.bedtool_path, fna=self.option('ref_genome'), bed=bed_path, fapath= fa_path)
        # print(cmd)
        # os.system(cmd)
        # return seq
    def run(self):
        self.start_listener()
        self.fire("start")
        self.get_fasta()
        time.sleep(5)
        self.set_db()


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        try:
            main_id = ObjectId(self.option('main_id'))
        except:
            pass
        main_info = dict(
            task_id=self.option('task_id'),
            project_type='prok_rna',
            params=self.option('params'),
            status="end",
            main_id=main_id,
            tool_type=self.option('tool_type'),
            relate_id=main_id,
            relate_name=self.option('relate_name'),
            file_path=self.s3_path
        )
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        main_id_ = db[collection_name].insert_one(SON(main_info)).inserted_id
        conn = db[collection_name]
        task_id = conn.find_one({'_id': main_id_})['task_id']
        conn.update({'_id': main_id_, "task_id": task_id}, {"$set": {'main_id': main_id_}}, upsert=True)
        super(ToolPrimerWorkflow, self).end()


    def set_db(self):
        api = self.api.api("prok_rna.tool_lab_api")
        self.s3_path = '{}/{}'.format(self._sheet.output, 'seq.txt')
        api.add_primer(self.option('main_id'))
        self.set_output()

    def set_output(self):
        self.end()