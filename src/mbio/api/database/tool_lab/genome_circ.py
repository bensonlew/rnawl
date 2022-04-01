# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20210202
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
import glob
import re
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from Bio import SeqIO


class GenomeCirc(ApiBase):
    def __init__(self, bind_object):
        super(GenomeCirc, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_genome_circ(self, main_id, genome_circ, genome_type, genome_fa, genome_struction, params=None, project_sn='tool_lab', task_id='tool_lab',):
        # add main table info
        if main_id is None:
            name = "Genome_circ" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Genome_circ',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_genome_circ', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        chr_list = list()
        text_detail = list()
        if genome_type == 'finish':
            genome_fasta = SeqIO.parse(genome_fa, 'fasta')
            for i in genome_fasta:
                name = i.name
                seq = i.seq
                chr_list.append(i.name)
                text_detail.append({'text': '{},{}'.format(name, str(len(seq))), 'circ_id': main_id, 'name': name,
                                    'x': "", 'y': '', 'group': name, 'type': 'text'})
            self.create_db_table('sg_genome_circ_text_detail', text_detail)
            text_dict = {"name": "name", "condition": {"type": "text"}}
            text_json = json.dumps(text_dict, sort_keys=True, separators=(',', ':'))
            self.update_db_record('sg_genome_circ', main_id, status='end', main_id=main_id, text_data=text_json)

            df = pd.read_table(genome_circ, header=0, sep='\t')
            df.fillna('', inplace=True)
            df['type'] = 'c_arrows'
            df['circ_id'] = main_id
            df['symbol'] = 'arrow2'
            circ_detail = df.to_dict('r')
            self.create_db_table('sg_genome_circ_detail', circ_detail)
            dict_a = {"chr": "chr", "start": "start", "end": "end", "direction": "direction", "group": "group", "label": "label", "symbol": "symbol", "condition": {"type": "c_arrows"}}
            dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
            chr_list.sort()
            chr_list = json.dumps(chr_list)
            self.update_db_record('sg_genome_circ', main_id, main_id=main_id, c_arrows_data=dict_b, chr_list=chr_list)
        if genome_type == 'scan':

            name = os.path.basename(genome_fa).strip().split('.')[0]
            chr_list = [name]
            genome_fasta = SeqIO.parse(genome_fa, 'fasta')
            length_total = 0
            for i in genome_fasta:
                length_total += len(i.seq)
            text_detail = [{'text': '{},{}'.format(name, str(length_total)), 'circ_id': main_id, 'name': name,
                            'x': "", 'y': '', 'group': name, 'type': 'text'}]
            self.create_db_table('sg_genome_circ_text_detail', text_detail)
            text_dict = {"name": "name", "condition": {"type": "text"}}
            text_json = json.dumps(text_dict, sort_keys=True, separators=(',', ':'))
            self.update_db_record('sg_genome_circ', main_id, status='end', main_id=main_id, text_data=text_json)

            df = pd.read_table(genome_circ, header=0, sep='\t')
            df.fillna('', inplace=True)
            df['type'] = 'c_arrows'
            df['circ_id'] = main_id
            df['symbol'] = 'arrow2'
            df['chr'] = name
            circ_detail = df.to_dict('r')
            self.create_db_table('sg_genome_circ_detail', circ_detail)
            dict_a = {"chr": "chr", "start": "start", "end": "end", "direction": "direction", "group": "group", "label": "label", "symbol": "symbol", "condition": {"type": "c_arrows"}}
            dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
            chr_list = json.dumps(chr_list, ensure_ascii=False)
            self.update_db_record('sg_genome_circ', main_id, main_id=main_id, c_arrows_data=dict_b, chr_list=chr_list)


