# -*- coding: utf-8 -*-
# __author__ = 'litangjian'
import csv
import datetime
import json
import os

import pandas as pd
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId

from mbio.api.database.lnc_rna.api_base import ApiBase


class ExpCorrsf(ApiBase):
    def __init__(self, bind_object):
        """数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到

        :param bind_object:
        """
        super(ExpCorrsf, self).__init__(bind_object)
        self._project_type = 'lnc_rna'

    @report_check
    def add_ExpCorrsf(self, upload_dir, work_dir, genes_info=None, main_id=None, corr_way=None,
                      padjust_way=None, params=None, project_sn=None, task_id=None):
        if main_id is None:
            # prepare main table info
            name = "ExpCoranalysis" + '_' + corr_way + '_' + padjust_way + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if isinstance(params, dict):
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                json_dir=work_dir + "/record.json",
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='ExpCorrsf main table',
                params=params,
                status="start",
                exp_level=''
            )
            main_id = self.create_db_table('sg_exp_corrsf', [main_info])

        if isinstance(main_id, (str, bytes, unicode)):
            main_id = ObjectId(main_id)

        genes_dic = None
        if genes_info:
            with open(genes_info) as in_handler:
                genes_dic = json.load(in_handler)

        record_path = os.path.join(work_dir, "record.json")
        with open(record_path) as in_handler:
            data = json.load(in_handler)

        nodes = []
        for i, dic in enumerate(data['nodes']):
            tmp_dict = {
                'node_id': i,
                'seq_id': dic['id'],
                'node_name': dic['name'],
                'group': dic['group'],
                'expcorr_id': main_id
            }
            if genes_dic:
                tmp_dict['type'] = genes_dic[dic['id']]

            nodes.append(tmp_dict)

        self.create_db_table('sg_exp_corrsf_node', nodes)

        # gene_id_1       gene_id_2       cor     p_value q_value
        detail_dic = {}
        with open(work_dir + "/express_correlation_info.xls") as in_handler:
            for dic in csv.DictReader(in_handler, delimiter='\t'):
                new_id_one = dic['gene_id_1'].strip()
                new_id_two = dic['gene_id_2'].strip()
                tmp_dict = {
                    'cor': round(float(dic['cor']), 6),
                    'p_value': round(float(dic['p_value']), 6),
                    'q_value': round(float(dic['q_value']), 6)
                }
                detail_dic[new_id_one + new_id_two] = tmp_dict
                detail_dic[new_id_two + new_id_one] = tmp_dict

        edges = []
        for dic in data['links']:
            source = int(dic['source'])
            source_dic = nodes[source]
            target = int(dic['target'])
            target_dic = nodes[target]
            t_id = target_dic['seq_id'].strip()
            s_id = source_dic['seq_id'].strip()

            corr_type = '-'
            if genes_dic:
                corr_type = genes_dic[s_id] + '-' + genes_dic[t_id]

            tmp_dict = detail_dic.get(t_id + s_id) or detail_dic.get(s_id + t_id, {})
            edges.append({
                "to_id": target,
                "to_name": target_dic['node_name'],
                "to_expcorr": t_id,
                "from_id": source,
                "from_name": source_dic['node_name'],
                "from_expcorr": s_id,
                "distance": dic['distance'],
                'expcorr_id': main_id,
                'cor': tmp_dict.get('cor'),
                'p_value': tmp_dict.get('q_value'),
                'q_value': tmp_dict.get('q_value'),
                'type': corr_type
            })

        self.create_db_table('sg_exp_corrsf_edge', edges)

        self.update_db_record('sg_exp_corrsf', main_id, status="end",
                              main_id=main_id, json_dir=os.path.join(upload_dir, "record.json"))


if __name__ == '__main__':
    anno = ExpCorrsf(None)
    work_dir = '/mnt/ilustre/users/sanger-dev/sg-users/litangjian'
    anno.add_ExpCorrsf(work_dir=work_dir, corr_way='spearman', padjust_way='fdr_bh')
