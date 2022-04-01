# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from collections import OrderedDict
import os
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import json
from mbio.api.database.denovo_rna_v2.api_base import ApiBase


class Ssr(ApiBase):
    def __init__(self, bind_object):
        super(Ssr, self).__init__(bind_object)
        self._project_type = 'denovo_rna_v2'

    @report_check
    def add_ssr(self, ssr_statistic_path, ssr_detail_path, ssr_class_path, name=None, params=None, stat_id=None, project_sn='denovo_rna_v2', task_id='denovo_rna_v2'):
        ssr_statistic_dict = OrderedDict()
        ssr_statistic_dict_p = OrderedDict()
        # params = {"rep_1": 10, "rep_2":6, "rep_3": 5, "rep_4":5, "rep_5": 5, "rep_6":5, "distance":100}
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        with open(ssr_statistic_path, 'r') as f2:
            header2 = f2.readline()
            for line in f2:
                line = line.strip().split('\t')
                ssr_statistic_dict[line[0]] = int(line[1])

        with open(ssr_statistic_path, 'r') as f2:
            for i in range(4):
                read_over = f2.readline()
            for line in f2:
                line = line.strip().split('\t')
                ssr_statistic_dict_p[line[0] + '_p'] = float("%.6f" % float(line[2]))

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'Ssr_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'start',
            'desc': 'ssr结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'ssr_statistic_dict':ssr_statistic_dict,
            'ssr_statistic_dict_p': ssr_statistic_dict_p
        }

        ssr_id = self.create_db_table('sg_ssr', [insert_data])
        self.add_ssr_detail(ssr_id, ssr_detail_path)
        self.add_ssr_class(ssr_id, ssr_class_path)
        self.update_db_record('sg_ssr', ssr_id, status="end", main_id=ObjectId(ssr_id))
        return ssr_id

    @report_check
    def add_ssr_detail(self, ssr_id, ssr_detail_path):
        if not isinstance(ssr_id, ObjectId):
            if isinstance(ssr_id, types.StringTypes):
                ssr_id = ObjectId(ssr_id)
            else:
                self.bind_object.set_error('ssr_id必须为ObjectId对象或其对应的字符串！', code="52001101")
        if not os.path.exists(ssr_detail_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(ssr_detail_path), code="52001102")
        data_list = []
        with open(ssr_detail_path, 'r') as f1:
            f1.readline()
            for line in f1:
                line = line.strip().split('\t')
                data = [
                    ('ssr_id', ssr_id),
                    ('seq_id', line[0]),
                    ('ssr_nr', int(line[1])),
                    ('ssr_type', line[2]),
                    ('ssr_motif', line[3]),
                    ('ssr_len', int(line[4])),
                    ('start', int(line[5])),
                    ('end', int(line[6])),

                    ('fp1', line[7] if line[7] else " "),
                    ('tm1_1', line[8] if line[8] else " "),
                    ('size_1_1', int(float(line[9])) if (line[9]) else " "),
                    ('rp1', line[10] if line[10] else " "),
                    ('tm1_2', line[11] if line[11] else " "),
                    ('size_1_2', int(float(line[12])) if line[12] else " "),
                    ('pro1_size', int(float(line[13])) if line[13] else " "),
                    ('start_1', int(float(line[14])) if line[14] else " "),
                    ('end_1', int(float(line[15])) if line[15] else " "),

                    ('fp2', line[16] if line[16] else " "),
                    ('tm2_1', line[17] if line[17] else " "),
                    ('size_2_1', int(float(line[18])) if line[18] else " "),
                    ('rp2', line[19] if line[19] else " "),
                    ('tm2_2', line[20] if line[20] else " "),
                    ('size_2_2', int(float(line[21])) if line[21] else " "),
                    ('pro2_size', int(float(line[22])) if line[22] else " "),
                    ('start_2', int(float(line[23])) if line[23] else " "),
                    ('end_2', int(float(line[24])) if line[24] else " "),

                    ('fp3', line[25] if line[25] else " "),
                    ('tm3_1', line[26] if line[26] else " "),
                    ('size_3_1', int(float(line[27])) if line[27] else " "),
                    ('rp3', line[28] if line[28] else " "),
                    ('tm3_2', line[29] if line[29] else " "),
                    ('size_3_2', int(float(line[30])) if line[30] else " "),
                    ('pro3_size', int(float(line[31])) if line[31] else " "),
                    ('start_3', int(float(line[32])) if line[32] else " "),
                    ('end_3', int(float(line[33])) if line[33] else " "),

                    ('position', line[34])
                ]
                data = SON(data)
                data_list.append(data)

        self.create_db_table('sg_ssr_detail', data_list)


    def add_ssr_class(self, ssr_id, ssr_class_path):
        if not isinstance(ssr_id, ObjectId):
            if isinstance(ssr_id, types.StringTypes):
                ssr_id = ObjectId(ssr_id)
            else:
                self.bind_object.set_error('ssr_id必须为ObjectId对象或其对应的字符串！', code="52001103")
        if not os.path.exists(ssr_class_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(ssr_class_path), code="52001104")

        def type_convert(x):
            if len(x) == 3:
                return 'Mono-nucleotide repeats'
            elif len(x) == 5:
                return 'Di-nucleotide repeats'
            elif len(x) == 7:
                return 'Tri-nucleotide repeats'
            elif len(x) == 9:
                return 'Tetra-nucleotide repeats'
            elif len(x) == 11:
                return 'Penta-nucleotide repeats'
            elif len(x) == 13:
                return 'Hexa-nucleotide repeats'

        data_list = []
        with open(ssr_class_path, 'r') as f3:
            # strip会去除第二行多余的\t，所以第二行会出现超出索引范围
            header3_1 = f3.readline()
            header3_2= f3.readline()
            for line in f3:
                line = line.strip().split('\t')
                data = [
                    ('ssr_id', ssr_id),
                    ('ssr_type', line[0]),
                    ('ssr_catlogy', type_convert(line[0])),
                    ('rep_num_1_5', int(line[1])),
                    ('rep_num_6_10', int(line[2])),
                    ('rep_num_11_15', int(line[3])),
                    ('rep_num_up_15', int(line[4]))

                   ]
                data = SON(data)
                data_list.append(data)

        self.create_db_table('sg_ssr_class', data_list)

if __name__ == '__main__':
    anno = Ssr(None)
    ssr_statistic_path = '/mnt/ilustre/users/sanger-dev/workspace/20171206/Single_test_ssr49650/Ssr/ssr_type.txt'
    ssr_detail_path = '/mnt/ilustre/users/sanger-dev/workspace/20171206/Single_test_ssr49650/Ssr/tmp.txt'
    ssr_class_path = '/mnt/ilustre/users/sanger-dev/workspace/20171206/Single_test_ssr49650/Ssr/ssr_repeats_class.txt'
    ssr_id = anno.add_ssr(ssr_statistic_path = ssr_statistic_path, ssr_detail_path=ssr_detail_path, ssr_class_path=ssr_class_path ,name=None, params=None)
