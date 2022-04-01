# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import datetime
import json
import os
import types
import unittest
import pandas as pd

from biocluster.api.database.base import report_check
from biocluster.config import Config
from bson.objectid import ObjectId
from bson.son import SON

from api_base import ApiBase


class LncTarget(ApiBase):
    def __init__(self, bind_object):
        super(LncTarget, self).__init__(bind_object)
        self.result_dir = ''
        self.result_file = {}
        self.trans_gene = {}
        self.trans_isgene = {}
        self.task_id = self.bind_object.sheet.id
        self.project_sn = self.bind_object.sheet.project_sn
        self.has_new = True
        self.anno_type = 'origin'
        self.species_name = ""
        self.target_lncrna = dict()
        self.tran2gene = dict()
        self.kegg_json = Config().SOFTWARE_DIR + "/database/KEGG/br08901.json"
        self.annot_type = "ref"
        self.target_name = ""
        # self._db_name = Config().MONGODB + '_ref_rna'

    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                raise Exception('assemble_id必须为ObjectId对象或其对应的字符串！')
        return object_id

    def import_target_cistrans_web(self, target_cistrans_id, cis_target_file, trans_target_file):
        '''
        导入cis把基因信息
        '''
        self.bind_object.logger.info("开始导入靶基因")
        task_id = self.task_id
        # result_dir = os.path.dirname()
        # print result_dir
        '''
        params_dict.update({
            "submit_location": "target_cis",
            "task_id": task_id,
            "task_type":2
        })
        '''

        # self.remove_table_by_main_record(main_table='target_cis', task_id=task_id, detail_table=['target_cis_detail'], detail_table_key='target_cis_id')
        # target_cis_id = self.add_target_cis(params=params_dict, name=None, species_name=species_name)
        self.add_target_cistrans_detail(target_cistrans_id, cis_target_file, trans_target_file)

    def import_diff_set(self, diff_sum_file, category, level, name, source, desc):

        seq_list = list()
        category_list = list()
        kind_list = list()
        exp = pd.read_table(diff_sum_file, header=0, index_col=0)
        print diff_sum_file
        seq_list = list(exp.index)
        category_list = [category for x in seq_list]
        diff_set_id = self.add_geneset(name, level, source, seq_list, category_list, desc)
        return diff_set_id



    def add_geneset(self, name, level, source, seq_list, category_list, desc):
        params = ""
        time_now = datetime.datetime.now()
        type = ", ".join(["{}: {}".format(x, category_list.count(x)) for x in set(category_list)])
        main_dict = {
            'task_id': self.task_id,
            'project_sn': self.project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'params': params,
            'status': 'start',
            'level': level,
            'length': len(seq_list),
            'source': source,
            'type': type,
            'is_use': 1,
            'version': 'v1',
            'desc': desc
        }
        main_id = self.create_db_table('geneset', [main_dict])
        arg_dict = {
            'seq_list': seq_list,
            'category_list': category_list
        }
        self.add_geneset_detail(arg_dict, main_id)
        self.update_db_record('geneset', main_id, insert_dict={'main_id': main_id, 'status': 'end'})
        return main_id

    def add_geneset_detail(self, arg_dict, geneset_id):
        self.create_db_table('geneset_detail', [arg_dict], {'geneset_id': geneset_id})


    def import_target_cistrans(self, params_dict, cis_target, trans_target, filter_lnc_file=None, filter_target_file=None, version="v1.0"):
        '''
        导入cis把基因信息
        '''
        self.bind_object.logger.info("开始导入靶基因")
        task_id = self.task_id
        result_dir = os.path.dirname(cis_target)
        print result_dir

        lncset_id = self.import_diff_set(filter_lnc_file, "lncRNA", "T", "diff_lncRNA_set", "DE_lncR_T_detail", "所有默认比较组获得的差异lncRNA集合")
        targetset_id = self.import_diff_set(filter_target_file, "mRNA", "G", "diff_gene_set", "DE_mR_G_detail", "所有默认比较组获得的差异基因集合")
        params_dict.update({
            "submit_location": "targe_cistrans",
            "task_id": task_id,
            "task_type": 2,
            "lncset_id": str(lncset_id),
            "targetset_id": str(targetset_id)
        })

        # self.remove_table_by_main_record(main_table='target_cistrans', task_id=task_id, detail_table=['target_cistrans_detail'], detail_table_key='target_cistrans_id')
        target_cistrans_id = self.add_target_cistrans(params=params_dict, version=version)
        self.add_target_cistrans_detail(target_cistrans_id, cis_target, trans_target)

    def add_target_cistrans_detail(self, target_cistrans_id, cis_target_file, trans_target_file=None):
        if not isinstance(target_cistrans_id, ObjectId):
            if isinstance(target_cistrans_id, types.StringTypes):
                target_cistrans_id = ObjectId(target_cistrans_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(trans_target_file):
            raise Exception('{}所指定的路径不存在，请检查！'.format(trans_target_file))
        if cis_target_file:
            if not os.path.exists(cis_target_file):
                raise Exception('{}所指定的路径不存在，请检查！'.format(cis_target_file))
        data_list = []

        with open(cis_target_file, 'r') as known_target:
            known_target.readline()
            for line in known_target:
                cols = line.strip().split("\t")

                '''
                if cols[4] == "known":
                    cols[4] = "ref"
                if cols[4] == "novel":
                    cols[4] = "new"
                '''

                data = [('lncrna_id', cols[0]),
                        ('gene_id', cols[1]),
                        ('mgene_id', cols[2]),
                        ('target_type', cols[3]),
                        ('lnc_type', cols[4]),
                        ('gene_name', cols[5]),
                        ('corr', float(cols[6])),
                        ('pvalue', float(cols[7])),
                        ('padjust', float(cols[8])),
                        ('lnc_position', cols[11]),
                        ('target_position', cols[12]),
                        ('location', cols[9]),
                        ('distance', int(cols[10]))
                        ]
                data.append(('target_cistrans_id', target_cistrans_id))
                data = SON(data)
                data_list.append(data)

        self.bind_object.logger.info("读取cis靶基因预测统计信息成功!")

        with open(trans_target_file, 'r') as new_target:
            head_line = new_target.readline()
            for line in new_target:
                cols = line.strip().split("\t")
                '''
                if cols[4] == "known":
                    cols[4] = "ref"
                if cols[4] == "novel":
                    cols[4] = "new"
                '''
                data = [('lncrna_id', cols[0]),
                        ('gene_id', cols[1]),
                        ('mgene_id', cols[2]),
                        ('target_type', cols[3]),
                        ('lnc_type', cols[4]),
                        ('gene_name', cols[5]),
                        ('corr', float(cols[6])),
                        ('pvalue', float(cols[7])),
                        ('padjust', float(cols[8])),
                        ]

                data.append(('target_cistrans_id', target_cistrans_id))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('target_cistrans_detail', data_list)
            self.db['target_cistrans'].update({"_id": target_cistrans_id},
                                              {"$set": {"status": "end", "main_id": target_cistrans_id}})
        except Exception as e:
            self.bind_object.set_error("导入cistrans靶基因预测统计信息出错!")
        else:
            self.bind_object.logger.info("导入cistrans靶基因预测统计信息成功!")

    def import_target_detail_web(self, target_id, new_target_file, known_target_file, params_dict, new_seq=None,
                                 known_seq=None, anno_type="latest", species_name=None, last_id_target=None):
        self.bind_object.logger.info("开始导入靶基因")
        task_id = self.task_id
        self.anno_type = anno_type
        result_dir = os.path.dirname(known_target_file)
        print result_dir

        if last_id_target:
            last_id_target = ObjectId(last_id_target)
        else:
            pass

        if last_id_target:
            self.bind_object.logger.info("删除表格为 {}".format(last_id_target))
            self.remove_table_by_main_record(main_table='target', _id=last_id_target,
                                             detail_table=['target_detail', 'target_stat'],
                                             detail_table_key='target_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        target_id, columns = self.add_target(params=params_dict, name=None, species_name=species_name,
                                             result_dir=result_dir, target_id=target_id)
        self.add_target_detail(target_id, new_target_file, known_target_file, new_seq, known_seq, columns=columns)

    @report_check
    def add_target_cis(self, params, name=None, species_name=None, result_dir=None, target_id=None):
        if target_id:
            if not isinstance(target_id, ObjectId):
                if isinstance(target_id, types.StringTypes):
                    target_id = ObjectId(target_id)
                else:
                    raise Exception('target_id必须为ObjectId对象或其对应的字符串！')
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        if name:
            self.target_name = name
        else:
            self.target_name = 'Target_cis_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': self.target_name,
            'type': self.anno_type,
            'params': params,
            'status': 'start',
            'desc': '注释统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'species_name': species_name,
            'result_dir': result_dir,
            'version': 'v1'
        }

        collection = self.db['target_cis']
        target_cis_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add target_cis!")
        return target_cis_id

    @report_check
    def add_target_cistrans(self, params, name=None, species_name=None, result_dir=None, target_id=None,
                            version="v1.0"):
        if target_id:
            if not isinstance(target_id, ObjectId):
                if isinstance(target_id, types.StringTypes):
                    target_id = ObjectId(target_id)
                else:
                    raise Exception('target_id必须为ObjectId对象或其对应的字符串！')
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        print "params is {}".format(params)
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        if name:
            self.target_name = name
        else:
            self.target_name = 'Target_cistrans_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': self.target_name,
            'type': self.anno_type,
            'params': params,
            'version': version,
            'status': 'start',
            'desc': '注释统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'species_name': species_name,
            'result_dir': result_dir,
            'version': 'v1'
        }

        collection = self.db['target_cistrans']
        target_cistrans_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add target_cistrans!")
        return target_cistrans_id

    @report_check
    def add_target_trans(self, params, name=None, species_name=None, result_dir=None, target_id=None):
        '''
        trans 靶基因主表
        '''
        if target_id:
            if not isinstance(target_id, ObjectId):
                if isinstance(target_id, types.StringTypes):
                    target_id = ObjectId(target_id)
                else:
                    raise Exception('target_id必须为ObjectId对象或其对应的字符串！')
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        if name:
            self.target_name = name
        else:
            self.target_name = 'Target_trans_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': self.target_name,
            'type': self.anno_type,
            'params': params,
            'status': 'start',
            'desc': '注释统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'species_name': species_name,
            'result_dir': result_dir,
            'version': 'v1'
        }

        collection = self.db['target_trans']
        target_trans_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add target_trans!")
        return target_trans_id

    def add_target_geneset(self, new_target_file, known_target_file, diff_summary, lnctype='G', mtype='G',
                           target_type='cis'):
        '''
        插入差异靶基因集
        '''
        diff_lnc_dict = dict()
        diff_target_dict = dict()
        groups = list()

        with open(diff_summary, 'r') as diff_f:
            header = diff_f.readline()
            groups = header.strip("\n").split("\t")[1:-1]
            for group in groups:
                diff_lnc_dict[group] = []
                diff_target_dict[group] = set()
            nums = diff_f.readline()
            for line in diff_f:
                lnc_rna = line.strip("\n").split("\t")[0]
                stats = line.strip("\n").split("\t")[1:-1]
                for n, stat in enumerate(stats):
                    if stat == "yes":
                        diff_lnc_dict[groups[n]].append(lnc_rna)
                    else:
                        pass

        with open(known_target_file, 'r') as known_target:
            known_target.readline()
            for line in known_target:
                cols = line.strip().split("\t")
                for group in groups:
                    if lnctype == 'G':
                        lnc = cols[1]
                    else:
                        lnc = cols[0]
                    if mtype == 'G':
                        m = cols[3]
                    else:
                        m = cols[2]

                    if lnc in diff_lnc_dict[group]:
                        diff_target_dict[group].add(m)

        with open(new_target_file, 'r') as new_target:
            head_line = new_target.readline()
            for line in new_target:
                cols = line.strip().split("\t")
                for group in groups:
                    if lnctype == 'G':
                        lnc = cols[1]
                    else:
                        lnc = cols[0]
                    if mtype == 'G':
                        m = cols[3]
                    else:
                        m = cols[2]

                    if lnc in diff_lnc_dict[group]:
                        diff_target_dict[group].add(m)

        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        for group in groups:
            geneset_main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=group + '_' + target_type + "_targets",
                type=mtype,
                desc='差异lncrna靶基因集',
                group_id="",
                gene_length=len(diff_target_dict[group]),
                is_use=0,
                version='v1'
            )

            geneset_detail_info = [{"seq_list": list(diff_target_dict[group])}]
            self.add_set(geneset_main_info, geneset_detail_info)

    @report_check
    def add_set(self, main_info, detail_info):
        '''
        插入单个基因集
        '''
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        main_info.update(dict(status="start", created_ts=created_ts))
        main_id = self.create_db_table('geneset', [main_info])
        self.create_db_table('geneset_detail', detail_info, tag_dict={"geneset_id": main_id})
        # self.update_db_record('geneset', main_id, status="end", is_use=0, main_id=main_id)
        task_id = main_info['task_id']
        record_dict = {"_id": main_id, "task_id": task_id}
        self.update_db_record('geneset', main_id, query_dict=record_dict, status="end", is_use=0, main_id=main_id,
                              params=task_id)
        return main_id

    @report_check
    def add_target_cis_detail(self, target_cis_id, known_target_file, new_target_file=None):
        if not isinstance(target_cis_id, ObjectId):
            if isinstance(target_cis_id, types.StringTypes):
                target_cis_id = ObjectId(target_cis_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(new_target_file):
            raise Exception('{}所指定的路径不存在，请检查！'.format(new_target_file))
        if new_target_file:
            if not os.path.exists(known_target_file):
                raise Exception('{}所指定的路径不存在，请检查！'.format(known_target_file))
        data_list = []

        with open(known_target_file, 'r') as known_target:
            known_target.readline()
            for line in known_target:
                cols = line.strip().split("\t")

                data = [('lncrna_id', cols[0]),
                        ('gene_id', cols[1]),
                        ('mrna_id', cols[2]),
                        ('mgene_id', cols[3]),
                        ('gene_name', cols[4]),
                        ('chr', cols[5]),
                        ('strand', cols[6]),
                        ('location', cols[7]),
                        ('distance', int(cols[8]))
                        ]
                data.append(('target_cis_id', target_cis_id))
                data.append(('type', 'known'))
                data = SON(data)
                data_list.append(data)

        with open(new_target_file, 'r') as new_target:
            head_line = new_target.readline()
            for line in new_target:
                cols = line.strip().split("\t")

                data = [('lncrna_id', cols[0]),
                        ('gene_id', cols[1]),
                        ('mrna_id', cols[2]),
                        ('mgene_id', cols[3]),
                        ('gene_name', cols[4]),
                        ('chr', cols[5]),
                        ('strand', cols[6]),
                        ('location', cols[7]),
                        ('distance', int(cols[8]))
                        ]
                data.append(('target_cis_id', target_cis_id))
                data.append(('type', 'novel'))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('target_cis_detail', data_list)
            self.db['target_cis'].update({"_id": target_cis_id},
                                         {"$set": {"status": "end", "main_id": target_cis_id}})
        except Exception as e:
            self.bind_object.set_error("导入cis靶基因预测统计信息出错!")
        else:
            self.bind_object.logger.info("导入cis靶基因预测统计信息成功!")

    @report_check
    def add_target_trans_detail(self, target_trans_id, known_target_file, new_target_file=None):
        if not isinstance(target_trans_id, ObjectId):
            if isinstance(target_trans_id, types.StringTypes):
                target_trans_id = ObjectId(target_trans_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')

        if not os.path.exists(new_target_file):
            self.bind_object.logger.info('{}所指定的路径不存在，请检查！'.format(new_target_file))
        if new_target_file:
            if not os.path.exists(known_target_file):
                self.bind_object.logger.info('{}所指定的路径不存在，请检查！'.format(known_target_file))
        data_list = []
        if os.path.exists(known_target_file):
            with open(known_target_file, 'r') as known_target:
                known_target.readline()
                for line in known_target:
                    cols = line.strip().split("\t")
                    data = [('lncrna_id', cols[0]),
                            ('gene_id', cols[1]),
                            ('mrna_id', cols[2]),
                            ('mgene_id', cols[3]),
                            ('gene_name', cols[4]),
                            ('l_start', cols[5]),
                            ('l_end', cols[6]),
                            ('m_start', cols[7]),
                            ('m_end', cols[8]),
                            ('energy', float(cols[9]))
                            ]
                    data.append(('target_trans_id', target_trans_id))
                    data.append(('type', 'known'))
                    data = SON(data)
                    data_list.append(data)

        if os.path.exists(new_target_file):
            with open(new_target_file, 'r') as new_target:
                head_line = new_target.readline()
                for line in new_target:
                    cols = line.strip().split("\t")
                    data = [('lncrna_id', cols[0]),
                            ('gene_id', cols[1]),
                            ('mrna_id', cols[2]),
                            ('mgene_id', cols[3]),
                            ('gene_name', cols[4]),
                            ('l_start', cols[5]),
                            ('l_end', cols[6]),
                            ('m_start', cols[7]),
                            ('m_end', cols[8]),
                            ('energy', float(cols[9]))
                            ]
                    data.append(('target_trans_id', target_trans_id))
                    data.append(('type', 'novel'))
                    data = SON(data)
                    data_list.append(data)

        try:
            self.create_db_table('target_trans_detail', data_list)
            self.db['target_trans'].update({"_id": target_trans_id},
                                           {"$set": {"status": "end", "main_id": target_trans_id}})
        except Exception as e:
            self.bind_object.set_error("导入trans靶基因预测统计信息出错!")
        else:
            self.bind_object.logger.info("导入trans靶基因预测统计信息成功!")

    def run_web(self, target_file, target_file2, annotation_mudule_dir, trans2gene, params_dict, task_id, stat_id,
                last_id, taxon='Animals', exp_level='transcript'):
        """
        annotation_mudule_dir
        annotation_ref_dir
        merge_gen
        merge_tran
        new_anno_path: 新序列注释的结果文件夹
        pfam_path:转录本的pfam_domain
        merge_tran_output: 转录本的merge_annot tool输出结果路径
        merge_gene_output: 基因的merge_annot tool"
        """
        self.bind_object.logger.info("开始导表数据路径为 {}".format(annotation_mudule_dir))
        self.import_target(target_file, target_file2)
        self.set_result_dir(annotation_mudule_dir)
        self.get_trans2gene(trans2gene)
        self.task_id = task_id

        stat_id = ObjectId(stat_id)
        if last_id:
            last_id = ObjectId(last_id)
        else:
            pass
        # stat
        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['ref_stat_path'],
                                        venn_path=self.result_file['ref_venn_path'], seq_type="ref",
                                        exp_level=exp_level)
        self.update_db_record('annotation_stat', stat_id, status="end", main_id=stat_id)
        # query

        if last_id:
            self.bind_object.logger.info("删除表格为 {}".format(last_id))
            self.remove_table_by_main_record(main_table='annotation_stat', _id=last_id,
                                             detail_table=['annotation_stat_detail'], detail_table_key='stat_id')
            self.bind_object.logger.info("删除表格成功 {}".format(last_id))

        query_old_id = self.get_table_by_main_record(main_table='annotation_query', task_id=task_id,
                                                     type=self.anno_type)
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        if self.annot_type == "ref":
            self.add_annotation_query_denovo_detail(query_id=query_id,
                                                    query_path=self.result_file['query_path' + "ref"], seq_type="ref",
                                                    anno_type="T", exp_level=exp_level)
        else:
            self.add_annotation_query_denovo_detail2(query_id=query_id,
                                                     query_path=self.result_file['query_path' + "ref"], seq_type="ref",
                                                     anno_type="T", exp_level=exp_level)
        self.update_db_record('annotation_query', query_id, status="end", main_id=query_id)

        if query_old_id:
            self.remove_table_by_main_record(main_table='annotation_query', _id=query_old_id['_id'],
                                             detail_table=['annotation_query_detail'], detail_table_key='query_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        # go
        cog_old_id = self.get_table_by_main_record(main_table='annotation_cog', task_id=task_id, type=self.anno_type)
        params_select_cog = dict(
            [(k, params_dict.get(k, None)) for k in ('cog_evalue', 'cog_similarity', 'cog_identity')])
        params_select_cog = json.dumps(params_select_cog, sort_keys=True, separators=(',', ':'))
        cog_id = self.add_annotation_cog(name=None, params=params_select_cog, result_dir=self.result_dir)
        if self.annot_type == "ref":
            self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_gene_sum_path' + 'ref'],
                                           seq_type="ref", anno_type="G")
        else:
            self.add_annotation_cog_detail2(cog_id=cog_id, cog_path=self.result_file['n_gene_sum_path' + 'ref'],
                                            seq_type="ref", anno_type="G")
        self.update_db_record('annotation_cog', cog_id, status="end", main_id=cog_id)

        if cog_old_id:
            self.remove_table_by_main_record(main_table='annotation_cog', _id=cog_old_id['_id'],
                                             detail_table=['annotation_cog_detail'], detail_table_key='cog_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        go_old_id = self.get_table_by_main_record(main_table='annotation_go', task_id=task_id, type=self.anno_type)
        params_select_nr = dict([(k, params_dict.get(k, None)) for k in ('nr_evalue', 'nr_similarity', 'nr_identity')])
        params_select_nr = json.dumps(params_select_nr, sort_keys=True, separators=(',', ':'))
        go_id = self.add_annotation_go(name=None, params=params_select_nr, result_dir=self.result_dir)

        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T",
                                        gos_path=self.result_file['gos_path_ref'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G",
                                    gos_path=self.result_file['gene_gos_path_ref'])
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2,
                                     level_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2,
                                      go_path=self.result_file['gene_stat_level2_ref'])
        self.update_db_record('annotation_go', go_id, status="end", main_id=go_id)

        if go_old_id:
            self.remove_table_by_main_record(main_table='annotation_go', _id=go_old_id['_id'],
                                             detail_table=['annotation_go_detail', 'annotation_go_graph',
                                                           'annotation_go_level', 'annotation_go_list'],
                                             detail_table_key='go_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        kegg_old_id = self.get_table_by_main_record(main_table='annotation_kegg', task_id=task_id, type=self.anno_type)
        # kegg
        params_select_kegg = dict(
            [(k, params_dict.get(k, None)) for k in ('kegg_evalue', 'kegg_similarity', 'kegg_identity')])
        params_select_kegg = json.dumps(params_select_kegg, sort_keys=True, separators=(',', ':'))
        kegg_id = self.add_annotation_kegg(name=None, params=None, result_dir=self.result_dir)

        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T",
                                           table_path=self.result_file['table_path_ref'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G",
                                       table_path=self.result_file['gene_table_path_ref'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G",
                                       level_path=self.result_file['gene_pathway_path_ref'],
                                       png_dir=self.result_file['gene_png_path_ref'])

        if kegg_old_id:
            self.remove_table_by_main_record(main_table='annotation_kegg', _id=kegg_old_id['_id'],
                                             detail_table=['annotation_kegg_categories', 'annotation_kegg_level',
                                                           'annotation_kegg_table', 'annotation_kegg_pic'],
                                             detail_table_key='kegg_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

    def stat(self, stat_path):
        with open(stat_path, "rb") as f:
            id_list = []
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                q_id = line[5]
                id_list.append(q_id)
        id_list = list(set(id_list))
        return id_list

    def get_group_dict(self, group_id):
        group_dict = dict()
        document = self.db['specimen_group'].find_one({'main_id': ObjectId(group_id)})
        for group, samples in zip(document['category_names'], document['specimen_names']):
            group_dict[group] = samples
        return group_dict


class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """

    def test_mongo(test):
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            "id": "whole_transcriptome",
            "project_sn": "whole_transcriptome",
            "type": "workflow",
            "name": "whole_transcriptome.whole_transcriptome_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wsheet)

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("whole_transcriptome.lnc_target")
        params_dict = {
            "up_dis": 10,
            "down_dis": 10,
            'group_dict': "{\"A\": [\"A1\", \"A2\"], \"C\": [\"C2\", \"C3\"], \"S\": [\"S1\"\"S3\"]}",
            'group_id': "All",
            "submit_location": "target_cistrans",
            'corr_cutoff': "0.9",
            'corr_way': "pearson",
            'padjust_way': "BH",
            "pvalue_type": "qvalue",
            "pvalue_cutoff": "0.05",
            "task_id": "whole_transcriptome",
            "task_type": 2
        }

        test_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/TargetCistrans/output"
        wf.test_api.import_target_cistrans(
            params_dict,
            os.path.join(test_dir + '/cis_annot.xls'),
            os.path.join(test_dir + '/trans_annot.xls'),
            filter_lnc_file = "/mnt/ilustre/users/sanger-dev/workspace/20191129/WholeTranscriptome_tsg_36305/DiffSplit/output/lncrna/summary.txt",
            filter_target_file = "/mnt/ilustre/users/sanger-dev/workspace/20191129/WholeTranscriptome_tsg_36305/Transfer/output/diff_exp_g/summary.txt"
        )


if __name__ == '__main__':
    unittest.main()
