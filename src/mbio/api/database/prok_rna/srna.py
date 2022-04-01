# !/usr/bin/python
# -*- coding: utf-8 -*-
import types
import os
import json
import unittest
import datetime
import sqlite3
import re
from collections import OrderedDict, defaultdict
from bson.son import SON
from bson.objectid import ObjectId
from mbio.api.database.prok_rna.api_base import ApiBase
from biocluster.api.database.base import Base, report_check
import glob
import copy


class Srna(ApiBase):
    def __init__(self, bind_object):
        super(Srna, self).__init__(bind_object)
        self._project_type = 'prok_rna'
        self.result_dir = self.get_result_dir()

    def get_result_dir(self):
        result_dir = self.bind_object.sheet.output if self.bind_object else None
        if result_dir:
            if result_dir.startswith('tsanger:'):
                result_dir = result_dir.replace('tsanger:','/mnt/ilustre/tsanger-data/')
            else:
                result_dir = result_dir.replace('sanger:','/mnt/ilustre/data/')
        return result_dir

    @report_check
    def add_main_table(self, collection_name, params=None, name=None, desc='sRNA'):
        """
        添加主表的导表函数
        :param collection_name: 主表的collection名字
        :param params: 主表的参数
        :param name: 主表的名字
        :return:
        """
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn if self.bind_object else "prok_rna",
            "task_id": self.bind_object.sheet.id if self.bind_object else "prok_rna",
            "status": "end",
            "name": name if name else collection_name + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'desc': desc,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')) if params else None
        }

        collection = self.db[collection_name]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_predict_rna_bed(self, predict_id, predict_path):
        if not isinstance(predict_id, ObjectId):
            if isinstance(predict_id, types.StringTypes):
                predict_id = ObjectId(predict_id)
            else:
                raise Exception('predict_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(predict_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(predict_path))
        data_list = list()
        length_list = list()
        with open(predict_path, 'r') as f:
            _ = f.readline()
            lines = f.readlines()
            for line in lines:
                line = line.strip('\n').split('\t')
                length_list.append(line[7])
                data = [
                    ('predict_id', predict_id),
                    ('location', line[0]),
                    ('start', int(line[1])),
                    ('end', int(line[2])),
                    ('srna', line[3]),
                    # ('Type', line[4]),
                    ('strand', line[5]),
                    ('antisense', line[6]),
                    ('length', int(line[7])),
                    ('sequence', line[8]),
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_srna_predict_detail']
            collection.insert_many(data_list)
            self.db['sg_srna_predict'].update({"_id": predict_id}, {"$set": {"result_dir": self.result_dir, "main_id": predict_id}})
        except Exception as e:
            if self.bind_object:
                self.bind_object.set_error("导入sRNA预测信息：%s出错!" % (predict_path))
            else:
                pass
        else:
            if self.bind_object:
                self.bind_object.logger.info("导入sRNA预测信息：%s 成功!" % predict_path)
            else:
                pass
#------添加画length图的表-----
        length_list = [int(x) for x in length_list]
        max_length = max(length_list)
        range2length = OrderedDict()
        n = 0
        # while n < max_length:
        #     s = str(n+1) + '-' +str(n+100)
        #     if not s in range2length:
        #         range2length[s] = 0
        #     for le in length_list:
        #         if n <= le <= n+100:
        #             range2length[s] += 1
        #     n += 100
        while n < 500:
            s = str(n + 1) + '-' + str(n + 50)
            if not s in range2length:
                range2length[s] = 0
                for le in length_list:
                    if n + 1 <= le <= n+50:
                        range2length[s] += 1
                n += 50
        over_max = 0
        for le in length_list:
            if le > 500:
                over_max += 1
            range2length.update({'>=501':over_max})
        self.db['sg_srna_predict'].update({"_id": predict_id},
                                              {"$set": {"length": range2length}})

    @report_check
    def add_srna_annot(self, annot_id, annot_path, rfam=False, bsrd=False):
        if not isinstance(annot_id, ObjectId):
            if isinstance(annot_id, types.StringTypes):
                annot_id = ObjectId(annot_id)
            else:
                raise Exception('predict_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(annot_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(annot_path))

        xls_list = list()
        list_list = list()
        for file in os.listdir(annot_path):
            if u'_vs_' in file and u'.xls' in file:
                xls_list.append(annot_path + '/' + file)
            if u'_vs_' in file and u'.list' in file:
                list_list.append(annot_path + '/' + file)
        if bsrd:
            if len(xls_list) != 5 or len(list_list) != 5:   # added BSRD annotation @20210622
                raise Exception('指定的路径下必要的文件不存在，请检查！')
        else:
            if len(xls_list) != 4 or len(list_list) != 4:   # added BSRD annotation @20210622
                raise Exception('指定的路径下必要的文件不存在，请检查！')

        data_list = list()
        for file in xls_list:
            annot_type = os.path.basename(file).strip('.xls').split('_vs_')[-1]
            if rfam and annot_type == 'rfam': #rfam annotation using infernal @20210622
                with open(file, 'r') as f:
                    f.readline()
                    for line in f:
                        line = line.strip('\n').split('\t')
                        data = [
                            ('annot_id', annot_id),
                            ('annot_type', 'rfam'),
                            ('score', float(line[16])),
                            ('evalue', float(line[17])),
                            ('family', line[1]),
                            ('srna', line[3]),
                            ('srna_len', abs(int(line[10]) - int(line[9])) + 1),
                            ('srna_start', int(line[9])),
                            ('srna_end', int(line[10])),
                            ('hit_name', line[2]),
                            ('hit_len', abs(int(line[8]) - int(line[7])) + 1),
                            ('hit_start', int(line[7])),
                            ('hit_end', int(line[8])),
                            ('description', line[-1]),
                        ]
                        data = SON(data)
                        data_list.append(data)
            else:
                with open(file, 'r') as f:
                    # _ = f.readline()
                    # lines = f.readlines()
                    f.readline()
                    for line in f:
                        line = line.strip('\n').split('\t')
                        data = [
                            ('annot_id', annot_id),
                            ('annot_type', annot_type),
                            ('score', float(line[0])),
                            ('evalue', float(line[1])),
                            ('hsp_len', int(line[2])),
                            ('identity', float(line[3])),
                            ('similarity', float(line[4])),
                            ('srna', line[5]),
                            ('srna_len', int(line[6])),
                            ('srna_start', int(line[7])),
                            ('srna_end', int(line[8])),
                            ('hit_name', line[9]),
                            ('hit_len', int(line[10])),
                            ('hit_start', int(line[11])),
                            ('hit_end', int(line[12])),
                            ('description', line[13]),
                        ]
                        data = SON(data)
                        data_list.append(data)
        try:
            collection = self.db['sg_srna_anno_detail']
            # collection.insert_many(data_list)
            self.create_db_table('sg_srna_anno_detail', data_list)
            self.db['sg_srna_anno'].update({"_id": annot_id}, {
                "$set": {"main_id": annot_id}})
        except Exception as e:
            if self.bind_object:
                self.bind_object.set_error("导入sRNA注释信息：%s出错!" % (annot_path))
            else:
                pass
        else:
            if self.bind_object:
                self.bind_object.logger.info("导入sRNA注释信息：%s 成功!" % annot_path)
            else:
                pass

        venn_dict = {'annot_id': annot_id,}
        for file in list_list:
            annot_type = os.path.basename(file).strip('.list').split('_vs_')[-1]
            with open(file, 'r') as f:
                srna_list = f.read().split('\n')
                if annot_type.lower() != 'bsrd' and len(srna_list) > 0:  # modified by zhangyitong on 20210825
                    venn_dict.update({annot_type.lower(): ";".join(srna_list)})
                # self.db['sg_srna_annot'].update({"_id": annot_id}, {"$set": {annot_type + "_list": ";".join(srna_list), annot_type + "_num": len(srna_list)}})
        try:
            self.db['sg_srna_anno_venn'].insert(venn_dict)
        except:
            self.bind_object.set_error("导入sRNA注释venn信息：%s出错!" % (annot_path))
        else:
            self.bind_object.logger.info("导入sRNA注释venn信息：%s 成功!" % annot_path)
        with open(annot_path + '/rfam_stat.xls', 'r') as rfam_s:

            rfam_list = list()
            for line in rfam_s.readlines():
                rfam_dict = {'annot_id': annot_id, 'for_pie': True}
                line = line.strip().split('\t')
                if len(line) > 2:
                    if line[0] == 'All_reads' or line[0] == 'All_matched':
                        rfam_dict.update({
                            'type': line[0],
                            'tn': float(line[1]),
                            'tp': float(line[2].strip('%')),
                            'for_pie': False
                        })
                    else:
                        rfam_dict.update({
                            'type': line[0],
                            'tn': float(line[1]),
                            'tp': float(line[2].strip('%')),
                        })
                    rfam_list.append(rfam_dict)
        try:
            self.db['sg_srna_anno_rfam'].insert_many(rfam_list)
        except:
            self.bind_object.set_error("导入sRNA注释rfam统计信息：%s出错!" % (annot_path))
        else:
            self.bind_object.logger.info("导入sRNA注释rfam统计信息：%s 成功!" % annot_path)
        #             rfam_dict[line[0]] = {'Total Number': line[1], 'Total Percent (%)': line[2].strip('%')}  #line[1] + '&' + line[2]
        # self.db['sg_srna_anno'].update({"_id": annot_id}, {
        #     "$set": {'rfam_info': rfam_dict}})

        stat_list = list()
        with open(annot_path + '/annotation_stat.xls') as annot_stat:
            _ = annot_stat.readline()
            for line in annot_stat.readlines():
                line = line.strip().split('\t')
                if len(line) >2:
                    stat = [
                        ('annot_id', annot_id),
                        ('srna', line[0]),
                        ('tarbase', line[1]),
                        ('map', line[2]),
                        ('rfam', line[3]),
                        ('siphi', line[4]),
                        ('sum', int(line[5])),
                    ]
                    # if bsrd:  # added BSRD annotation @20210622
                    #     stat.extend([('bsrd', line[5]), ('sum', int(line[6]))])
                    # else:
                    #     stat.append(('sum', int(line[5])))
                    stat = SON(stat)
                    stat_list.append(stat)
        try:
            collection = self.db['sg_srna_anno_stat']
            if len(stat_list) > 0:
                collection.insert_many(stat_list)
        except Exception as e:
            self.bind_object.set_error("导入sRNA注释信息：%s出错!" % (annot_path))
        else:
            self.bind_object.logger.info("导入sRNA注释信息：%s 成功!" % annot_path)

    @report_check
    def add_srna_fold(self, fold_id, fold_path):
        if not isinstance(fold_id, ObjectId):
            if isinstance(fold_id, types.StringTypes):
                fold_id = ObjectId(fold_id)
            else:
                raise Exception('predict_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(fold_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(fold_path))
        data_list = list()
        with open(fold_path + '/RNAfold_stat.txt', 'r') as f:
            _ = f.readline()
            lines = f.readlines()
            for line in lines:
                line = line.strip('\n').split('\t')
                data = [
                    ('fold_id', fold_id),
                    ('srna', line[0]),
                    ('mfe', float(line[1])),
                    ('diversity', float(line[2])),
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_srna_fold_detail']
            collection.insert_many(data_list)
            self.db['sg_srna_fold'].update({"_id": fold_id},
                                                  {"$set": {"result_dir": self.result_dir, "main_id": fold_id}})
        except Exception as e:
            if self.bind_object:
                self.bind_object.set_error("导入sRNA二级结构预测信息：%s出错!" % (fold_path))
            else:
                pass
        else:
            if self.bind_object:
                self.bind_object.logger.info("导入sRNA二级结构预测信息：%s 成功!" % fold_path)
            else:
                pass

    @report_check
    def add_srna_target_filter_bydiff(self, target_id, target_path, diff_path):
        diff_set = set()
        for diff in glob.glob(diff_path + '/*.DE.list'):
            with open(diff, 'r') as diff_r:
                for line in diff_r:
                    if line.strip():
                        diff_set.add(line.strip().split('\t')[0])
        if not isinstance(target_id, ObjectId):
            if isinstance(target_id, types.StringTypes):
                target_id = ObjectId(target_id)
            else:
                raise Exception('predict_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(target_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(target_path))

        plex_list = list()
        with open(target_path + '/RNAplex_merge', 'r') as f:
            diff_list = list(diff_set)
            gene_count = defaultdict(int)
            _ = f.readline()
            for line in f:
                line = line.strip('\n').split('\t')
                if len(line) > 2:
                    if line[0] in diff_list or line[1] in diff_list:
                        gene_count[line[0]] += 1 if line[0] in diff_list else 0
                        gene_count[line[1]] += 1 if line[1] in diff_list else 0
                        plex = [
                            ('target_id', target_id),
                            # ('target_type', 'RNAplex'),
                            ('srna', line[0]),
                            ('target', line[1]),
                            # ('alignment', line[2]),
                            ('energy', float(line[3])),
                            ('query_start', int(line[4])),
                            ('query_end', int(line[5])),
                            ('target_start', int(line[6])),
                            ('target_end', int(line[7])),
                        ]
                        plex = SON(plex)
                        plex_list.append(plex)
                        if gene_count[line[0]] == 100:
                            diff_list.remove(line[0])
                        if gene_count[line[1]] == 100:
                            diff_list.remove(line[1])

        try:
            # collection = self.db['sg_srna_target_plex']
            # collection.insert_many(plex_list)
            self.create_db_table('sg_srna_target_plex', plex_list)
        except Exception as e:
            self.bind_object.set_error("导入sRNA靶基因预测信息-plex：%s出错!" % (target_path))
        else:
            self.bind_object.logger.info("导入sRNA靶基因预测信息-plex：%s 成功!" % target_path)

        hybrid_list = list()
        with open(target_path + '/RNAhybrid_merge', 'r') as f:
            diff_list = list(diff_set)
            gene_count = defaultdict(int)
            _ = f.readline()
            for line in f:
                line = line.strip('\n').split('\t')
                if len(line) > 2:
                    if line[0] in diff_list or line[1] in diff_list:
                        gene_count[line[0]] += 1 if line[0] in diff_list else 0
                        gene_count[line[1]] += 1 if line[1] in diff_list else 0
                        hybrid = [
                            ('target_id', target_id),
                            # ('target_type', 'RNAhybrid'),
                            ('srna', line[0]),
                            ('target', line[1]),
                            ('len_mir', int(line[2])),
                            ('len_tar', int(line[3])),
                            ('energy', float(line[4])),
                            ('pvalue', float(line[5])),
                        ]
                        hybrid = SON(hybrid)
                        hybrid_list.append(hybrid)
                        if gene_count[line[0]] == 100:
                            diff_list.remove(line[0])
                        if gene_count[line[1]] == 100:
                            diff_list.remove(line[1])
        try:
            # collection = self.db['sg_srna_target_hybrid']
            # collection.insert_many(hybrid_list)
            self.create_db_table('sg_srna_target_hybrid', hybrid_list)
        except Exception as e:
            self.bind_object.set_error("导入sRNA靶基因预测信息-hybrid：%s出错!" % (target_path))
        else:
            self.bind_object.logger.info("导入sRNA靶基因预测信息-hybrid：%s 成功!" % target_path)

        stat_list = list()
        with open(target_path + '/combine_RNAplex_RNAhybrid', 'r') as comb_r:
            diff_list = list(diff_set)
            gene_count = defaultdict(int)
            _ = comb_r.readline()
            for line in comb_r:
                line = line.strip().split('\t')
                if len(line) > 2:
                    if line[0] in diff_list or line[1] in diff_list:
                        gene_count[line[0]] += 1 if line[0] in diff_list else 0
                        gene_count[line[1]] += 1 if line[1] in diff_list else 0
                        stat = [
                            ('target_id', target_id),
                            ('srna', line[0]),
                            ('target', line[1]),
                            ('database', line[2]),
                        ]
                        stat = SON(stat)
                        stat_list.append(stat)
                        if gene_count[line[0]] == 100:
                            diff_list.remove(line[0])
                        if gene_count[line[1]] == 100:
                            diff_list.remove(line[1])
        try:
            # collection = self.db['sg_srna_target_stat']
            # collection.insert_many(stat_list)
            self.create_db_table('sg_srna_target_stat', stat_list)
            self.db['sg_srna_target'].update({"_id": target_id},
                                           {"$set": {"main_id": target_id}})
        except Exception as e:
            self.bind_object.set_error("导入sRNA靶基因预测统计信息：%s出错!" % (target_path))
        else:
            self.bind_object.logger.info("导入sRNA靶基因预测统计信息：%s 成功!" % target_path)

    @report_check
    def add_srna_target(self, target_id, target_path):
        if not isinstance(target_id, ObjectId):
            if isinstance(target_id, types.StringTypes):
                target_id = ObjectId(target_id)
            else:
                raise Exception('predict_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(target_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(target_path))

        plex_list = list()
        with open(target_path + '/RNAplex_merge', 'r') as f:
            _ = f.readline()
            for row, line in enumerate(f):
                line = line.strip('\n').split('\t')
                plex = [
                    ('target_id', target_id),
                    # ('target_type', 'RNAplex'),
                    ('srna', line[0]),
                    ('target', line[1]),
                    # ('alignment', line[2]),
                    ('energy', float(line[3])),
                    ('query_start', int(line[4])),
                    ('query_end', int(line[5])),
                    ('target_start', int(line[6])),
                    ('target_end', int(line[7])),
                ]
                plex = SON(plex)
                plex_list.append(plex)
                if (row+1) % 100000 ==0:
                    try:
                        self.create_db_table('sg_srna_target_plex', plex_list)
                    except Exception as e:
                        self.bind_object.set_error("导入sRNA靶基因预测信息-plex：%s出错!" % (target_path))
                    else:
                        self.bind_object.logger.info("导入sRNA靶基因预测信息-plex：%s 成功,已经导入了%s行!" % (target_path,str(row)))
                    plex_list = list()
        if plex_list:
            try:
                collection = self.db['sg_srna_target_plex']
                # collection.insert_many(plex_list)
                self.create_db_table('sg_srna_target_plex', plex_list)
            except Exception as e:
                self.bind_object.set_error("导入sRNA靶基因预测信息-plex：%s出错!" % (target_path))
            else:
                self.bind_object.logger.info("导入sRNA靶基因预测信息-plex：%s 成功!" % target_path)

        hybrid_list = list()
        with open(target_path + '/RNAhybrid_merge', 'r') as f:
            _ = f.readline()
            for row, line in enumerate(f):
                line = line.strip('\n').split('\t')
                hybrid = [
                    ('target_id', target_id),
                    # ('target_type', 'RNAhybrid'),
                    ('srna', line[0]),
                    ('target', line[1]),
                    ('len_mir', int(line[2])),
                    ('len_tar', int(line[3])),
                    ('energy', float(line[4])),
                    ('pvalue', float(line[5])),
                ]
                hybrid = SON(hybrid)
                hybrid_list.append(hybrid)
                if (row+1) % 100000 == 0:
                    try:
                        self.create_db_table('sg_srna_target_hybrid', hybrid_list)
                    except Exception as e:
                        self.bind_object.set_error("导入sRNA靶基因预测信息-hybrid：%s出错!" % (target_path))
                    else:
                        self.bind_object.logger.info("导入sRNA靶基因预测信息-hybrid：%s 成功!,已经导入了%s行!" % (target_path,str(row)))
                    hybrid_list = list()
        if hybrid_list:
            try:
                collection = self.db['sg_srna_target_hybrid']
                # collection.insert_many(hybrid_list)
                self.create_db_table('sg_srna_target_hybrid', hybrid_list)
            except Exception as e:
                self.bind_object.set_error("导入sRNA靶基因预测信息-hybrid：%s出错!" % (target_path))
            else:
                self.bind_object.logger.info("导入sRNA靶基因预测信息-hybrid：%s 成功!" % target_path)

        stat_list = list()
        with open(target_path + '/combine_RNAplex_RNAhybrid', 'r') as f:
            _ = f.readline()
            for row, line in enumerate(f):
                line = line.strip().split('\t')
                stat = [
                    ('target_id', target_id),
                    ('srna', line[0]),
                    ('target', line[1]),
                    ('database', line[2]),
                ]
                stat = SON(stat)
                stat_list.append(stat)
                if (row + 1) % 100000 == 0:
                    try:
                        self.create_db_table('sg_srna_target_stat', stat_list)
                        self.db['sg_srna_target'].update({"_id": target_id},
                                                         {"$set": {"main_id": target_id}})
                    except Exception as e:
                        self.bind_object.set_error("导入sRNA靶基因预测统计信息：%s出错!" % (target_path))
                    else:
                        self.bind_object.logger.info("导入sRNA靶基因预测统计信息：%s 成功,已经导入了%s行!" % (target_path,str(row)))
                    stat_list = list()
        if stat_list:
            try:
                collection = self.db['sg_srna_target_stat']
                # collection.insert_many(stat_list)
                self.create_db_table('sg_srna_target_stat', stat_list)
                self.db['sg_srna_target'].update({"_id": target_id},
                                               {"$set": {"main_id": target_id}})
            except Exception as e:
                self.bind_object.set_error("导入sRNA靶基因预测统计信息：%s出错!" % (target_path))
            else:
                self.bind_object.logger.info("导入sRNA靶基因预测统计信息：%s 成功!" % target_path)

    def run(self, predict_path, annot_path, fold_path, target_path, params):
        predict_id = self.add_main_table('sg_srna_predict', params=params)
        self.add_predict_rna_bed(predict_id, predict_path)

        annot_id = self.add_main_table('sg_srna_anno', params=params)
        self.add_srna_annot(annot_id, annot_path)

        fold_id = self.add_main_table('sg_srna_fold', params=params)
        self.add_srna_fold(fold_id, fold_path)

        target_id = self.add_main_table('sg_srna_target', params=params)
        self.add_srna_target(target_id, target_path)

#-------添加将碱基与蛋白序列文件建mysql库的操作-------20180901
    def parse_seq_file(self, seq_file, match=re.compile(r'>([^\s]+)').match):
        """
        generator for parsing sequence file
        :param seq_file: fasta sequence file
        :param match: regexp pattern for parsing line startswith '>'
        :return: seq_tag, sequence
        """
        with open(seq_file, 'r') as f:
            j = 0
            seq_tag = ""
            sequence = ""
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    j += 1
                    if j > 1:
                        yield seq_tag, sequence
                        #print seq_tag, sequence
                    seq_tag = match(line).groups()[0]
                    sequence = ''
                else:
                    sequence += line.strip()
            else:
                yield seq_tag, sequence

    def build_seq_database(self, db_path, file, table_name = 'seq_gene', type_ = 'gene'):
        """
        build sequence db
        """
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        # add seq detail
        cursor.execute('DROP TABLE IF EXISTS {}'.format(table_name))
        cursor.execute('CREATE TABLE {} ({}_id text, {}_seq text, length text)'.format(table_name, type_, type_))
        match = re.compile(r'>(.*?)\s+.*').match
        parser = self.parse_seq_file(file, match)
        for par in parser:
            #print pep[0]
            par_id = par[0]
            par_seq = par[1]
            cursor.execute("INSERT INTO {} VALUES (\"{}\",\"{}\",\"{}\")".format(
                table_name, par_id, par_seq, str(len(par_seq))))
        conn.commit()
        conn.close()

class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """
    def test_mongo(test):
        from mbio.workflows.prok_rna.prokrna_test_api import ProkrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            "id": "prok_rna_srna",
            #+ str(random.randint(1,10000)),
            #"id": "denovo_rna_v2",
            "project_sn": "prok_rna_srna",
            #+ str(random.randint(1,10000)),
            "type": "workflow",
            "name": "prok_rna.prokrna_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = ProkrnaTestApiWorkflow(wsheet)

        predict_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180801/Single_Srna_2589_fyt/Srna/output/rockhopper/genome.predicted_RNA.bed.xls"
        fold_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180801/Single_Srna_2589_fyt/Srna/output/srna_fold"
        annot_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180801/Single_Srna_2589_fyt/Srna/output/srna_annot"
        target_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180801/Single_Srna_2589_fyt/Srna/output/srna_target"
        db_path = '/mnt/ilustre/users/sanger-dev/workspace/20180831/Single_Srna_8837_fyt/Srna/output/rockhopper'
        genome_path = '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/pipline/ref/GCF_000009345.1_ASM934v1_genomic.fna'

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("prok_rna.srna")
        params = {
            "name": 'try_all',
            "database": 'all',
            "software": 'RNAplex,RNAhybrid',
        }

        # predict_id = wf.test_api.add_main_table('sg_srna_predict_bed', params=params)
        # wf.test_api.add_predict_rna_bed(predict_id, predict_dir)
        #
        # annot_id = wf.test_api.add_main_table('sg_srna_annot', params=params)
        # wf.test_api.add_srna_annot(annot_id, annot_dir)
        #
        # fold_id = wf.test_api.add_main_table('sg_srna_fold', params=params)
        # wf.test_api.add_srna_fold(fold_id, fold_dir)
        #
        # target_id = wf.test_api.add_main_table('sg_srna_target', params=params)
        # wf.test_api.add_srna_target(target_id, target_dir)

        #wf.test_api.run(predict_dir, annot_dir, fold_dir, target_dir, params)
        # wf.test_api.build_seq_database(db_path + '/cds.fa.db.sqlite3', db_path + '/cds.fa', table_name = 'seq_gene', type_ = 'gene')
        # wf.test_api.build_seq_database(db_path + '/cds.faa.db.sqlite3', db_path + '/cds.faa', table_name = 'seq_protein', type_ = 'protein')
        wf.test_api.build_seq_database(db_path + '/genome_fna.db.sqlite3', genome_path, table_name = 'seq_genome', type_ = 'genome')

if __name__ == '__main__':
    unittest.main()
