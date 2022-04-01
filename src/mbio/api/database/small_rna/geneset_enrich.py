# -*- coding: utf-8 -*-
# __author__ = 'qindanhua, qinjincheng'

from mbio.api.database.small_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId
import types
import os
from bson.son import SON
import math
import re

class GenesetEnrich(ApiBase):
    def __init__(self, bind_object):
        super(GenesetEnrich, self).__init__(bind_object)

    @report_check
    def add_go_enrich_detail(self, go_enrich_id, go_enrich_dir):
        """
        GO富集详情导表函数
        :param go_enrich_id: 主表ID
        :param go_enrich_dir: 结果文件（不是文件夹）
        :return:
        """
        if not isinstance(go_enrich_id, ObjectId):
            if isinstance(go_enrich_id, types.StringTypes):
                go_enrich_id = ObjectId(go_enrich_id)
            else:
                raise Exception('go_enrich_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(go_enrich_dir):
            raise Exception('{}所指定的路径不存在。请检查！'.format(go_enrich_dir))
        data_list = []
        go_type = []
        with open(go_enrich_dir, 'r') as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                if line[2] == 'p':
                    continue
                data = [
                    ('go_enrich_id', go_enrich_id),
                    ('go_id', line[0]),
                    ('go_type', line[1]),
                    ('enrichment', line[2]),
                    ('discription', line[3]),
                    ('ratio_in_study', line[4]),
                    ('ratio_in_pop', line[5]),
                    ('p_uncorrected', float(line[6])),
                    ('p_corrected', float(line[9])),
                    ('enrich_factor', float(line[4].split("/")[0]) / float(line[5].split("/")[0])),
                    ('depth', int(line[7])),
                    ('study_count', int(line[4].split("/")[0])),
                    ('pop_count', int(line[5].split("/")[0])),
                    ('seq_list', line[-1]),
                    ('seq_str', line[-1].split(";"))
                ]
                go_type.append(line[1])
                data = SON(data)
                data_list.append(data)
        if data_list:
            # 插入-logpvalue -logpcorrected 值相关字段
            pvalues = [dict(son)['p_uncorrected'] for son in data_list]
            if len([x for x in pvalues if x > 0]) > 0:
                pvalues_min = min([x for x in pvalues if x > 0]) / 10
            else:
                pvalues_min = 0.0001
            pvalues_min = - math.log10(pvalues_min)
            log10x = [-math.log10(x) if x > 0 else pvalues_min for x in pvalues]
            for i in range(0, len(log10x)):
                data_list[i]['neg_log10p_uncorrected'] = log10x[i]

            pvalues = [dict(son)['p_corrected'] for son in data_list]
            if len([x for x in pvalues if x > 0]) > 0:
                pvalues_min = min([x for x in pvalues if x > 0]) / 10
            else:
                pvalues_min = 0.0001
            pvalues_min = - math.log10(pvalues_min)
            log10x = [-math.log10(x) if x > 0 else min(pvalues) / 10 for x in pvalues]
            for i in range(0, len(log10x)):
                data_list[i]['neg_log10p_corrected'] = log10x[i]

            try:
                collection = self.db['sg_geneset_go_enrich_detail']
                collection.insert_many(data_list)
                coll = self.db['sg_geneset_go_enrich']
                go_type = list(set(go_type))
                coll.update({'_id': go_enrich_id}, {'$set': {'categories': go_type}})
                # main_collection = self.db['sg_geneset_go_enrich']
                # main_collection.update({"_id": ObjectId(go_enrich_id)}, {"$set": {"status": "end"}})
            except Exception as e:
                print("导入go富集信息：%s出错:%s" % (go_enrich_dir, e))
            else:
                print("导入go富集信息：%s成功!" % (go_enrich_dir))
        else:
            raise Exception('GO富集没有结果')
        # remove needless line of which col3 is p
        self.bind_object.logger.debug('start modifying {} for removing needless line'.format(go_enrich_dir))
        need_line_list = list()
        with open(go_enrich_dir) as fr:
            need_line_list.append(fr.readline())
            for line in fr:
                if line.split('\t')[2] == 'e':
                    need_line_list.append(line)
        with open(go_enrich_dir, 'w') as fw:
            fw.writelines(need_line_list)
        self.bind_object.logger.debug('succeed in modifying {}'.format(go_enrich_dir))

    @report_check
    def add_kegg_enrich_detail(self, kegg_enrich_id, kegg_enrich_table, geneset_list_path, all_list_path):
        """
        KEGG富集详情表导表函数
        :param kegg_enrich_id: 主表id
        :param kegg_enrich_table: 结果表
        :return:
        """
        if not isinstance(kegg_enrich_id, ObjectId):
            if isinstance(kegg_enrich_id, types.StringTypes):
                kegg_enrich_id = ObjectId(kegg_enrich_id)
            else:
                raise Exception('kegg_enrich_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(kegg_enrich_table):
            raise Exception('kegg_enrich_table所指定的路径:{}不存在，请检查！'.format(kegg_enrich_table))
        data_list = []
        # geneset_length = len(open(geneset_list_path, "r").readlines())
        # all_list_length = len(open(all_list_path, "r").readlines())
        kegg_type1 = []
        with open(kegg_enrich_table, 'rb') as r:
            for line in r:
                if re.match(r'\w', line):
                    line = line.strip('\n').split('\t')
                    insert_data = {
                        'kegg_enrich_id': kegg_enrich_id,
                        'term': line[1],
                        'database': line[2],
                        'id': line[3].split("path:")[1] if "path:" in line[3] else line[3],
                        'discription': line[1],
                        'study_count': int(line[0]),
                        "background_number": line[5].split("/")[1],
                        'ratio_in_study': line[4],
                        'ratio_in_pop': line[5],
                        'enrich_factor': float(line[0])/float(line[5].split("/")[0]),
                        'pvalue': float(line[6]),
                        'corrected_pvalue': float(line[7]) if not line[7] == "None" else "None",
                        'seq_list': line[8],
                        'hyperlink': line[9],
                        'kegg_type': "".join([x[0] for x in line[11].split(' ')])
                    }
                    kegg_type1.append("".join([x[0] for x in line[11].split(' ')]))
                    data_list.append(insert_data)
            if data_list:
                # 插入-logpvalue -logpcorrected 值相关字段
                pvalues = [dict(son)['pvalue'] for son in data_list]
                if len([x for x in pvalues if x > 0]) > 0:
                    pvalues_min = min([x for x in pvalues if x > 0])/10
                else:
                    pvalues_min = 0.0001
                pvalues_min = - math.log10(pvalues_min)
                log10x = [-math.log10(x) if x>0 else pvalues_min for x in pvalues]
                for i in range(0,len(log10x)):
                    data_list[i]['neg_log10p_uncorrected'] = log10x[i]

                pvalues = [dict(son)['corrected_pvalue'] for son in data_list]
                if len([x for x in pvalues if x > 0]) > 0:
                    pvalues_min = min([x for x in pvalues if x > 0])/10
                else:
                    pvalues_min = 0.0001
                pvalues_min = - math.log10(pvalues_min)

                log10x = [-math.log10(x) if x>0 else pvalues_min for x in pvalues]
                for i in range(0,len(log10x)):
                    data_list[i]['neg_log10p_corrected'] = log10x[i]

                kegg_type1 = list(set(kegg_type1))
                try:
                    collection = self.db['sg_geneset_kegg_enrich_detail']
                    collection.insert_many(data_list)
                    coll = self.db['sg_geneset_kegg_enrich']

                    coll.update({'_id': kegg_enrich_id}, {'$set': {'categories': kegg_type1}})
                    # main_collection = self.db['sg_geneset_kegg_enrich']
                    # main_collection.update({"_id": ObjectId(enrich_id)}, {"$set": {"status": "end"}})
                except Exception as e:
                    self.bind_object.set_error("导入kegg富集统计表：%s信息出错:%s" % (kegg_enrich_table, e))
                else:
                    self.bind_object.logger.info("导入kegg富集统计表:%s信息成功!" % kegg_enrich_table)
            else:
                coll = self.db['sg_geneset_kegg_enrich']
                coll.update({'_id': kegg_enrich_id}, {'$set': {'desc': 'no_result'}})
                # self.bind_object.logger.info("kegg富集统计表没结果：" % kegg_enrich_table)
                raise Exception("kegg富集统计表没结果")

    @report_check
    def add_kegg_enrich_pic(self, main_table_id, level_path, png_dir):
        # 导入图片信息数据
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                kegg_id = ObjectId(main_table_id)
            else:
                raise Exception('main_table_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(level_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(level_path))
        if not os.path.exists(png_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(png_dir))
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                pid = re.sub('path:', '', line[3])
                if os.path.exists(png_dir + '/' + line[3] + '.html.mark'):
                    with open(png_dir + '/' + line[3] + '.html.mark', 'r') as mark_f:
                        for line_mark in mark_f.readlines():
                            # print len(line_mark.strip("\n").split("\t"))
                            if len(line_mark.strip("\n").split("\t")) == 10:
                                [png, shape, bg_color, fg_color, coords, title, kos, href, gene_list, geneset_list] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                            else:
                                continue

                            insert_data = {
                                'kegg_enrich_id': kegg_id,
                                'pathway_id': line[3],
                                'shape': shape,
                                'bg_colors': bg_color,
                                'fg_colors': fg_color,
                                'coords': coords,
                                'href': href,
                                'kos': kos,
                                'title': title,
                                'gene_list': gene_list.split("|"),
                                'geneset_list': geneset_list.split("|")
                            }

                            if bg_color != "" and len(bg_color.split(",")) > 0:
                                insert_data.update({'bg_type': len(bg_color.split(","))})
                            if fg_color != "" and len(fg_color.split(",")) > 0:
                                insert_data.update({'fg_type': len(fg_color.split(","))})
                            data_list.append(insert_data)
                else:
                    self.bind_object.logger.info("kegg 图片{} 不存在html标记!".format(line[3]))

        if data_list:
            try:
                collection = self.db['sg_geneset_kegg_enrich_pic']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入kegg注释图片信息：%s、%s出错!" % (level_path, png_dir))
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))
