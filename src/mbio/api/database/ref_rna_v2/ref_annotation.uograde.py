# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import re
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
import json
import unittest
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from mbio.api.database.ref_rna_v2.api_base import ApiBase
import sys


class RefAnnotation(ApiBase):
    def __init__(self, bind_object):

        super(RefAnnotation, self).__init__(bind_object)
        self.kegg_json = Config().SOFTWARE_DIR + "/database/KEGG/br08901.json"
        self.trans_gene = dict()


    def run(self, task_id, go_dir):
        """
        """

        conn = self.db["sg_annotation_go"]
        main_id = conn.find_one({'task_id': task_id,  'type': "origin"})['main_id']

        main_id = ObjectId(main_id)
        print main_id

        conn_detail = self.db["sg_annotation_go_detail"]
        query_dict = {"seq_type": "all", "anno_type": "G", "go_id": main_id}

        result = conn_detail.find_one(query_dict)
        print "one_result is {}".format(result)
        if result:
            conn_detail.delete_many(query_dict)
        else:
            pass
            # print 'fail to find {} by {}'.format(table_name, kwargs)

        go_id = main_id
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=2, r_go_path=go_dir + '/go_lev2_gene.stat.xls')
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=3, r_go_path=go_dir + '/go_lev3_gene.stat.xls')
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=4, r_go_path=go_dir + '/go_lev4_gene.stat.xls')


        result = conn_detail.find_one(query_dict)
        print result

    def run2(self, task_id, go_dir):
        """
        """

        conn = self.db["sg_annotation_go"]
        main_id = conn.find_one({'task_id': task_id,  'type': "origin"})['main_id']

        main_id = ObjectId(main_id)
        print main_id

        conn_detail = self.db["sg_annotation_go_detail"]
        query_dict = {"seq_type": "new", "anno_type": "G", "go_id": main_id}

        result = conn_detail.find_one(query_dict)

        print "one_result is {}".format(result)

        if result:
            conn_detail.delete_many(query_dict)
        else:
            pass
            # print 'fail to find {} by {}'.format(table_name, kwargs)

        merge_annot_path = go_dir
        go_id = main_id

        seq_type = "new"
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=merge_annot_path + "/newannot_class/go/go_lev2_gene.stat.xls")
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=merge_annot_path + "/newannot_class/go/go_lev3_gene.stat.xls")
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=merge_annot_path + "/newannot_class/go/go_lev4_gene.stat.xls")


        query_dict = {"seq_type": "ref", "anno_type": "G", "go_id": main_id}

        result = conn_detail.find_one(query_dict)

        print "one_result is {}".format(result)

        if result:
            conn_detail.delete_many(query_dict)
        else:
            pass

        seq_type = "ref"
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=merge_annot_path + "/refannot_class/go/go_lev2_gene.stat.xls")
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path= merge_annot_path + "/refannot_class/go/go_lev3_gene.stat.xls")
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=merge_annot_path + "/refannot_class/go/go_lev4_gene.stat.xls")

        result = conn_detail.find_one(query_dict)
        print result



    def run_pfam(self, trans2gene, trans2gene_ref, pfam_id, pfam_path, pfam_path_ref, query_id, query_path, query_path_ref, gene_pfam_path, gene_pfam_pathref,  gene_pfam_pathall , ta):
        '''
        pass
        '''
        pfam_id = ObjectId(pfam_id)
        query_id = ObjectId(query_id)


        self.get_trans2gene(trans2gene, trans2gene_ref)

        self.result_file = dict()
        self.result_file['pfam_path'] = pfam_path
        self.result_file['pfam_path' + "ref"] = pfam_path_ref
        self.result_file['pfam_path' + "all"] = ta

        self.result_file['gene_pfam_path'] = gene_pfam_path
        self.result_file['gene_pfam_pathref'] = gene_pfam_pathref
        self.result_file['gene_pfam_pathall'] = gene_pfam_pathall
        exp_level = "transcript"
        self.has_new = True

        print "删除 {}".format('sg_annotation_pfam_detail')
        self.remove_db_record('sg_annotation_pfam_detail',
                              query_dict={'pfam_id': pfam_id})

        print "删除 {}".format('sg_annotation_pfam_bar')
        self.remove_db_record('sg_annotation_pfam_bar',
                              query_dict={'pfam_id': pfam_id})

        if self.has_new:
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T", exp_level=exp_level)
        self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'  + "ref"], seq_type="ref", anno_type="T", exp_level=exp_level)

        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path'], seq_type="new", anno_type="G")
        if exp_level.lower() == "transcript":
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path' + "ref"], seq_type="ref", anno_type="T")
        self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path' + "ref"], seq_type="ref", anno_type="G")
        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path' + "all"], seq_type="all", anno_type="T")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path' + "all"], seq_type="all", anno_type="G")




        self.remove_db_record('sg_annotation_query_detail',
                              query_dict={'query_id': query_id})


        self.result_file['query_path'] = query_path
        self.result_file['query_path' + "ref"] = query_path_ref

        if self.has_new:
            self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path'], seq_type = "new", anno_type="T", exp_level=exp_level)
        self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path' + "ref"], seq_type = "ref", anno_type="T", exp_level=exp_level)
        self.update_db_record('sg_annotation_query',query_id, status="end", main_id=query_id)


    def get_kegg_map2class(self):
        '''
        返回kegg mapid 与 classI classII 对应关系字典
        key: map00010
        value: (classi_name, classii_name)
        '''
        map2class = dict()
        with open(self.kegg_json, "rb") as f:
            root = json.load(f)
        classI = root['children']

        for class1 in classI:
            class1_name = class1['name']
            for class2 in class1['children']:
                class2_name = class2['name']
                paths = class2['children']
                for path in paths:
                    path_name = "map" + str(path['name']).split(" ")[0]
                    map2class[path_name] = (class1_name, class2_name)
        return map2class


    def add_annotation_pfam_detail(self, pfam_id, pfam_path, seq_type, anno_type, exp_level):
        """
        pfam_path: pfam_domain
        """
        if not isinstance(pfam_id, ObjectId):
            if isinstance(pfam_id, types.StringTypes):
                pfam_id = ObjectId(pfam_id)
            else:
                self.bind_object.set_error('pfam_id必须为ObjectId对象或其对应的字符串！', code="537018149")
        if not os.path.exists(pfam_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pfam_path), code="537018150")
        data_list = []
        with open(pfam_path, "r") as f:
            lines = f.readlines()
            last_seq_id = ''
            last_pfam_id = ''
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[2] != last_seq_id or line[3] != last_pfam_id:
                    # 过滤不在参考范围的转录本注释
                    if exp_level.lower() == "gene" or line[0] not in self.trans_isgene:
                        continue
                    data = [
                        ('pfam_id', pfam_id),
                        ('seq_type', seq_type),
                        # ('anno_type', anno_type),
                        ('transcript_id', line[0]),
                        ('pfam', line[2]),
                        ('domain', line[3]),
                        ('description', line[4]),
                        ('protein_id', line[1]),
                        ('e_value', float(line[9])),
                        ('length', int(line[6])-int(line[5])),
                        ('protein_start', int(line[5])),
                        ('protein_end', int(line[6])),
                        ('pfam_start', int(line[7])),
                        ('pfam_end', int(line[8])),
                    ]
                    if  anno_type == 'T':
                        try:
                            data.append(('gene_id', self.trans_gene[line[0]]))
                            data.append(('is_gene', self.trans_isgene[line[0]]))
                        except:
                            data.append(('gene_id', line[0]))
                            data.append(('is_gene', True))
                        data = SON(data)
                        data_list.append(data)
                    else:
                        data = SON(data)
                        data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_pfam_detail']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入pfam注释信息:%s失败！" , variables=( pfam_path), code="537018151")
            else:
                print  "导入pfam注释信息:成功"



    def add_annotation_pfam_bar(self, pfam_id, pfam_path, seq_type, anno_type):
        pfam = []
        domain = {}
        with open(pfam_path, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[3] not in pfam:
                    pfam.append(line[3])
                    domain[line[3]] = 1
                else:
                    domain[line[3]] += 1
        if not isinstance(pfam_id, ObjectId):
            if isinstance(pfam_id, types.StringTypes):
                pfam_id = ObjectId(pfam_id)
            else:
                self.bind_object.set_error('pfam_id必须为ObjectId对象或其对应的字符串！', code="537018152")
        if not os.path.exists(pfam_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pfam_path), code="537018153")
        data_list = []
        dom = zip(domain.values(), domain.keys())
        dom_sort = sorted(dom, reverse=True)
        for num,dom in dom_sort:
            data = [
                ('pfam_id', pfam_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('domain', dom),
                ('num', num)
            ]
            data = SON(data)
            data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_pfam_bar']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入pfam注释信息:%s失败！" , variables=( pfam_path), code="537018154")
            else:
                pass
                # self.bind_object.logger.info("导入pfam注释信息:%s成功！" % pfam_path)


    def add_annotation_query_denovo_detail(self, query_id, query_path, seq_type, anno_type, exp_level):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="537018203")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="537018204")
        data_list = []
        map2class = self.get_kegg_map2class()
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                is_gene = True
                if self.trans_isgene.has_key(line[1]) and self.trans_isgene[line[1]] == False:
                    is_gene = False
                if exp_level.lower() == "gene" and is_gene == False:
                    continue
                data = [
                    ('query_id', query_id),
                    #('anno_type', anno_type),
                    ('seq_type', seq_type),
                    ('transcript_id', line[1]),
                    ('gene_id', line[0]),
                    ('is_gene', is_gene),
                    ('gene_name', line[3]),
                ]
                try:
                    data.append(('length', line[4]))
                except:
                    data.append(('length', None))
                try:
                    data.append(('description', line[5]))
                except:
                    data.append(('description', None))
                try:
                    data.append(('cog', line[6]))
                    data.append(('cog_description', line[7]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                # try:
                #     data.append(('nog', line[4]))
                #     data.append(('nog_description', line[6]))
                # except:
                #     data.append(('nog', None))
                #     data.append(('nog_description', None))
                try:
                    data.append(('ko_id', line[8]))
                except:
                    data.append(('ko_id', None))
                try:
                    data.append(('ko_name', line[9]))
                except:
                    data.append(('ko_name', None))

                try:
                    data.append(('pathways', line[10]))
                    paths = [pathway.split("(")[0] for pathway in line[10].split("; ")]
                    if line[10] != "":
                        pathway_class1 = [map2class[path][0] for path in paths if path in map2class]
                        pathway_class2 = [map2class[path][1] for path in paths if path in map2class]
                    data.append(('pathways_class1', ";".join(pathway_class1)))
                    data.append(('pathways_class2', ";".join(pathway_class2)))

                except:
                    data.append(('pathways', None))
                try:
                    data.append(('pfam', line[11]))
                except:
                    data.append(('pfam', None))
                try:
                    data.append(('go', line[12]))
                except:
                    data.append(('go', None))
                try:
                    data.append(('nr', line[13]))
                except:
                    data.append(('nr', None))
                try:
                    data.append(('swissprot', line[14]))
                except:
                    data.append(('swissprot', None))
                try:
                    data.append(('enterz', line[15]))
                except:
                    data.append(('enterz', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入转录本注释统计信息：%s出错!" , variables=(query_path), code="537018205")
        else:
            print "导入转录本注释统计信息：成功"

    def add_annotation_go_all(self, go_id, seq_type, anno_type, level, r_go_path):
        """
        r_go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls(ref)
        n_go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls(new)
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                pass
                # self.bind_object.set_error('go_id必须为ObjectId对象或其对应的字符串！', code="537018177")
        if not os.path.exists(r_go_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(r_go_path), code="537018178")
        data_list1, data_list2, query_ids = list(), list(), list()
        funlist, termlist = {}, {}
        with open(r_go_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                fun = line[0] + "|||" + line[1] + "|||" + line[2]
                term = line[0] + "|||" + line[1]
                if level == 3:
                    fun += "|||" + line[3] + "|||" + line[4]
                    term = line[0] + "|||" + line[3]
                if level == 4:
                    fun += "|||" + line[3] + "|||" + line[4]
                    fun += "|||" + line[5] + "|||" + line[6]
                    term = line[0] + "|||" + line[5]
                funlist[fun] = line[-1].split(";")
                if term not in termlist:
                    termlist[term] = set(line[-1].split(";"))
                else:
                    for q in line[-1].split(";"):
                        if q not in termlist[term]:
                            termlist[term].add(q)
                query_ids.extend(line[-1].split(";"))
        query_ids = list(set(query_ids))
        for term in sorted(termlist):
            terms = term.split("|||")
            seq_list = termlist[term]
            percent = float(len(seq_list)) / len(query_ids)
            data = [
                ('go_id', go_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('level', level),
                ('term_type', terms[0]),
                ('go_term', terms[1]),
                ('seq_number', len(seq_list)),
                ('percent', round(percent, 4))
                #('seq_list', ";".join(seq_list))
            ]
            data = SON(data)
            data_list1.append(data)
        if data_list1:
            try:
                collection = self.db['sg_annotation_go_graph']
                collection.insert_many(data_list1)
            except Exception, e:
                # self.bind_object.set_error("导入go注释画图all信息出错：%s" % (r_go_path))
                print "导入go注释画图all出错：%s" % (r_go_path)
            else:
                pass
                # self.bind_object.logger.info("导入go注释画图all信息成功：%s" % (r_go_path))
        for fun in funlist:
            terms = fun.split("|||")
            data = [
                ('go_id', go_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('level', level),
                ('goterm', terms[0]),
                ('goterm_2', terms[1]),
                ('goid_2', terms[2])
            ]
            if level >= 3:
                data.append(('goterm_3', terms[3]))
                data.append(('goid_3', terms[4]))
            if level == 4:
                data.append(('goterm_4', terms[5]))
                data.append(('goid_4', terms[6]))
            seq_list = funlist[fun]
            percent = float(len(funlist[fun])) / len(query_ids)
            data.append(('seq_number', len(seq_list)))
            data.append(('percent', round(percent, 4)))
            if level == 2:
                data.append(('seq_list', ";".join(seq_list)))
            data = SON(data)
            data_list2.append(data)
        if data_list2:
            try:
                collection = self.db['sg_annotation_go_detail']
                collection.insert_many(data_list2)
            except Exception, e:
                pass
                # self.bind_object.set_error("导入go注释all出错：%s" , variables=(r_go_path), code="537018179")
            else:
                pass
                # self.bind_object.logger.info("导入go注释all信息成功：%s" % (r_go_path))

    def add_annotation_go_detail(self, go_id, seq_type, anno_type, level, go_path):
        """
        go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id必须为ObjectId对象或其对应的字符串！', code="537018165")
        if not os.path.exists(go_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(go_path), code="537018166")
        data_list = list()
        with open(go_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('go_id', go_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('level', level),
                    ('goterm', line[0]),
                    ('goterm_2', line[1]),
                    ('goid_2', line[2]),
                    ('seq_number', int(line[-3])),
                    ('percent', round(float(line[-2]), 4)),
                    #('seq_list', line[-1])
                ]
                if level == 2:
                    data.append(('seq_list', line[-1]))
                if level >= 3:
                    data.append(('goterm_3', line[3]))
                    data.append(('goid_3', line[4]))
                if level == 4:
                    data.append(('goterm_4', line[5]))
                    data.append(('goid_4', line[6]))
                data = SON(data)
                data_list.append(data)
            if data_list:
                try:
                    collection = self.db['sg_annotation_go_detail']
                    collection.insert_many(data_list)
                except Exception, e:
                    self.bind_object.set_error("导入go注释信息：%s出错!" , variables=(go_path), code="537018167")
                else:
                    print "导入go注释信息：%s成功!" % (go_path)


    def get_trans2gene(self, trans2gene, trans2gene_ref=None):
        self.new_gene_set = set()
        self.known_gene_set = set()

        self.new_trans_set = set()
        self.known_trans_set = set()
        self.trans_isgene = dict()
        if trans2gene and os.path.exists(trans2gene):
            pass
        elif trans2gene_ref and os.path.exists(trans2gene_ref):
            pass
        else:
            pass
        if trans2gene:
            with open(trans2gene, 'rb') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    self.trans_gene[line[0]] = line[1]
                    self.new_trans_set.add(line[0])
                    if line[2] == "yes":
                        self.trans_isgene[line[0]] = True
                        self.new_gene_set.add(line[1])
                    else:
                        self.trans_isgene[line[0]] = False
        if trans2gene_ref:
            with open(trans2gene_ref, 'rb') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    self.trans_gene[line[0]] = line[1]
                    self.known_trans_set.add(line[0])
                    if line[2] == "yes":
                        self.known_gene_set.add(line[1])
                        self.trans_isgene[line[0]] = True
                    else:
                        self.trans_isgene[line[0]] = False



if __name__ == '__main__':
    trans2gene = sys.argv[1]
    trans2gene_ref = sys.argv[2]
    pfam_id = sys.argv[3]
    pfam_path = sys.argv[4]
    pfam_path_ref = sys.argv[5]
    query_id = sys.argv[6]
    query_path = sys.argv[7]
    query_path_ref = sys.argv[8]

    gn = sys.argv[9]
    gr = sys.argv[10]
    ga = sys.argv[11]
    ta = sys.argv[12]

    task_id = sys.argv[1]
    go_dir = sys.argv[2]
    ref_api = RefAnnotation(None)

    ref_api.run_pfam(trans2gene, trans2gene_ref, pfam_id, pfam_path, pfam_path_ref, query_id, query_path, query_path_ref, gn, gr, ga, ta)
    # ref_api.run2(task_id, go_dir)
