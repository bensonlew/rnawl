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
from api_base import ApiBase
import gzip
from Bio import SeqIO


class SmallTarget(ApiBase):
    def __init__(self, bind_object):
        super(SmallTarget, self).__init__(bind_object)
        self.result_dir = ''
        self.result_file = {}
        self.trans_gene = {}
        self.trans_isgene = {}
        self.task_id = self.bind_object.sheet.id
        self.project_sn = self.bind_object.sheet.project_sn
        self.has_new = True
        self.anno_type = 'origin'
        self.species_name = ""
        self.target_smallrna = dict()
        self.tran2gene = dict()
        self.kegg_json = Config().SOFTWARE_DIR + "/database/KEGG/br08901.json"
        self.annot_type = "ref"
        self.target_name = ""
        self.version = "v1.1"
        #self._db_name = Config().MONGODB + '_ref_rna'

    def import_target(self, target, target_new=None, version="v1"):
        '''
        导入 靶基因信息
        '''
        self.version = version
        with open(target, 'rb') as target_f:
            target_f.readline()
            for line in target_f:
                if line.startswith("#"):
                    pass
                else:
                    mirna, target = line.strip().split("\t")[:2]
                    if target in self.target_smallrna:
                        self.target_smallrna[target].append(mirna)
                    else:
                        self.target_smallrna[target] = [mirna]
        if target_new:
            with open(target_new, 'rb') as target_f:
                target_f.readline()
                for line in target_f:
                    if line.startswith("#"):
                        pass
                    else:
                        mirna, target = line.strip().split("\t")[:2]
                        if target in self.target_smallrna:
                            self.target_smallrna[target].append(mirna)
                        else:
                            self.target_smallrna[target] = [mirna]

    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                raise Exception('assemble_id必须为ObjectId对象或其对应的字符串！')
        return object_id

    '''
    def get_trans2gene(self, trans2gene):
        if trans2gene:
            with open(trans2gene, 'rb') as f:
                lines = f.readlines()
                for line in lines:
                    cols = line.strip().split("\t")
                    self.tran2gene[cols[1]] = cols[0]
    '''

    def get_trans2gene(self, trans2gene_ref):
        '''
        根据靶基因长度判断哪一条转录本作为unigene的代表序列
        '''

        if trans2gene_ref:
            gene2length = dict()
            gene2trans = dict()
            with open(trans2gene_ref, 'rb') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    self.trans_gene[line[0]] = line[1]
                    if line[2] == "yes":
                        self.trans_isgene[line[0]] = True
                    else:
                        self.trans_isgene[line[0]] = False

        self.tran2gene = self.trans_gene

        self.bind_object.logger.info("读入基因转录本对应关系文件结束")

    def import_target_detail(self, target_file, params_dict, new_seq=None, known_seq=None, anno_type="origin", species_name=None, target_dir=None, version="v1"):
        self.version = version
        self.bind_object.logger.info("开始导入靶基因")
        task_id = self.task_id
        result_dir = os.path.dirname(target_file)
        print result_dir
        params_dict.update({
            "submit_location": "targe",
            "task_id": task_id,
            "task_type": 2
        })

        self.remove_table_by_main_record(main_table='target_mi', task_id=task_id, detail_table=['target_mi_detail', 'target_mi_stat'], detail_table_key='target_id')
        miset_id = self.import_diff_mi_id(known_seq, new_seq)
        targetset_id = self.import_diff_target_id(target_file)
        params_dict.update({
            "miset_id": str(miset_id),
            "targetset_id": str(targetset_id)
        })

        target_id, columns = self.add_target(params=params_dict, name=None, species_name=species_name, result_dir=result_dir)

        method_set = set()
        for method in [
            'm_miranda', 'm_targetscan', 'm_psrobot', 'm_targetfinder', 'm_rnahybrid',
            'l_miranda', 'l_targetscan', 'l_psrobot', 'l_targetfinder', 'l_rnahybrid',
            'c_miranda', 'c_targetscan', 'c_psrobot', 'c_targetfinder', 'c_rnahybrid']:
            if params_dict[method] == "yes":
                method_set.add(method)

        for target_sub_file in ['c_known', 'c_novel', 'l_known', 'l_novel', 'm_known', 'm_novel']:
            if target_sub_file.startswith("m_"):
                method_set_sub = [m.split("_")[1] for m in method_set if m.startswith("m_")]
            if target_sub_file.startswith("l_"):
                method_set_sub = [m.split("_")[1] for m in method_set if m.startswith("l_")]
            if target_sub_file.startswith("c_"):
                method_set_sub = [m.split("_")[1] for m in method_set if m.startswith("c_")]
            self.add_target_detail(target_id, target_file, target_sub_file, new_seq, known_seq, method_set = method_set_sub, target_dir=target_dir)

    def import_diff_mi_id(self, known_seq, new_seq):
        seq_list = list()
        category_list = list()
        kind_list = list()
        for seq in SeqIO.parse(known_seq, "fasta"):
            seq_list.append(seq.id)
            category_list.append("miRNA")
            kind_list.append("ref")
        for seq in SeqIO.parse(new_seq, "fasta"):
            seq_list.append(seq.id)
            category_list.append("miRNA")
            kind_list.append("new")

        source = "DE_miR_T_detail"
        level = "T"
        name = "diff_mirna_set"
        desc = "所有默认比较组获得的差异miRNA集合"
        diff_mirnaset_id = self.add_geneset(name, level, source, seq_list, category_list, kind_list, desc)
        return diff_mirnaset_id

    def import_diff_target_id(self, target_file):
        seq_list = list()
        category_list = list()
        kind_list = list()
        print target_file
        for target_sub_file in ['c_known', 'c_novel', 'l_known', 'l_novel', 'm_known', 'm_novel']:
            seq_file = os.path.join(target_file, target_sub_file, "target.fa")
            category = {'m': 'mRNA', 'c': 'circRNA', 'l': 'lncRNA'}.get(target_sub_file.split('_')[0])
            kind = {'known': 'ref', 'novel': 'new'}.get(target_sub_file.split('_')[1])
            if os.path.exists(seq_file):
                for seq in SeqIO.parse(seq_file, "fasta"):
                    # 去除序列可能的冗余
                    if seq.id in seq_list:
                        continue
                    seq_list.append(seq.id)
                    category_list.append(category)
                    kind_list.append(kind)

        source = "merge"
        level = "T"
        name = "diff_non_miRNA_set"
        desc = "所有默认比较组获得的差异ncRNA和mRNA集合"
        diff_target_id = self.add_geneset(name, level, source, seq_list, category_list, kind_list, desc)
        return diff_target_id

    def add_geneset(self, name, level, source, seq_list, category_list, kind_list, desc):
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
            'category_list': category_list,
            'kind_list': kind_list
        }
        self.add_geneset_detail(arg_dict, main_id)
        self.update_db_record('geneset', main_id, insert_dict={'main_id': main_id, 'status': 'end'})
        return main_id

    def add_geneset_detail(self, arg_dict, geneset_id):
        self.create_db_table('geneset_detail', [arg_dict], {'geneset_id': geneset_id})

    def import_target_detail_web(self, target_id, target_file, method_set, new_seq=None, known_seq=None, anno_type="latest", species_name=None, last_id_target=None, target_dir=None, version="v1"):
        self.bind_object.logger.info("开始导入靶基因")
        target_id = ObjectId(target_id)


        for target_sub_file in ['c_known', 'c_novel', 'l_known', 'l_novel', 'm_known', 'm_novel']:
            if target_sub_file.startswith("m_"):
                method_set_sub = [m.split("_")[1] for m in method_set if m.startswith("m_")]
            if target_sub_file.startswith("l_"):
                method_set_sub = [m.split("_")[1] for m in method_set if m.startswith("l_")]
            if target_sub_file.startswith("c_"):
                method_set_sub = [m.split("_")[1] for m in method_set if m.startswith("c_")]
            self.add_target_detail(target_id, target_file, target_sub_file, new_seq, known_seq, method_set = method_set_sub, target_dir=target_dir)

        if last_id_target:
            self.bind_object.logger.info("删除表格为 {}".format(last_id_target))
            self.remove_table_by_main_record(main_table='target_mi', _id=last_id_target, detail_table=['target_mi_detail', 'target_mi_stat'], detail_table_key='target_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")


        method_uniq = set([m.split("_")[1] for m in method_set])
        columns, columns_select = self.get_columns(method_uniq)

        print "columns is"
        print columns, columns_select
        self.db['target_mi'].update({"_id": target_id},
                                    {"$set": {"status": "end", "columns": columns, "columns_select": columns_select, "main_id": target_id}})


    def check_soft(self, params, soft):
        if 'm_' + soft in params and params['m_' + soft] == "yes":
            return True
        elif 'c_' + soft in params and params['c_' + soft] == "yes":
            return True
        elif 'l_' + soft in params and params['l_' + soft] == "yes":
            return True
        else:
            return False


    @report_check
    def add_target(self, params, name=None, species_name=None, result_dir=None, target_id=None):
        if target_id:
            if not isinstance(target_id, ObjectId):
                if isinstance(target_id, types.StringTypes):
                    target_id = ObjectId(target_id)
                else:
                    raise Exception('target_id必须为ObjectId对象或其对应的字符串！')
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        columns = []
        columns_select = []
        if self.check_soft(params, 'miranda'):
            columns.extend(['start_miranda', 'end_miranda', 'score_miranda', 'energy_miranda'])

            columns.extend(['paired_miranda'])
            columns_select.extend(['score_miranda', 'energy_miranda'])
        if self.check_soft(params, 'targetscan'):
            columns.extend(['utr_start_targetscan', 'utr_end_targetscan', 'msa_start_targetscan', 'msa_end_targetscan', 'sead_match_targetscan'])
            columns.extend(['pct', 'context_score'])
            columns_select.extend(['pct'])
        if self.check_soft(params, 'psrobot'):
            columns.extend(['start_psrobot', 'end_psrobot', 'score_psrobot'])
            columns.extend(['paired_psrobot'])
            columns_select.extend(['score_psrobot'])

        if self.check_soft(params, 'targetfinder'):
            columns.extend(['start_targetfinder', 'end_targetfinder', 'score_targetfinder'])
            columns.extend(['paired_targetfinder'])
            columns_select.extend(['score_targetfinder'])
        if self.check_soft(params, 'rnahybrid'):
            columns.extend(['start_rnahybrid', 'end_rnahybrid', 'energy_rnahybrid', 'pvalue_rnahybrid'])
            columns.extend(['paired_rnahybrid'])
            columns_select.extend(['pvalue_rnahybrid', 'energy_rnahybrid'])

        if target_id:
            self.db['target_mi'].update({"_id": target_id},
                    {"$set": {"columns": columns, "columns_select": columns_select, "main_id": target_id}})
            collection = self.db["target_mi"]
            result = collection.find_one({"_id": target_id})
            self.target_name = result['name']
            return target_id, columns

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        if name:
            self.target_name = name
        else:
            self.target_name = 'SmallTarget_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': self.target_name,
            'type': self.anno_type,
            'params': params,
            "version": "v1",
            'status': 'start',
            'desc': '注释统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'columns': columns,
            'columns_select': columns_select,
            'species_name': species_name,
            'result_dir': result_dir
            # 'seq_type': seq_type,
        }

        collection = self.db['target_mi']
        stat_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add target_mi!")
        return stat_id, columns


    @report_check
    def get_columns(self, methods):
        columns = []
        columns_select = []
        if 'miranda' in methods:
            columns.extend(['start_miranda', 'end_miranda', 'score_miranda', 'energy_miranda'])
            if self.version in ['v1.1']:
                columns.extend(['paired_miranda'])
            columns_select.extend(['score_miranda', 'energy_miranda'])
        if 'targetscan' in methods:
            columns.extend(['utr_start_targetscan', 'utr_end_targetscan', 'msa_start_targetscan', 'msa_end_targetscan', 'sead_match_targetscan'])
            if self.version in ['v1.1']:
                columns.extend(['pct', 'context_score'])
                columns_select.extend(['pct'])
        if 'psrobot' in methods:
            columns.extend(['start_psrobot', 'end_psrobot', 'score_psrobot'])
            if self.version in ['v1.1']:
                columns.extend(['paired_psrobot'])
            columns_select.extend(['score_psrobot'])

        if 'targetfinder' in methods:
            columns.extend(['start_targetfinder', 'end_targetfinder', 'score_targetfinder'])
            if self.version in ['v1.1']:
                columns.extend(['paired_targetfinder'])
            columns_select.extend(['score_targetfinder'])
        if 'rnahybrid' in methods:
            columns.extend(['start_rnahybrid', 'end_rnahybrid', 'energy_rnahybrid', 'pvalue_rnahybrid'])
            if self.version in ['v1.1']:
                columns.extend(['paired_rnahybrid'])
            columns_select.extend(['pvalue_rnahybrid', 'energy_rnahybrid'])

        return columns, columns_select

    def add_target_geneset(self, new_target_file, known_target_file, diff_summary):
        diff_small_dict = dict()
        diff_target_dict = dict()
        groups = list()

        with open(diff_summary, 'r') as diff_f:
            header = diff_f.readline()
            groups = header.strip("\n").split("\t")[1:-1]
            for group in groups:
                diff_small_dict[group] = []
                diff_target_dict[group] = set()
            nums= diff_f.readline()
            for line in diff_f:
                small_rna = line.strip("\n").split("\t")[0]
                stats = line.strip("\n").split("\t")[1:-1]
                for n,stat in enumerate(stats):
                    if stat == "yes":
                        diff_small_dict[groups[n]].append(small_rna)
                    else:
                        pass

        with open(known_target_file, 'r') as known_target:
            known_target.readline()
            for line in known_target:
                cols = line.strip().split("\t")
                for group in groups:
                    if cols[0] in diff_small_dict[group]:
                        diff_target_dict[group].add(cols[2])

        with open(new_target_file, 'r') as new_target:
            head_line = new_target.readline()
            for line in new_target:
                cols = line.strip().split("\t")
                for group in groups:
                    if cols[0] in diff_small_dict[group]:
                        diff_target_dict[group].add(cols[2])

        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        for group in groups:
            geneset_main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=group + "_targets",
                type="G",
                desc='差异smallrna靶基因集',
                group_id="",
                gene_length=len(diff_target_dict[group]),
                is_use=0
            )

            geneset_detail_info = [{"seq_list": list(diff_target_dict[group])}]
            if len(diff_target_dict[group]) != 0:
                self.add_set(geneset_main_info, geneset_detail_info)

    @report_check
    def add_set(self, main_info, detail_info):
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        main_info.update(dict(status="start", created_ts=created_ts))
        main_id = self.create_db_table('sg_geneset', [main_info])
        self.create_db_table('sg_geneset_detail', detail_info, tag_dict={"geneset_id": main_id})
        # self.update_db_record('sg_geneset', main_id, status="end", is_use=0, main_id=main_id)
        task_id = main_info['task_id']
        record_dict = {"_id": main_id, "task_id": task_id}
        self.update_db_record('sg_geneset', main_id, query_dict=record_dict, status="end", is_use=0, main_id=main_id, params=task_id)
        return main_id


    def get_target_detail(self, target_file, soft):
        target_dict = dict()
        region = False
        if soft == "rnahybrid":
            with gzip.open(target_file, 'rb') as f:
                for line in f:
                    if line.startswith("target: "):
                        tar = line.strip().split()[1]
                    if line.startswith("miRNA : "):
                        mi = line.strip().split()[-1]
                    if line.startswith("target "):
                        region = True

                    if region:
                        if mi + "|" + tar in target_dict:
                            target_dict[mi + "|" + tar].append(line.strip("\n"))
                        else:
                            target_dict[mi + "|" + tar] = [line.strip("\n")]
                    if line.startswith("miRNA  "):
                        region = False
        if soft == "miranda":
             with gzip.open(target_file, 'rb') as f:
                for line in f:
                    if line.startswith("Performing Scan:"):
                        tar = line.strip().split()[4]
                        mi = line.strip().split()[2]

                    if line.startswith("   Query:"):
                        region = True

                    if region:
                        if mi + "|" + tar in target_dict:
                            target_dict[mi + "|" + tar].append(line.strip("\n"))
                        else:
                            target_dict[mi + "|" + tar] = [line.strip("\n")]
                    if line.startswith("   Ref:"):
                        region = False
        if soft == "targetfinder":
            with gzip.open(target_file, 'rb') as f:
                num = 0
                for line in f:
                    num += 1
                    if num % 6 == 1:
                        tar = line.strip().split()[2]
                        mi = line.strip().split()[0]


                    if num % 6 in [3,4,5]:
                        if mi + "|" + tar in target_dict:
                            target_dict[mi + "|" + tar].append(line.strip("\n"))
                        else:
                            target_dict[mi + "|" + tar] = [line.strip("\n")]

        if soft == "psrobot":
            with gzip.open(target_file, 'rb') as f:
                num = 0
                for line in f:
                    num += 1
                    if num % 7 == 1:
                        tar = line.strip().split()[3]
                        mi = line.strip().split()[0].lstrip(">")

                    if num % 7 in [3,4,5]:
                        if mi + "|" + tar in target_dict:
                            target_dict[mi + "|" + tar].append(line.strip("\n"))
                        else:
                            target_dict[mi + "|" + tar] = [line.strip("\n")]

        return target_dict


    def change_align_list2html(self, data, soft):
        if soft == "miranda":
            a = data[0][16:-3]
            b = data[1][16:]
            b = b.replace(" ", "*")
            c = data[2][16:-3]


            b_list = []
            for i in range(0, len(c)):
                if a[i] + c[i] in ['AU', 'UA', 'GC', 'CG' , 'at', 'ta', 'au', 'ua', 'gc', 'cg']:
                    b_list.append("|")
                elif a[i] + c[i] in ['gt', 'gu', 'tg', 'ug']:
                    b_list.append(":")
                else:
                    b_list.append(b[i])
            b = "".join(b_list)
            pair_list = [
                "miRNA:  3' {} 5'".format(a),
                "           {}   ".format(b),
                "target: 5' {} 3'".format(c),
            ]

            return "<table class='targe_detail_paired'><tr><td>miRNA:</td><td>3'</td><td >{}</td><td>5'</td></tr><tr><td></td><td></td><td>{}</td></tr><tr><td>target:</td><td>5'</td><td>{}</td><td>3'</td></tr></table>".format(a, b, c), pair_list

        if soft == "rnahybrid":
            line1 = data[0][16:-3]
            line2 = data[1][16:-3]
            line3 = data[2][16:-3]
            line4 = data[3][16:-3]
            a_list = list()
            b_list = list()
            c_list = list()
            for i in range(0, len(line1)):
                if line1[i] != " ":
                    a_list.append(line1[i])
                elif line2[i] != " ":
                    a_list.append(line2[i])
                else:
                    a_list.append("-")

            for i in range(0, len(line1)):
                if line3[i] != " ":
                    c_list.append(line3[i])
                elif line4[i] != " ":
                    c_list.append(line4[i])
                else:
                    c_list.append("-")

            for i in range(0, len(a_list)):
                if a_list[i] + c_list[i] in ['AU', 'UA', 'GC', 'CG']:
                    b_list.append("|")
                elif a_list[i] + c_list[i] in ['GU', 'UG']:
                    b_list.append(":")
                else:
                    b_list.append("*")

            a = "".join(a_list)
            b = "".join(b_list)
            c = "".join(c_list)

            pair_list = [
                "target: 5' {} 3'".format(a),
                "           {}   ".format(b),
                "miRNA:  3' {} 5'".format(c)
            ]
            return "<table class='targe_detail_paired'><tr><td>target:</td><td>5'</td><td>{}</td><td>3'</td></tr><tr><td></td><td></td><td>{}</td></tr><tr><td>miRNA:</td><td>3'</td><td >{}</td><td>5'</td></tr></table>".format(a, b, c), pair_list

        if soft == "psrobot":
            a = data[0][18:].split()[0]
            b = data[1][18:]
            c = data[2][18:].split()[0]
            pair_list = [
                "target: 5' {} 3'".format(a),
                "           {}   ".format(b),
                "miRNA:  3' {} 5'".format(c)
            ]
            return "<table class='targe_detail_paired'><tr><td>target:</td><td>5'</td><td>{}</td><td>3'</td></tr><tr><td></td><td></td><td>{}</td></tr><tr><td>miRNA:</td><td>3'</td><td >{}</td><td>5'</td></tr></table>".format(a, b, c), pair_list
        if soft == "targetfinder":
            a = data[0][11:-3]
            b = data[1][11:-3].replace(".", ":").replace(" ", "*")
            c = data[2][11:-3]
            pair_list = [
                "miRNA:  3' {} 5'".format(a),
                "           {}   ".format(b),
                "target: 5' {} 3'".format(c),
            ]
            return "<table class='targe_detail_paired'><tr><td>miRNA:</td><td>3'</td><td >{}</td><td>5'</td></tr><tr><td></td><td></td><td>{}</td></tr><tr><td>target:</td><td>5'</td><td>{}</td><td>3'</td></tr></table>".format(a, b, c), pair_list

    @report_check
    def add_target_detail(self, target_id, target_file, target_sub_file, new_seq=None, known_seq=None, method_set = set(), target_dir=None):
        if not isinstance(target_id, ObjectId):
            if isinstance(target_id, types.StringTypes):
                target_id = ObjectId(target_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')

        columns, columns_select = self.get_columns(method_set)
        columns = [x for x in columns if not x.startswith("paired")]

        known_target_file = os.path.join(target_file, target_sub_file, "known_target.xls")
        new_target_file = os.path.join(target_file, target_sub_file, "novol_target.xls")
        target_dir = os.path.join(target_file, target_sub_file)

        category = "mRNA"
        if target_sub_file.startswith("l"):
            category = "lncRNA"
        if target_sub_file.startswith("c"):
            category = "circRNA"

        target_kind = "ref"
        if target_sub_file.endswith("novel"):
            target_kind = "new"

        if not os.path.exists(new_target_file):
            # circ 可能不预测
            return
            # raise Exception('{}所指定的路径不存在，请检查！'.format(new_target_file))
        if not os.path.exists(known_target_file):
            return
            # raise Exception('{}所指定的路径不存在，请检查！'.format(known_target_file))

        # header = ['small_rna', 'target'] + columns + ['mir_tar_base']
        data_list = []
        known_smallrna_list = []
        novol_smallrna_list = []
        known_target_list = []
        novol_target_list = []
        header = ['small_rna', 'target', 'gene', 'name'] + columns + ['mirtarbase']

        paired_known = dict()
        paired_novol = dict()
        for soft in ["miranda", "rnahybrid", "targetscan", "targetfinder", "psrobot"]:
            if os.path.exists(target_dir + "/known_" + soft + "_detail.txt.gz"):
                paired_known[soft] =  self.get_target_detail(target_dir + "/known_" + soft + "_detail.txt.gz", soft)
                paired_novol[soft] =  self.get_target_detail(target_dir + "/novol_" + soft + "_detail.txt.gz", soft)
            else:
                pass

        with open(known_target_file, 'r') as known_target:
            known_target.readline()
            for line in known_target:
                cols = line.strip("\n").split("\t")
                cols = cols[:4] + cols[5:]
                for n in range(4, len(cols)):
                    try:
                        cols[n] = float(cols[n])
                    except Exception:
                        pass
                data = zip(header, cols)
                data.append(('target_id', target_id))
                data.append(('mi_kind', 'ref'))
                data.append(('target_kind', target_kind))
                data.append(('category', category))
                data = SON(data)
                for soft, pd_dict in paired_known.items():
                    if cols[0] + '|' + cols[1] in pd_dict:
                        data["paired_" + soft], data["pairlist_" + soft]  = self.change_align_list2html(pd_dict[cols[0] + '|' + cols[1]], soft)
                    else:
                        pass
                data_list.append(data)
                known_smallrna_list.append(cols[0])
                known_target_list.append(cols[2])


        with open(new_target_file, 'r') as new_target:
            head_line = new_target.readline()
            for line in new_target:
                cols = line.strip("\n").split("\t")
                cols = cols[:4] + cols[5:]
                for n in range(4, len(cols)):
                    try:
                        cols[n] = float(cols[n])
                    except Exception:
                        pass
                data = zip(header, cols)
                data.append(('target_id', target_id))
                data.append(('mi_kind', 'new'))
                data.append(('target_kind', target_kind))
                data.append(('category', category))
                data = SON(data)
                for soft, pd_dict in paired_novol.items():
                    if cols[0] + '|' + cols[1] in pd_dict:
                        data["paired_" + soft], data["pairlist_" + soft] = self.change_align_list2html(pd_dict[cols[0] + '|' + cols[1]], soft)
                    else:
                        pass
                data_list.append(data)
                novol_smallrna_list.append(cols[0])
                novol_target_list.append(cols[2])

        all_known_small_list = []
        all_novol_small_list = []

        if known_seq:
            for seq in SeqIO.parse(known_seq, "fasta"):
                all_known_small_list.append(seq.id)
        if new_seq:
            for seq in SeqIO.parse(new_seq, "fasta"):
                all_novol_small_list.append(seq.id)

        data_stat_list = [
            SON([('type', 'Known miRNA'),
                 ('target_id', target_id),
                 ('all_mirna', len(set(all_known_small_list))),
                 ('mirna_target', len(set(known_smallrna_list))),
                 ('target', len(set(known_target_list))),
            ]),
            SON([('type', 'Novel miRNA'),
                 ('target_id', target_id),
                 ('all_mirna', len(set(all_novol_small_list))),
                 ('mirna_target', len(set(novol_smallrna_list))),
                 ('target', len(set(novol_target_list))),
            ]),
            SON([('type', 'Total'),
                 ('target_id', target_id),
                 ('all_mirna', len(set(all_known_small_list + all_novol_small_list))),
                 ('mirna_target', len(set(known_smallrna_list + novol_smallrna_list))),
                 ('target', len(set(known_target_list + novol_target_list))),
            ])
        ]

        try:
            self.create_db_table('target_mi_detail', data_list)
            self.create_db_table('target_mi_stat', data_stat_list)
            self.db['target_mi'].update({"_id": target_id},
                                           {"$set": {"status": "end", "main_id": target_id}})
        except Exception as e:
            self.bind_object.set_error("导入靶基因预测统计信息出错!")
        else:
            self.bind_object.logger.info("导入靶基因预测统计信息成功!")


    def run(self, target_file, params_dict, taxon='Animals', exp_level='transcript', version="v1.0"):
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
        self.bind_object.logger.info("开始到表情数据路径为 {}".format(target_file))
        self.import_target(target_file)
        # self.set_result_dir(annotation_mudule_dir)
        # self.get_trans2gene(trans2gene)

        task_id = self.task_id



    def run_web(self, target_file, target_file2, annotation_mudule_dir, trans2gene, params_dict, task_id, stat_id, last_id, taxon='Animals', exp_level='transcript'):
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




    def run_webroot(self, result_dir, trans2gene, trans2gene_ref, params_dict, task_id, stat_id, last_id, taxonomy, exp_level):
        """
        用于注释重运行导表
        result_dir: 新序列注释的结果文件夹
        trans2gene: 转录本基因对应关系文件
        params_dict: 参数
        stat_id: 统计结果主表ID
        last_id: 上次重运行ID
        taxonomy: 物种分类
        """
        self.bind_object.logger.info("开始导表webroot数据路径为 {}".format(result_dir))
        self.set_result_dir(result_dir)
        self.get_trans2gene(trans2gene, trans2gene_ref)
        # self.get_trans2gene(trans2gene)
        self.bind_object.logger.info("开始导表task_id为 {}".format(self.task_id))

        stat_id = ObjectId(stat_id)
        if last_id:
            last_id = ObjectId(last_id)
        else:
            pass
        self.task_id = task_id
        #task_id = self.task_id


class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """
    def test_mongo(test):
        # whole_transcriptome/small_target
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet


        data = {
            "id": "whole_transcriptome",
            "project_sn": "whole_transcriptome",
            "type": "workflow",
            "name": 'whole_transcriptome.whole_transcriptome_test_api',
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wsheet)

        target = "/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/TargetMirna"

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("whole_transcriptome.small_target")

        wf.test_api.species_name = "Homo_sapiens"

        params_target = {
            "miranda": "yes",
            "targetscan": "yes",
            "rnahybrid": "yes",
            "min_support": 2,
        }

        for method in [
            'm_miranda', 'm_targetscan', 'm_psrobot', 'm_targetfinder', 'm_rnahybrid',
            'l_miranda', 'l_targetscan', 'l_psrobot', 'l_targetfinder', 'l_rnahybrid',
             'c_miranda', 'c_targetscan', 'c_psrobot', 'c_targetfinder', 'c_rnahybrid']:
            params_target.update({method: "yes"})

        for par in ['m_miranda_score', 'm_miranda_energy', 'm_miranda_strict',
                    'm_rnahybird_num', 'm_rnahybird_energy', 'm_rnahybird_pvalue',
                    'm_ps_robot_score', 'm_targetfinder_score',
                    'l_miranda_score', 'l_miranda_energy', 'l_miranda_strict',
                    'l_rnahybird_num', 'l_rnahybird_energy', 'l_rnahybird_pvalue',
                    'l_ps_robot_score', 'l_targetfinder_score'
                    'c_miranda_score', 'c_miranda_energy', 'c_miranda_strict',
                    'c_rnahybird_num', 'c_rnahybird_energy', 'c_rnahybird_pvalue',
                    'c_ps_robot_score', 'c_targetfinder_score'
        ]:
            if hasattr(data, par):
                params_target.update({
                    par: "1"
                })

        target = "/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/TargetMirna"
        new_seq = '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer1/output/srna/novel_mirna/novel_mature_seq.fa'
        known_seq =  '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer1/output/srna/known_mirna/mature.fa'
        target_file = "/mnt/ilustre/users/sanger-dev/workspace/20191016/Single_whole_7047_9776/TargetMirna/output"
        wf.test_api.import_target_detail(target_file, params_target, new_seq, known_seq, anno_type="origin", species_name="Homo_sapiens")


if __name__ == '__main__':
    unittest.main()
