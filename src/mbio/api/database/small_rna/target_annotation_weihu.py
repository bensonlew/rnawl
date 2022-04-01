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
from mbio.api.database.small_rna.api_base import ApiBase
import gzip
from Bio import SeqIO


class TargetAnnotationWeihu(ApiBase):
    def __init__(self, bind_object):
        super(TargetAnnotationWeihu, self).__init__(bind_object)
        self.result_dir = ''
        self.result_file = {}
        self.trans_gene = {}
        self.trans_isgene = {}
        self.task_id = self.bind_object.sheet.id
        self.has_new = True
        self.anno_type = 'origin'
        self.species_name = ""
        self.target_smallrna = dict()
        self.tran2gene = dict()
        self.kegg_json = Config().SOFTWARE_DIR + "/database/KEGG/br08901.json"
        self.annot_type = "ref"
        self.target_name = ""
        self.version = "v1"
        #self._db_name = Config().MONGODB + '_ref_rna'


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


    def set_result_dir(self, annot_path):
        '''
        根据注释模块结果导入结果路径
        '''
        annot_path =  annot_path + "/"
        self.result_dir = annot_path
        self.bind_object.logger.info("导入**** {}".format(annot_path))
        self.result_file['ref_stat_path'] = os.path.join(annot_path, "anno_stat/all_annotation_statistics.xls")
        self.result_file['ref_venn_path'] = os.path.join(annot_path, "anno_stat/venn")

        ## 判断注释来源于哪个版本
        if os.path.exists(os.path.join(annot_path, "anno_stat/all_anno_detail.xls")):
            self.result_file['query_path' + "ref"] = os.path.join(annot_path, "anno_stat/all_anno_detail.xls")
            self.result_file['gos_path' + '_ref'] = annot_path + "go/query_gos.list"
            self.result_file['gene_gos_path' + '_ref'] = annot_path + "anno_stat/go_stat/gene_gos.list"

            self.result_file['table_path' + '_ref'] = annot_path + "kegg/kegg_table.xls"
            self.result_file['gene_table_path' + '_ref'] = annot_path + "anno_stat/kegg_stat/gene_kegg_table.xls"
            self.result_file['n_gene_sum_path' + 'ref'] = annot_path + "anno_stat/cog_stat/gene_cog_summary.xls"
            self.result_file['gene_stat_level2_ref'] = annot_path + "anno_stat/go_stat/gene_go12level_statistics.xls"
            self.result_file['gene_pathway_path_ref'] = annot_path + "anno_stat/kegg_stat/gene_pathway_table.xls"
        elif os.path.exists(os.path.join(annot_path, "anno_stat/trans_anno_detail.xls")):
            self.result_file['query_path' + "ref"] = os.path.join(annot_path, "anno_stat/trans_anno_detail.xls")
            self.annot_type = "denovo"
            self.result_file['gos_path' + '_ref'] = annot_path + "go/query_gos.list"
            self.result_file['gene_gos_path' + '_ref'] = annot_path + "anno_stat/go_stat/gene_gos.list"

            self.result_file['table_path' + '_ref'] = annot_path + "kegg/kegg_table.xls"
            self.result_file['gene_table_path' + '_ref'] = annot_path + "anno_stat/kegg_stat/gene_kegg_table.xls"
            self.result_file['n_gene_sum_path' + 'ref'] = annot_path + "anno_stat/cog_stat/gene_cog_summary.xls"
            self.result_file['gene_stat_level2_ref'] = annot_path + "anno_stat/go_stat/gene_go12level_statistics.xls"
            self.result_file['gene_pathway_path_ref'] = annot_path + "anno_stat/kegg_stat/gene_pathway_table.xls"
        else:
            self.annot_type = "annot2"
            self.result_file['ref_stat_path'] = os.path.join(annot_path, "all_stat.xls")
            self.result_file['ref_venn_path'] = os.path.join(annot_path)

            self.result_file['query_path' + "ref"] = os.path.join(annot_path, "all_annot.xls")

            self.result_file['gos_path' + '_ref'] = annot_path + "go/go_list_tran.xls"
            self.result_file['gene_gos_path' + '_ref'] = annot_path + "go/go_list_gene.xls"

            self.result_file['table_path' + '_ref'] = annot_path + "kegg/kegg_gene_tran.xls"
            self.result_file['gene_table_path' + '_ref'] = annot_path + "kegg/kegg_gene_gene.xls"
            self.result_file['n_gene_sum_path' + 'ref'] = annot_path + "cog/summary.G.tsv"
            # self.result_file['n_gene_sum_path' + 'new'] = annot_path + "cog/summary.G.tsv"
            self.result_file['gene_stat_level2_ref'] = annot_path + "go/go_lev2_gene.stat.xls"
            self.result_file['gene_pathway_path_ref'] = annot_path + "kegg/kegg_pathway_gene.xls"

        if os.path.exists(annot_path + "anno_stat/kegg_stat/gene_pathway"):
            self.result_file['gene_png_path_ref'] = annot_path + "anno_stat/kegg_stat/gene_pathway"
        elif os.path.exists(annot_path + "anno_stat/kegg_stat/pathways"):
            self.result_file['gene_png_path_ref'] = annot_path + "anno_stat/kegg_stat/pathways"
        elif os.path.exists(annot_path + "kegg/kegg_pathway_gene_dir"):
            self.result_file['gene_png_path_ref'] = annot_path + "kegg/kegg_pathway_gene_dir"


            
        for key, value in self.result_file.items():
            if self.has_new == False:
                if key.endswith("_ref") or key.startswith("ref_"):
                    if os.path.exists(value):
                        pass
                    else:
                        raise Exception('{}对应的结果文件{} 不存在，请检查'.format(key, value))
            else:
                if os.path.exists(value):
                    pass
                else:
                    raise Exception('{}对应的结果文件{} 不存在，请检查'.format(key, value))
        self.bind_object.logger.info("数据路径正确，文件完整 {}")

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


    def import_target_detail(self, new_target_file, known_target_file, params_dict, new_seq=None, known_seq=None, anno_type="origin", species_name=None, target_dir=None, version="v1"):
        saved_args = locals()
        self.bind_object.logger.info("import_target_detail {}".format(saved_args))
        # print("saved_args is", saved_args)
        self.version = version
        self.species_name = species_name
        self.bind_object.logger.info("开始导入靶基因")
        task_id = self.task_id
        result_dir = os.path.dirname(known_target_file)
        # print result_dir
        params_dict.update({
            "submit_location": "targe",
            "task_id": task_id,
            "task_type":2
        })

        self.remove_table_by_main_record(main_table='sg_target', task_id=task_id, detail_table=['sg_target_detail', 'sg_target_stat'], detail_table_key='target_id')
        target_id, columns = self.add_target(params=params_dict, name=None, species_name=species_name, result_dir=result_dir)
        columns = [x for x in columns if not x.startswith("paired")]
        self.add_target_detail(target_id, new_target_file, known_target_file, new_seq, known_seq, columns = columns, target_dir=target_dir)

    def import_target_detail_web(self, target_id, new_target_file, known_target_file, params_dict, new_seq=None, known_seq=None, anno_type="latest", species_name=None, last_id_target=None, target_dir=None, version="v1"):
        self.version = version
        self.bind_object.logger.info("开始导入靶基因")
        self.species_name = species_name
        task_id = self.task_id
        self.anno_type = anno_type
        result_dir = os.path.dirname(known_target_file)
        # print result_dir

        if last_id_target:
            last_id_target = ObjectId(last_id_target)
        else:
            pass

        if last_id_target:
            self.bind_object.logger.info("删除表格为 {}".format(last_id_target))
            self.remove_table_by_main_record(main_table='sg_target', _id=last_id_target, detail_table=['sg_target_detail', 'sg_target_stat'], detail_table_key='target_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")


        target_id, columns = self.add_target(params=params_dict, name=None, species_name=species_name, result_dir=result_dir, target_id=target_id)
        columns = [x for x in columns if not x.startswith("paired")]
        self.add_target_detail(target_id, new_target_file, known_target_file, new_seq, known_seq, columns = columns, target_dir=target_dir)


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
        if 'miranda' in params and params['miranda'] == "yes":
            columns.extend(['start_miranda', 'end_miranda', 'score_miranda', 'energy_miranda'])
            if self.version >= 'v1.1':
                columns.extend(['paired_miranda'])
            columns_select.extend(['score_miranda', 'energy_miranda'])
        if 'targetscan' in params and params['targetscan'] == "yes":
            columns.extend(['utr_start_targetscan', 'utr_end_targetscan', 'msa_start_targetscan', 'msa_end_targetscan', 'sead_match_targetscan'])
            if self.version >= 'v1.1':
                columns.extend(['pct', 'context_score'])
                columns_select.extend(['pct'])
        if 'psrobot' in params and params['psrobot'] == "yes":
            columns.extend(['start_psrobot', 'end_psrobot', 'score_psrobot'])
            if self.version >= 'v1.1':
                columns.extend(['paired_psrobot'])
            columns_select.extend(['score_psrobot'])

        if 'targetfinder' in params and params['targetfinder'] == "yes":
            columns.extend(['start_targetfinder', 'end_targetfinder', 'score_targetfinder'])
            if self.version >= 'v1.1':
                columns.extend(['paired_targetfinder'])
            columns_select.extend(['score_targetfinder'])
        if 'rnahybrid' in params and params['rnahybrid'] == "yes":
            columns.extend(['start_rnahybrid', 'end_rnahybrid', 'energy_rnahybrid', 'pvalue_rnahybrid'])
            if self.version >= 'v1.1':
                columns.extend(['paired_rnahybrid'])
            columns_select.extend(['pvalue_rnahybrid', 'energy_rnahybrid'])

        if self.version >= 'v2' and self.species_name in ['Homo_sapiens', 'Mus_musculus']:
            columns.extend(['starbase'])

        if target_id:
            self.db['sg_target'].update({"_id": target_id},
                    {"$set": {"columns": columns, "columns_select": columns_select, "main_id": target_id}})
            collection = self.db["sg_target"]
            result = collection.find_one({"_id": target_id})
            self.target_name = result['name']
            return target_id, columns

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        if name:
            self.target_name = name
        else:
            self.target_name = 'Target_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': self.target_name,
            'type': self.anno_type,
            'params': params,
            "version": "v1.1",
            'status': 'start',
            'desc': '注释统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'columns': columns,
            'columns_select': columns_select,
            'species_name': species_name,
            'result_dir': result_dir
            # 'seq_type': seq_type,
        }

        collection = self.db['sg_target']
        stat_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_target!")
        return stat_id, columns

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
                    if "yes" in stat:
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
                    if line.startswith("target 5"):
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
    def add_target_detail(self, target_id, new_target_file, known_target_file, new_seq=None, known_seq=None, columns = list(), target_dir=None):
        if not isinstance(target_id, ObjectId):
            if isinstance(target_id, types.StringTypes):
                target_id = ObjectId(target_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(new_target_file):
            raise Exception('{}所指定的路径不存在，请检查！'.format(new_target_file))
        if not os.path.exists(known_target_file):
            raise Exception('{}所指定的路径不存在，请检查！'.format(known_target_file))

        # header = ['small_rna', 'target'] + columns + ['mir_tar_base']
        data_list = []
        known_smallrna_list = []
        novol_smallrna_list = []
        known_target_list = []
        novol_target_list = []
        if columns[-1] == "starbase":
            columns = columns[: -1]
        header = ['small_rna', 'target', 'gene', 'name'] + columns + ['mirtarbase', 'starbase']

        # print "columns is {}".format(columns)
        # print "header is {}".format(header)

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
                data.append(('type', 'known'))
                data = SON(data)
                if self.version >= "v1.1":
                    for soft, pd_dict in paired_known.items():
                        if cols[0] + '|' + cols[1] in pd_dict:
                            data["paired_" + soft], data["pairlist_" + soft] = self.change_align_list2html(pd_dict[cols[0] + '|' + cols[1]], soft)
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
                data.append(('type', 'novel'))
                data = SON(data)
                if self.version >= "v1.1":
                    for soft, pd_dict in paired_novol.items():
                        if cols[0] + '|' + cols[1] in pd_dict:
                            data["paired_" + soft] , data["pairlist_" + soft] = self.change_align_list2html(pd_dict[cols[0] + '|' + cols[1]], soft)
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
            self.create_db_table('sg_target_detail', data_list)
            self.create_db_table('sg_target_stat', data_stat_list)
            self.db['sg_target'].update({"_id": target_id},
                                           {"$set": {"status": "end", "main_id": target_id}})
        except Exception as e:
            self.bind_object.set_error("导入靶基因预测统计信息出错!")
        else:
            self.bind_object.logger.info("导入靶基因预测统计信息成功!")


    def run(self, target_file, target_file2, annotation_mudule_dir, trans2gene, params_dict, taxon='Animals', exp_level='transcript', version="v1.0"):
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
        self.bind_object.logger.info("开始到表情数据路径为 {}".format(annotation_mudule_dir))
        self.import_target(target_file, target_file2, version=version)
        self.set_result_dir(annotation_mudule_dir)
        self.get_trans2gene(trans2gene)

        task_id = self.task_id

        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        # self.remove_table_by_main_record(main_table='sg_annotation_stat', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_stat_detail'], detail_table_key='stat_id')
        stat_id = self.add_annotation_stat(name=None, params=params, seq_type="new", database="nr,swissprot,pfam,cog,go,kegg", taxon=taxon, result_dir=self.result_dir, exp_level=exp_level)
        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['ref_stat_path'], venn_path=self.result_file['ref_venn_path'], seq_type = "ref", exp_level=exp_level)
        self.update_db_record('sg_annotation_stat', stat_id, status="end", main_id=stat_id)

        self.remove_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        if self.annot_type == "ref":
            self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path' + "ref"], seq_type = "ref", anno_type="T", exp_level=exp_level)
        else:
            self.add_annotation_query_denovo_detail2(query_id=query_id, query_path=self.result_file['query_path' + "ref"], seq_type = "ref", anno_type="T", exp_level=exp_level)
        self.update_db_record('sg_annotation_query',query_id, status="end", main_id=query_id)



    def run_web(self, target_file, target_file2, annotation_mudule_dir, trans2gene, params_dict, task_id, stat_id, last_id, taxon='Animals', exp_level='transcript', version='v1'):
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
        self.import_target(target_file, target_file2, version=version)
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

        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['ref_stat_path'], venn_path=self.result_file['ref_venn_path'], seq_type = "ref", exp_level=exp_level)
        self.update_db_record('sg_annotation_stat', stat_id, status="end", main_id=stat_id)
        # query

        if last_id:
            self.bind_object.logger.info("删除表格为 {}".format(last_id))
            self.remove_table_by_main_record(main_table='sg_annotation_stat', _id=last_id, detail_table=['sg_annotation_stat_detail'], detail_table_key='stat_id')
            self.bind_object.logger.info("删除表格成功 {}".format(last_id))

        query_old_id = self.get_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.anno_type)
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        if self.annot_type == "ref":
            self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path' + "ref"], seq_type = "ref", anno_type="T", exp_level=exp_level)
        else:
            self.add_annotation_query_denovo_detail2(query_id=query_id, query_path=self.result_file['query_path' + "ref"], seq_type = "ref", anno_type="T", exp_level=exp_level)
        self.update_db_record('sg_annotation_query', query_id, status="end", main_id=query_id)

        if query_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_query', _id=query_old_id['_id'], detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        # go
        cog_old_id = self.get_table_by_main_record(main_table='sg_annotation_cog', task_id=task_id, type=self.anno_type)
        params_select_cog = dict([(k,params_dict.get(k,None)) for k in ('cog_evalue', 'cog_similarity', 'cog_identity')])
        params_select_cog = json.dumps(params_select_cog, sort_keys=True, separators=(',', ':'))
        cog_id = self.add_annotation_cog(name=None, params=params_select_cog, result_dir=self.result_dir)
        if self.annot_type == "ref":
            self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_gene_sum_path' + 'ref'], seq_type="ref", anno_type="G")
        else:
            self.add_annotation_cog_detail2(cog_id=cog_id, cog_path=self.result_file['n_gene_sum_path' + 'ref'], seq_type="ref", anno_type="G")
        self.update_db_record('sg_annotation_cog', cog_id, status="end", main_id=cog_id)

        if cog_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_cog', _id=cog_old_id['_id'], detail_table=['sg_annotation_cog_detail'], detail_table_key='cog_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")


        go_old_id = self.get_table_by_main_record(main_table='sg_annotation_go', task_id=task_id, type=self.anno_type)
        params_select_nr = dict([(k,params_dict.get(k,None)) for k in ('nr_evalue', 'nr_similarity', 'nr_identity')])
        params_select_nr = json.dumps(params_select_nr, sort_keys=True, separators=(',', ':'))
        go_id = self.add_annotation_go(name=None, params=params_select_nr, result_dir=self.result_dir)


        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", gos_path=self.result_file['gos_path_ref'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", gos_path=self.result_file['gene_gos_path_ref'])
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2_ref'])
        self.update_db_record('sg_annotation_go', go_id, status="end", main_id=go_id)

        if go_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_go', _id=go_old_id['_id'], detail_table=['sg_annotation_go_detail', 'sg_annotation_go_graph','sg_annotation_go_level', 'sg_annotation_go_list'], detail_table_key='go_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        kegg_old_id = self.get_table_by_main_record(main_table='sg_annotation_kegg', task_id=task_id, type=self.anno_type)
        # kegg
        params_select_kegg = dict([(k,params_dict.get(k,None)) for k in ('kegg_evalue', 'kegg_similarity', 'kegg_identity')])
        params_select_kegg = json.dumps(params_select_kegg, sort_keys=True, separators=(',', ':'))
        kegg_id = self.add_annotation_kegg(name=None, params=None, result_dir=self.result_dir)

        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path_ref'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path_ref'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path_ref'], png_dir=self.result_file['gene_png_path_ref'])

        if kegg_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_kegg', _id=kegg_old_id['_id'], detail_table=['sg_annotation_kegg_categories', 'sg_annotation_kegg_level','sg_annotation_kegg_table', 'sg_annotation_kegg_pic'], detail_table_key='kegg_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

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
        self.bind_object.logger.info("开始导表task_id为 {}".format(task_id))
        # task_id = "denovo_rna_v2"
        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['ref_stat_path'], venn_path=self.result_file['ref_venn_path'], seq_type = "ref", exp_level=exp_level)
        if self.has_new:
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['new_stat_path'], venn_path=self.result_file['new_venn_path'], seq_type = "new", exp_level=exp_level)
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['all_stat_path'], venn_path=self.result_file['all_venn_path'], seq_type = "all", exp_level=exp_level)

        self.update_db_record('sg_annotation_stat', stat_id, status="end", main_id=stat_id)
        if last_id:
            self.bind_object.logger.info("删除表格为 {}".format(last_id))
            self.remove_table_by_main_record(main_table='sg_annotation_stat', _id=last_id, detail_table=['sg_annotation_stat_detail'], detail_table_key='stat_id')
            self.bind_object.logger.info("删除表格成功 {}".format(last_id))

        params_select_nr = dict([(k,params_dict.get(k,None)) for k in ('nr_evalue', 'nr_similarity', 'nr_identity')])
        params_select_nr = json.dumps(params_select_nr, sort_keys=True, separators=(',', ':'))

        nr_old_id = self.get_table_by_main_record(main_table='sg_annotation_nr', task_id=task_id, type=self.anno_type)
        self.bind_object.logger.info("查找表格{}".format(nr_old_id))
        nr_id = self.add_annotation_nr(name=None, params=params_select_nr, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_blast_nr_detail(blast_id=nr_id, seq_type="new", anno_type="T", database='nr', blast_path=self.result_file["nr" + "trans_blast_path"], exp_level=exp_level)
        self.add_annotation_blast_nr_detail(blast_id=nr_id, seq_type="ref", anno_type="T", database='nr', blast_path=self.result_file["nr" + "trans_blast_path" + "ref"], exp_level=exp_level)
        self.bind_object.logger.info("插入新表 {}".format(nr_id))
        self.update_db_record('sg_annotation_nr', nr_id, status="end", main_id=nr_id)


        if nr_old_id:
            self.bind_object.logger.info( "删除旧表 {}".format(nr_old_id['_id']))
            self.remove_table_by_main_record(main_table='sg_annotation_nr', _id=nr_old_id['_id'], detail_table=['sg_annotation_nr_detail', 'sg_annotation_nr_pie'], detail_table_key='nr_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        swissprot_old_id = self.get_table_by_main_record(main_table='sg_annotation_swissprot', task_id=task_id, type=self.anno_type)
        self.bind_object.logger.info("查找表格成{}".format(swissprot_old_id))
        params_select_swissprot = dict([(k,params_dict.get(k,None)) for k in ('swissprot_evalue', 'swissprot_similarity', 'swissprot_identity')])
        params_select_swissprot = json.dumps(params_select_swissprot, sort_keys=True, separators=(',', ':'))
        swissprot_id = self.add_annotation_swissprot(name=None, params=params_select_swissprot, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_blast_swissprot_detail(blast_id=swissprot_id, seq_type="new", anno_type="T", database='swissprot', blast_path=self.result_file["swissprot" + "trans_blast_path"], exp_level=exp_level)
        self.add_annotation_blast_swissprot_detail(blast_id=swissprot_id, seq_type="ref", anno_type="T", database='swissprot', blast_path=self.result_file["swissprot" + "trans_blast_path"  + "ref"], exp_level=exp_level)
        self.update_db_record('sg_annotation_swissprot', swissprot_id, status="end", main_id=swissprot_id)
        if swissprot_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_swissprot', _id=swissprot_old_id['_id'],  detail_table=['sg_annotation_swissprot_detail', 'sg_annotation_swissprot_pie'], detail_table_key='swissprot_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        pfam_old_id = self.get_table_by_main_record(main_table='sg_annotation_pfam', task_id=task_id, type=self.anno_type)
        params_select_pfam = dict([('pfam_evalue',params_dict['pfam_evalue'])])
        params_select_pfam = json.dumps(params_select_pfam, sort_keys=True, separators=(',', ':'))
        pfam_id = self.add_annotation_pfam(name=None, params=params_select_pfam, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T", exp_level=exp_level)
        self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'  + "ref"], seq_type="ref", anno_type="T", exp_level=exp_level)
        #self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path'], seq_type="new", anno_type="G")
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

        self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T", exp_level=exp_level)
        self.update_db_record('sg_annotation_pfam', pfam_id, status="end", main_id=pfam_id)

        if pfam_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_pfam', _id=pfam_old_id['_id'], detail_table=['sg_annotation_pfam_detail', 'sg_annotation_pfam_bar'], detail_table_key='pfam_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        query_old_id = self.get_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.anno_type)
        self.remove_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path'], seq_type = "new", anno_type="T", exp_level=exp_level)
        self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path' + "ref"], seq_type = "ref", anno_type="T", exp_level=exp_level)
        self.update_db_record('sg_annotation_query',query_id, status="end", main_id=query_id)

        if query_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_query', _id=query_old_id['_id'], detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        cog_old_id = self.get_table_by_main_record(main_table='sg_annotation_cog', task_id=task_id, type=self.anno_type)
        params_select_cog = dict([(k,params_dict.get(k,None)) for k in ('cog_evalue', 'cog_similarity', 'cog_identity')])
        params_select_cog = json.dumps(params_select_cog, sort_keys=True, separators=(',', ':'))
        cog_id = self.add_annotation_cog(name=None, params=params_select_cog, result_dir=self.result_dir)
        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_sum_path'], seq_type="new", anno_type="T")
            self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_gene_sum_path'], seq_type="new", anno_type="G")
        #r_sum_path = ref_anno_path + "/cog/cog_summary.xls"
        if exp_level.lower() == "transcript":
            self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_sum_path' + 'ref'], seq_type="ref", anno_type="T")
        #r_gene_sum_path = ref_anno_path + "/anno_stat/cog_stat/gene_cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_gene_sum_path' + 'ref'], seq_type="ref", anno_type="G")
        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=self.result_file['n_sum_path' + 'ref'], n_cog_path=self.result_file['n_sum_path'], seq_type="all", anno_type="T")
            self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=self.result_file['n_gene_sum_path' + 'ref'], n_cog_path=self.result_file['n_gene_sum_path'], seq_type="all", anno_type="G")
        self.update_db_record('sg_annotation_cog', cog_id, status="end", main_id=cog_id)
        if cog_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_cog', _id=cog_old_id['_id'], detail_table=['sg_annotation_cog_detail'], detail_table_key='cog_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        go_old_id = self.get_table_by_main_record(main_table='sg_annotation_go', task_id=task_id, type=self.anno_type)

        go_id = self.add_annotation_go(name=None, params=params_select_nr, result_dir=self.result_dir)
        if self.has_new:
            seq_type = "new"
            if exp_level.lower() == "transcript":
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2'])
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2'])
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3'])
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4'])
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", gos_path=self.result_file['gos_path'])
            self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4'])
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", gos_path=self.result_file['gene_gos_path'])
        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2_ref'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2_ref'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3_ref'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4_ref'])
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", gos_path=self.result_file['gos_path_ref'])
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4_ref'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3_ref'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4_ref'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", gos_path=self.result_file['gene_gos_path_ref'])

        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=2, r_go_path=self.result_file['stat_level2_all'])
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=3, r_go_path=self.result_file['stat_level3_all'])
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=4, r_go_path=self.result_file['stat_level4_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=2, r_go_path=self.result_file['gene_stat_level2_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=3, r_go_path=self.result_file['gene_stat_level3_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=4, r_go_path=self.result_file['gene_stat_level4_all'])
        self.update_db_record('sg_annotation_go', go_id, status="end", main_id=go_id)
        if go_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_go', _id=go_old_id['_id'], detail_table=['sg_annotation_go_detail', 'sg_annotation_go_graph','sg_annotation_go_level', 'sg_annotation_go_list'], detail_table_key='go_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        kegg_old_id = self.get_table_by_main_record(main_table='sg_annotation_kegg', task_id=task_id, type=self.anno_type)
        params_select_kegg = dict([(k,params_dict.get(k,None)) for k in ('kegg_evalue', 'kegg_similarity', 'kegg_identity')])
        params_select_kegg = json.dumps(params_select_kegg, sort_keys=True, separators=(',', ':'))
        kegg_id = self.add_annotation_kegg(name=None, params=params_select_kegg, result_dir=self.result_dir)
        if self.has_new:
            seq_type = "new"
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", categories_path=self.result_file['layer_path'])
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path'], png_dir=self.result_file['png_path'])
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path'])
                self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path'], png_dir=self.result_file['png_path'])
            self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", categories_path=self.result_file['gene_layer_path'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path'], png_dir=self.result_file['gene_png_path'])
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path'], png_dir=self.result_file['gene_png_path'])
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path'])

        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", categories_path=self.result_file['layer_path_ref'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path_ref'], png_dir=self.result_file['png_path_ref'])
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path_ref'], png_dir=self.result_file['png_path_ref'])
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path_ref'])
        self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", categories_path=self.result_file['gene_layer_path_ref'])
        self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path_ref'], png_dir=self.result_file['gene_png_path_ref'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path_ref'], png_dir=self.result_file['gene_png_path_ref'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path_ref'])

        if self.has_new:
            r_cate_path = self.result_file['layer_path_ref']
            n_cate_path = self.result_file['layer_path']
            r_gene_cate_path = self.result_file['gene_layer_path_ref']
            n_gene_cate_path = self.result_file['gene_layer_path']
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="T", r_cate_path=r_cate_path, n_cate_path=n_cate_path)
            self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="G", r_cate_path=r_gene_cate_path, n_cate_path=n_gene_cate_path)
            pathway_path = self.result_file['pathway_path' + '_all']
            png_path = self.result_file['png_path' + '_all']
            gene_pathway_path = self.result_file['gene_pathway_path' + '_all']
            gene_png_path =  self.result_file['gene_png_path' + '_all']
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="T", level_path=pathway_path, png_dir=png_path)
                self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="all", anno_type="T", level_path=pathway_path, png_dir=png_path)
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="G", level_path=gene_pathway_path, png_dir=gene_png_path)
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="all", anno_type="G", level_path=gene_pathway_path, png_dir=gene_png_path)

        self.update_db_record('sg_annotation_kegg', kegg_id, status="end", main_id=kegg_id)
        if kegg_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_kegg', _id=kegg_old_id['_id'], detail_table=['sg_annotation_kegg_categories', 'sg_annotation_kegg_level','sg_annotation_kegg_table'], detail_table_key='kegg_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

    def add_annotation(self, name=None, params=None, ref_anno_path=None, new_anno_path=None, pfam_path=None, merge_tran_output=None, merge_gene_output=None):
        """
        ref_anno_path: 已知序列注释的结果文件夹
        new_anno_path: 新序列注释的结果文件夹
        pfam_path:转录本的pfam_domain
        merge_tran_output: 转录本的merge_annot tool输出结果路径
        merge_gene_output: 基因的merge_annot tool输出结果路径
        """
        new_stat_path = new_anno_path + "/anno_stat/all_annotation_statistics.xls"
        new_venn_path = new_anno_path + "/anno_stat/venn"
        stat_id = self.add_annotation_stat(name=None, params=params, seq_type="new", database="nr,swissprot,pfam,cog,go,kegg")
        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=new_stat_path, venn_path=new_venn_path)
        blast_id = self.add_annotation_blast(name=None, params=params, stat_id=stat_id)
        blast_path = new_anno_path + "/anno_stat/blast"
        if os.path.exists(blast_path):
            for db in ["nr", "swissprot"]:
                trans_blast_path = blast_path + "/" + db + '.xls'
                gene_blast_path = blast_path + '/gene_' + db + '.xls'
                self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="transcript", database=db, blast_path=trans_blast_path)
                self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="gene", database=db, blast_path=gene_blast_path)
        else:
            raise Exception("没有blast的结果文件夹")
        nr_id = self.add_annotation_nr(name=None, params=params, stat_id=stat_id)
        nr_path = new_anno_path + "/anno_stat/blast_nr_statistics"
        if os.path.exists(nr_path):
            evalue_path = nr_path + "/nr_evalue.xls"
            similar_path = nr_path + "/nr_similar.xls"
            gene_evalue_path = nr_path + "/gene_nr_evalue.xls"
            gene_similar_path = nr_path + "/gene_nr_similar.xls"
            self.add_annotation_nr_pie(nr_id=nr_id, evalue_path=evalue_path, similar_path=similar_path, seq_type="new", anno_type="transcript")
            self.add_annotation_nr_pie(nr_id=nr_id, evalue_path=gene_evalue_path, similar_path=gene_similar_path,  seq_type="new", anno_type="gene")
        else:
            raise Exception("新序列NR注释结果文件不存在")
        swissprot_id = self.add_annotation_swissprot(name=None, params=params, stat_id=stat_id)
        swissprot_path = new_anno_path + "/anno_stat/blast_swissprot_statistics"
        if os.path.exists(swissprot_path):
            evalue_path = swissprot_path + "/swissprot_evalue.xls"
            similar_path = swissprot_path + "/swissprot_similar.xls"
            gene_evalue_path = swissprot_path + "/gene_swissprot_evalue.xls"
            gene_similar_path = swissprot_path + "/gene_swissprot_similar.xls"
            self.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=evalue_path, similar_path=similar_path, seq_type="new", anno_type="transcript")
            self.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=gene_evalue_path, similar_path=gene_similar_path, seq_type="new", anno_type="gene")
        else:
            raise Exception("新序列Swiss-Prot注释结果文件不存在")
        pfam_id = self.add_annotation_pfam(name=None, params=params, stat_id=stat_id)
        gene_pfam_path = new_anno_path + "/anno_stat/pfam_stat/gene_pfam_domain"
        if os.path.exists(pfam_path) and os.path.exists(gene_pfam_path):
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=pfam_path, seq_type="new", anno_type="transcript")
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=gene_pfam_path, seq_type="new", anno_type="gene")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=pfam_path, seq_type="new", anno_type="transcript")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=gene_pfam_path, seq_type="new", anno_type="gene")
        else:
            raise Exception("pfam注释结果文件不存在")
        ref_stat_path = ref_anno_path + "/anno_stat/all_annotation_statistics.xls"
        ref_venn_path = ref_anno_path + "/anno_stat/venn"
        if os.path.exists(ref_stat_path) and os.path.exists(ref_venn_path):
            stat_id = self.add_annotation_stat(name=None, params=params, seq_type="ref" , database="cog,go,kegg")
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=ref_stat_path, venn_path=ref_venn_path)
        else:
            raise Exception("已知序列注释统计文件和venn图文件夹不存在")
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id)
        query_path = ref_anno_path + "/anno_stat/trans_anno_detail.xls"
        gene_query_path = ref_anno_path + "/anno_stat/gene_anno_detail.xls"
        self.add_annotation_query_detail(query_id=query_id, query_path=query_path, anno_type="transcript")
        self.add_annotation_gene_query_detail(query_id=query_id, query_path=gene_query_path, anno_type="gene")
        query_path = new_anno_path + "/anno_stat/trans_anno_detail.xls"
        gene_query_path = new_anno_path + "/anno_stat/gene_anno_detail.xls"
        self.add_annotation_query_detail(query_id=query_id, query_path=query_path, anno_type="transcript")
        self.add_annotation_gene_query_detail(query_id=query_id, query_path=gene_query_path, anno_type="gene")
        cog_id = self.add_annotation_cog(name=name, params=params)
        r_sum_path = ref_anno_path + "/cog/cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=r_sum_path, seq_type="ref", anno_type="transcript")
        r_gene_sum_path = ref_anno_path + "/anno_stat/cog_stat/gene_cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=r_gene_sum_path, seq_type="ref", anno_type="gene")
        n_sum_path = new_anno_path + "/cog/cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=n_sum_path, seq_type="new", anno_type="transcript")
        n_gene_sum_path = new_anno_path + "/anno_stat/cog_stat/gene_cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=n_gene_sum_path, seq_type="new", anno_type="gene")
        self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=r_sum_path, n_cog_path=n_sum_path, seq_type="all", anno_type="transcript")
        self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=r_gene_sum_path, n_cog_path=n_gene_sum_path, seq_type="all", anno_type="gene")

        def add_go(go_id, go_path, gene_go_path, seq_type):
            if os.path.exists(go_path) and os.path.exists(gene_go_path):
                stat_level2 = go_path + "/go12level_statistics.xls"
                stat_level3 = go_path + "/go123level_statistics.xls"
                stat_level4 = go_path + "/go1234level_statistics.xls"
                gene_stat_level2 = gene_go_path + "/gene_go12level_statistics.xls"
                gene_stat_level3 = gene_go_path + "/gene_go123level_statistics.xls"
                gene_stat_level4 = gene_go_path + "/gene_go1234level_statistics.xls"
                gos_path = go_path + "/query_gos.list"
                gene_gos_path = gene_go_path + "/gene_gos.list"
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, level_path=stat_level2)
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, level_path=gene_stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, go_path=stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=3, go_path=stat_level3)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=4, go_path=stat_level4)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, go_path=gene_stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=3, go_path=gene_stat_level3)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=4, go_path=gene_stat_level4)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, go_path=stat_level2)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=3, go_path=stat_level3)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=4, go_path=stat_level4)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, go_path=gene_stat_level2)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=3, go_path=gene_stat_level3)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=4, go_path=gene_stat_level4)
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="transcript", gos_path=gos_path)
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="gene", gos_path=gene_gos_path)
            else:
                raise Exception("GO注释的结果文件不存在")
        go_id = self.add_annotation_go(name=name, params=params)
        go_path = ref_anno_path + "/go"
        gene_go_path = ref_anno_path + "/anno_stat/go_stat"
        add_go(go_id=go_id, go_path=go_path, gene_go_path=gene_go_path, seq_type="ref")
        go_path = new_anno_path + "/go"
        gene_go_path = new_anno_path + "/anno_stat/go_stat"
        add_go(go_id=go_id, go_path=go_path, gene_go_path=gene_go_path, seq_type="new")
        r_stat_level2 = ref_anno_path + "/go/go12level_statistics.xls"
        r_stat_level3 = ref_anno_path + "/go/go123level_statistics.xls"
        r_stat_level4 = ref_anno_path + "/go/go1234level_statistics.xls"
        n_stat_level2 = new_anno_path + "/go/go12level_statistics.xls"
        n_stat_level3 = new_anno_path + "/go/go123level_statistics.xls"
        n_stat_level4 = new_anno_path + "/go/go1234level_statistics.xls"
        r_gene_stat_level2 = ref_anno_path + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        r_gene_stat_level3 = ref_anno_path + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        r_gene_stat_level4 = ref_anno_path + "/anno_stat/go_stat/gene_go1234level_statistics.xls"
        n_gene_stat_level2 = new_anno_path + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        n_gene_stat_level3 = new_anno_path + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        n_gene_stat_level4 = new_anno_path + "/anno_stat/go_stat/gene_go1234level_statistics.xls"
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=2, r_go_path=r_stat_level2, n_go_path=n_stat_level2)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=3, r_go_path=r_stat_level3, n_go_path=n_stat_level3)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=4, r_go_path=r_stat_level4, n_go_path=n_stat_level4)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=2, r_go_path=r_gene_stat_level2, n_go_path=n_gene_stat_level2)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=3, r_go_path=r_gene_stat_level3, n_go_path=n_gene_stat_level3)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=4, r_go_path=r_gene_stat_level4, n_go_path=n_gene_stat_level4)

        def add_kegg(kegg_id, kegg_path, gene_kegg_path, seq_type):
            if os.path.exists(kegg_path) and os.path.exists(gene_kegg_path):
                layer_path = kegg_path + "/kegg_layer.xls"
                pathway_path = kegg_path + "/pathway_table.xls"
                png_path = kegg_path + "/pathways"
                table_path = kegg_path + "/kegg_table.xls"
                gene_layer_path = gene_kegg_path + "/gene_kegg_layer.xls"
                gene_pathway_path = gene_kegg_path + "/gene_pathway_table.xls"
                gene_png_path = gene_kegg_path + "/gene_pathway"
                gene_table_path = gene_kegg_path + "/gene_kegg_table.xls"
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", categories_path=layer_path)
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", categories_path=gene_layer_path)
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", level_path=pathway_path, png_dir=png_path)
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", level_path=gene_pathway_path, png_dir=gene_png_path)
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", table_path=table_path)
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", table_path=gene_table_path)
            else:
                raise Exception("KEGG注释文件不存在")
        kegg_id = self.add_annotation_kegg(name=None, params=params)
        kegg_path = ref_anno_path + "/kegg"
        gene_kegg_path = ref_anno_path + "/anno_stat/kegg_stat"
        add_kegg(kegg_id=kegg_id, kegg_path=kegg_path, gene_kegg_path=gene_kegg_path, seq_type="ref")
        kegg_path = new_anno_path + "/kegg"
        gene_kegg_path = new_anno_path + "/anno_stat/kegg_stat"
        add_kegg(kegg_id=kegg_id, kegg_path=kegg_path, gene_kegg_path=gene_kegg_path, seq_type="new")
        r_cate_path = ref_anno_path + "/kegg/kegg_layer.xls"
        n_cate_path = new_anno_path + "/kegg/kegg_layer.xls"
        r_gene_cate_path = ref_anno_path + "/anno_stat/kegg_stat/gene_kegg_layer.xls"
        n_gene_cate_path = new_anno_path + "/anno_stat/kegg_stat/gene_kegg_layer.xls"
        self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="transcript", r_cate_path=r_cate_path, n_cate_path=n_cate_path)
        self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="gene", r_cate_path=r_gene_cate_path, n_cate_path=n_gene_cate_path)
        pathway_path = merge_tran_output + "/pathway_table.xls"
        png_path = merge_tran_output + "/all_pathways"
        gene_pathway_path = merge_gene_output + "/pathway_table.xls"
        gene_png_path = merge_gene_output + "/all_pathways"
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="transcript", level_path=pathway_path, png_dir=png_path)
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="gene", level_path=gene_pathway_path, png_dir=gene_png_path)


    @report_check
    def add_annotation_stat(self, name=None, params=None, seq_type=None, database=None, result_dir=None, taxon='Animals', exp_level="transcript"):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if self.anno_type == "latest":
            anno_str = "new"
        else:
            anno_str = self.anno_type

        insert_data = {
            'exp_levl':exp_level.lower(),
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationStat_' + anno_str + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            "taxonomy": taxon,
            'species_name': self.species_name,
            'status': 'start',
            'desc': '注释统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            # 'seq_type': seq_type,
            'database': database
        }
        collection = self.db['sg_annotation_stat']
        stat_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_stat!")
        return stat_id

    @report_check
    def add_annotation_stat_detail(self, stat_id, stat_path, venn_path, seq_type, exp_level):
        """
        database: 进行统计的数据库
        stat_path: all_annotation_statistics.xls
        venn_path: venn图目录
        """
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(stat_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(stat_path))
        if not os.path.exists(venn_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(venn_path))
        data_list = []
        tail = "gene"
        if self.annot_type == "annot2":
            database_venn = {
                'NR': 'nr/nr_venn_{}.txt'.format(tail),
                'Swiss-Prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
                'Swiss-prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
                'Pfam': 'pfam/pfam_venn_{}.txt'.format(tail),
                'KEGG': 'kegg/kegg_venn_{}.txt'.format(tail),
                'GO': 'go/go_venn_{}.txt'.format(tail),
                'COG': 'cog/cog_venn_{}.txt'.format(tail),
            }
        else:
            database_venn = {
                'NR': 'nr_venn.txt',
                'Swiss-Prot': 'swissprot_venn.txt',
                'Swiss-prot': 'swissprot_venn.txt',
                'Pfam': 'pfam_venn.txt',
                'KEGG': 'kegg_venn.txt',
                'GO': 'go_venn.txt',
                'COG': 'cog_venn.txt',
            }

        target_trans = self.target_smallrna.keys()
        target_genes = list(set([self.tran2gene[x] for x in target_trans if x in self.tran2gene]))

        with open(stat_path, 'r') as f:
            lines = f.readlines()
            all_list = []
            for db in ['NR', 'Swiss-Prot', 'Pfam', 'KEGG', 'GO', 'COG']:
                data = [
                    ('stat_id', stat_id),
                    ('type', db),
                    ('seq_type', 'gene')
                ]

                venn_list, gene_venn_list = None, None
                database = ["nr", "swissprot", "pfam", "kegg", "go", "string", "cog"]
                if database_venn.has_key(db):
                    venn = venn_path + "/" + database_venn[db]
                    if os.path.exists(venn):
                        with open(venn, "rb") as f:
                            venn_list = [x.strip() for x in f.readlines()]
                        # print "****"
                        # print venn_list
                        target_venn_list = list(set(target_trans).intersection(set(venn_list)))
                        # print target_venn_list
                        #print self.tran2gene
                        target_gene_list = [self.tran2gene[x] for x in target_venn_list if x in self.tran2gene]
                        # print target_gene_list
                        target_gene_list = list(set(target_gene_list))
                        #print target_gene_list
                        all_list.extend(target_gene_list)

                        data.append(("gene_list", target_gene_list))
                        data.append(("gene", len(target_gene_list)))
                        if len(target_genes) == 0:
                            data.append(('gene_percent', round(0)))
                        else:
                            data.append(('gene_percent', round(float(len(target_gene_list)/float(len(target_genes))), 4)))

                    else:
                        raise Exception("{}对应的venn.txt文件{} {}不存在".format(line[0], venn, gene_venn))
                data = SON(data)
                data_list.append(data)
            if len(target_genes) == 0:
                pct = 0
            else:
                pct = round(float(len(set(all_list))/float(len(target_genes))), 4)
            data_all = [
                ('stat_id', stat_id),
                ('type', 'Total_anno'),
                ('seq_type', 'gene'),
                ("gene", len(set(all_list))),
                ('gene_percent', pct)
            ]
            data_list.append(SON(data_all))
            data_all2 = [
                ('stat_id', stat_id),
                ('type', 'Total'),
                ('seq_type', 'gene'),
                ("gene", len(set(target_genes))),
                ('gene_percent', '1.0')
            ]
            data_list.append(SON(data_all2))
        try:
            collection = self.db['sg_annotation_stat_detail']
            collection.insert_many(data_list)
        except Exception, e:
            raise Exception("导入注释统计信息失败:%s" % stat_path)
        else:
            self.bind_object.logger.info("导入注释统计信息成功：%s" % (stat_path))

    @report_check
    def add_stat_detail(self, old_stat_id, stat_id, nr_evalue, gene_nr_evalue, sw_evalue, gene_sw_evalue):
        """
        注释重运行时注释统计导表sg_annotation_stat_detail
        """
        if not isinstance(old_stat_id, ObjectId):
            if isinstance(old_stat_id, types.StringTypes):
                old_stat_id = ObjectId(old_stat_id)
            else:
                raise Exception('old_stat_id必须为ObjectId对象或其对应的字符串！')
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        collection = self.db["sg_annotation_stat_detail"]
        results = collection.find({"stat_id": old_stat_id})
        data_list, data = [], []
        for result in results:
            db = result["type"]
            if db == "total":
                total_tran = result["transcript"]
                total_gene = result["gene"]
            # 增加数据库类型，区分有注释的基因和有分类的基因
            # if db in ["pfam", "total_anno", "total_anno_nsp","total_class", "total"]:
            if db in ["pfam", "total_anno",  "total"]:
                data = [
                    ('stat_id', stat_id),
                    ('type', result["type"]),
                    ('transcript', result["transcript"]),
                    ('gene', result["gene"]),
                    ('transcript_percent', result["transcript_percent"]),
                    ('gene_percent', result["gene_percent"]),
                    ('gene_list', result["gene_list"]),
                    ('transcript_list', result["transcript_list"])
                ]
                data = SON(data)
                data_list.append(data)
        nr_ids = self.stat(stat_path=nr_evalue)
        gene_nr_ids = self.stat(stat_path=gene_nr_evalue)
        data = [
            ('stat_id', stat_id),
            ('type', "nr"),
            ('transcript', len(nr_ids)),
            ('gene', len(gene_nr_ids)),
            ('transcript_percent', round(float(len(nr_ids))/total_tran, 4)),
            ('gene_percent', round(float(len(gene_nr_ids))/total_gene, 4)),
            ('gene_list', ",".join(gene_nr_ids)),
            ('transcript_list', ",".join(nr_ids))
        ]
        data = SON(data)
        data_list.append(data)
        sw_ids = self.stat(stat_path=sw_evalue)
        gene_sw_ids = self.stat(stat_path=gene_sw_evalue)
        data = [
            ('stat_id', stat_id),
            ('type', "swissprot"),
            ('transcript', len(sw_ids)),
            ('gene', len(gene_sw_ids)),
            ('transcript_percent', round(float(len(sw_ids))/total_tran, 4)),
            ('gene_percent', round(float(len(gene_sw_ids))/total_gene, 4)),
            ('gene_list', ",".join(gene_sw_ids)),
            ('transcript_list', ",".join(sw_ids))
        ]
        data = SON(data)
        data_list.append(data)
        try:
            collection = self.db['sg_annotation_stat_detail']
            collection.insert_many(data_list)
        except:
            raise Exception("导入注释统计信息出错")
        else:
            self.bind_object.logger.info("导入注释统计信息成功")

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

    @report_check
    def add_annotation_blast(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationBlast_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'blast最佳比对结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_blast']
        blast_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_blast!")
        return blast_id

    @report_check
    def add_annotation_blast_detail(self, blast_id, seq_type, anno_type, database, blast_path):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                raise Exception('blast_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(blast_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(blast_path))
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            flag = None
            for line in lines[1:]:
                line = line.strip().split('\t')
                query_name = line[5]
                hit_name = line[10]
                if flag == query_name:
                    pass
                else:
                    flag = query_name
                    data = {
                        'blast_id': blast_id,
                        'seq_type': seq_type,
                        'anno_type': anno_type,
                        'database': database,
                        'score': float(line[0]),
                        'e_value': float(line[1]),
                        'hsp_len': int(line[2]),
                        'identity_rate': round(float(line[3]), 4),
                        'similarity_rate': round(float(line[4]), 4),
                        'query_id': line[5],
                        'q_len': int(line[6]),
                        'q_begin': line[7],
                        'q_end': line[8],
                        'q_frame': line[9],
                        'hit_name': line[10],
                        'hit_len': int(line[11]),
                        'hsp_begin': line[12],
                        'hsp_end': line[13],
                        'hsp_frame': line[14],
                        'description': line[15]
                    }
                    if  anno_type == 'T':
                        try:
                            data.update({'gene_id': self.trans_gene[line[5]]})
                            data.update({'is_gene': self.trans_isgene[line[5]]})
                        except:
                            data.update({'gene_id': line[5]})
                            data.update({'is_gene': True})
                    collection = self.db['sg_annotation_blast_detail']
                    collection.insert_one(data).inserted_id
        self.bind_object.logger.info("导入blast信息：%s成功!" % (blast_path))

    @report_check
    def add_annotation_blast_nr_detail(self, blast_id, seq_type, anno_type, database, blast_path, exp_level):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                raise Exception('blast_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(blast_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(blast_path))
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            flag = None
            for line in lines[1:]:
                line = line.strip().split('\t')
                query_name = line[5]
                hit_name = line[10]
                if flag == query_name:
                    pass
                else:
                    flag = query_name
                    data = {
                        'nr_id': blast_id,
                        'seq_type': seq_type,
                        'database': database,
                        'score': float(line[0]),
                        'e_value': float(line[1]),
                        'hsp_len': int(line[2]),
                        'identity_rate': round(float(line[3]), 4),
                        'similarity_rate': round(float(line[4]), 4),
                        'transcript_id': line[5],
                        'q_len': int(line[6]),
                        'hit_name': line[10],
                        'description': line[15]
                    }
                    if anno_type == 'T':
                        try:
                            data.update({'gene_id': self.trans_gene[line[5]]})
                            data.update({'is_gene': self.trans_isgene[line[5]]})
                        except:
                            data.update({'gene_id': line[5]})
                            data.update({'is_gene': True})
                    else:
                        pass
                    if exp_level.lower() == "gene" and self.trans_isgene[line[5]] == False:
                        continue
                    else:
                        data_list.append(SON(data))
        collection = self.db['sg_annotation_nr_detail']
        if data_list:
            try:
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入nrblast信息：%s成功!" % (blast_path))
            except Exception, e:
                raise Exception("导入注释统计信息失败:%s" % blast_path)


    @report_check
    def add_annotation_blast_swissprot_detail(self, blast_id, seq_type, anno_type, database, blast_path, exp_level):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                raise Exception('blast_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(blast_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(blast_path))
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            flag = None
            for line in lines[1:]:
                line = line.strip().split('\t')
                query_name = line[5]
                hit_name = line[10]
                if flag == query_name:
                    pass
                else:
                    flag = query_name
                    data = {
                        'swissprot_id': blast_id,
                        'seq_type': seq_type,
                        # 'anno_type': anno_type,
                        'database': database,
                        'score': float(line[0]),
                        'e_value': float(line[1]),
                        'hsp_len': int(line[2]),
                        'identity_rate': round(float(line[3]), 4),
                        'similarity_rate': round(float(line[4]), 4),
                        'transcript_id': line[5],
                        'q_len': int(line[6]),
                        # 'q_begin': line[7],
                        # 'q_end': line[8],
                        # 'q_frame': line[9],
                        'hit_name': line[10],
                        # 'hit_len': int(line[11]),
                        # 'hsp_begin': line[12],
                        # 'hsp_end': line[13],
                        # 'hsp_frame': line[14],
                        'description': line[15]
                    }
                    if anno_type == 'T':
                        try:
                            data.update({'gene_id': self.trans_gene[line[5]]})
                            data.update({'is_gene': self.trans_isgene[line[5]]})
                        except:
                            data.update({'gene_id': line[5]})
                            data.update({'is_gene': True})
                    else:
                        pass
                    if exp_level.lower() == "gene" and self.trans_isgene[line[5]] == False:
                        continue
                    else:
                        data_list.append(SON(data))
        collection = self.db['sg_annotation_swissprot_detail']
        if data_list:
            try:
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入swissprot blast信息：%s成功!" % (blast_path))
            except Exception, e:
                raise Exception("导入swissprot注释统计信息失败:%s" % blast_path)

    @report_check
    def add_annotation_nr(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationNr_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'nr注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_nr']
        nr_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_nr!")
        return nr_id

    @report_check
    def add_annotation_nr_pie(self, nr_id, species_path, evalue_path, similar_path, seq_type, anno_type):
        if not isinstance(nr_id, ObjectId):
            if isinstance(nr_id, types.StringTypes):
                nr_id = ObjectId(nr_id)
            else:
                raise Exception('nr_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(evalue_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(evalue_path))
        if not os.path.exists(similar_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(similar_path))
        evalue, evalue_list,species_list_a, species_a, similar, similar_list = [], [], [], [], [], []
        data_list = []
        with open(evalue_path, "r") as f1, open(similar_path, "r") as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
            for line1 in lines1[1:]:
                line1 = line1.strip().split('\t')
                value = {"key": line1[0], "value": int(line1[1])}
                try:
                    value_list = {"key": line1[0], "value": line1[2]}
                except:
                    value_list = {"key": line1[0], "value": None}
                evalue.append(value)
                evalue_list.append(value_list)
            for line2 in lines2[1:]:
                line2 = line2.strip().split('\t')
                similarity = {"key": line2[0], "value": int(line2[1])}
                try:
                    similarity_list = {"key": line2[0], "value": line2[2]}
                except:
                    similarity_list = {"key": line2[0], "value": None}
                similar.append(similarity)
                similar_list.append(similarity_list)

        with open(species_path, "r") as f3:
            lines3 = f3.readlines()
            for line3 in lines3[1:15]:
                line3 = line3.strip().split("\t")
                if anno_type == "T":
                    species = {"key": line3[0], "value": int(line3[1])}
                    try:
                        species_list = {"key": line3[0], "value": line3[5]}
                    except:
                        species_list = {"key": line3[0], "value": None}
                elif anno_type == "G":
                    species = {"key": line3[0], "value": int(line3[2])}
                    try:
                        species_list = {"key": line3[0], "value": line3[6]}
                    except:
                        species_list = {"key": line3[0], "value": None}
                species_a.append(species)
                species_list_a.append(species_list)
            other = 0
            for line3 in lines3[16:]:
                line3 = line3.strip().split("\t")
                if anno_type == "T":
                    other += int(line3[1])
                elif anno_type == "G":
                    other += int(line3[2])
            species = {"key": 'other', "value": other}
            species_list = {"key": 'other', "value": None}
            species_a.append(species)
            species_list_a.append(species_list)

            data = [
                ('nr_id', nr_id),
                #('seq_type', seq_type),
                ('anno_type', anno_type),
                ('e_value', evalue),
                ('similar', similar),
                ('evalue_list', evalue_list),
                ('similar_list', similar_list),
                ('species', species_a),
                ('species_list', species_list_a),
            ]

        data = SON(data)
        data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_nr_pie']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入nr库注释作图信息evalue,similar：%s、%s出错!" % (evalue_path, similar_path))
            else:
                self.bind_object.logger.info("导入nr库注释作图信息evalue,similar：%s、%s成功!" % (evalue_path, similar_path))

    @report_check
    def add_annotation_swissprot(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationSwissprot_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'swissprot注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_swissprot']
        swissprot_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_swissprot!")
        return swissprot_id

    @report_check
    def add_annotation_swissprot_pie(self, swissprot_id, evalue_path, similar_path, seq_type, anno_type):
        """
        """
        if not isinstance(swissprot_id, ObjectId):
            if isinstance(swissprot_id, types.StringTypes):
                swissprot_id = ObjectId(swissprot_id)
            else:
                raise Exception('swissprot_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(evalue_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(evalue_path))
        if not os.path.exists(similar_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(similar_path))
        evalue, evalue_list, similar, similar_list = [], [], [], []
        data_list = []
        with open(evalue_path, "r") as f1, open(similar_path, "r") as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
            for line1 in lines1[1:]:
                line1 = line1.strip().split('\t')
                value = {"key": line1[0], "value": int(line1[1])}
                try:
                    value_list = {"key": line1[0], "value": line1[2]}
                except:
                    value_list = {"key": line1[0], "value": None}
                evalue.append(value)
                evalue_list.append(value_list)
            for line2 in lines2[1:]:
                line2 = line2.strip().split('\t')
                similarity = {"key": line2[0], "value": int(line2[1])}
                try:
                    similarity_list = {"key": line2[0], "value": line2[2]}
                except:
                    similarity_list = {"key": line2[0], "value": None}
                similar.append(similarity)
                similar_list.append(similarity_list)
        data = [
            ('swissprot_id', swissprot_id),
            # ('seq_type', seq_type),
            ('anno_type', anno_type),
            ('e_value', evalue),
            ('similar', similar),
            ('evalue_list', evalue_list),
            ('similar_list', similar_list),
        ]
        data = SON(data)
        data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_swissprot_pie']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入swissprot库注释作图信息evalue,similar：%s、%s出错!" % (evalue_path, similar_path))
            else:
                self.bind_object.logger.info("导入swissprot库注释作图信息evalue,similar：%s、%s成功!" % (evalue_path, similar_path))

    @report_check
    def add_annotation_pfam(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationPfam_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'pfam注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_pfam']
        pfam_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_pfam!")
        return pfam_id

    def add_annotation_pfam_detail(self, pfam_id, pfam_path, seq_type, anno_type, exp_level):
        """
        pfam_path: pfam_domain
        """
        if not isinstance(pfam_id, ObjectId):
            if isinstance(pfam_id, types.StringTypes):
                pfam_id = ObjectId(pfam_id)
            else:
                raise Exception('pfam_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(pfam_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(pfam_path))
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
                raise Exception("导入pfam注释信息:%s失败！" % pfam_path)
            else:
                self.bind_object.logger.info("导入pfam注释信息:%s成功" % pfam_path)

    @report_check
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
                raise Exception('pfam_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(pfam_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(pfam_path))
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
                raise Exception("导入pfam注释信息:%s失败！" % pfam_path)
            else:
                self.bind_object.logger.info("导入pfam注释信息:%s成功！" % pfam_path)

    @report_check
    def add_annotation_cog(self, name=None, params=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationCog_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'cog注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        collection = self.db['sg_annotation_cog']
        cog_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_cog!")
        return cog_id

    @report_check
    def add_annotation_cog_detail2(self, cog_id, cog_path, seq_type, anno_type):
        '''
        cog_path: cog_summary.xls
        seq_type: ref/new
        anno_type: transcript/gene
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                raise Exception('cog_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(cog_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(cog_path))
        data_list = list()
        with open(cog_path, 'r') as f:
            lines = f.readlines()
            for line in lines[2:]:
                line = line.strip().split('\t')
                data = [
                    ('cog_id', cog_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('type', line[0]),
                    ('function_categories', line[1]),
                    ('cog', int(len(line[4].split(";")))),
                 ]
                try:
                    data.append(('cog_list', line[4]))
                except:
                    data.append(('cog_list', None))
                # try:
                #     data.append(('nog_list', line[5]))
                # except:
                #     data.append(('nog_list', None))
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_cog_detail']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入cog注释信息：%s出错!" % (cog_path))
            else:
                self.bind_object.logger.info("导入cog注释信息：%s成功!" % (cog_path))


    @report_check
    def add_annotation_cog_detail(self, cog_id, cog_path, seq_type, anno_type):
        '''
        cog_path: cog_summary.xls
        seq_type: ref/new
        anno_type: transcript/gene
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                raise Exception('cog_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(cog_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(cog_path))
        data_list = list()
        with open(cog_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('cog_id', cog_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('type', line[0]),
                    ('function_categories', "[" + line[2] + "]" + " " + line[1]),
                    ('cog', int(line[3])),
                 ]
                try:
                    data.append(('cog_list', line[4]))
                except:
                    data.append(('cog_list', None))
                # try:
                #     data.append(('nog_list', line[5]))
                # except:
                #     data.append(('nog_list', None))
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_cog_detail']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入cog注释信息：%s出错!" % (cog_path))
            else:
                self.bind_object.logger.info("导入cog注释信息：%s成功!" % (cog_path))

    @report_check
    def add_annotation_cog_detail3(self, cog_id, cog_path, seq_type, anno_type):
        '''
        cog_path: cog_table
        seq_type: ref/new
        anno_type: transcript/gene
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                raise Exception('cog_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(cog_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(cog_path))
        data_list = list()
        with open(cog_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('cog_id', cog_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('type', line[0]),
                    ('function_categories', "[" + line[2] + "]" + " " + line[1]),
                    ('cog', int(line[3])),
                 ]
                try:
                    data.append(('cog_list', line[4]))
                except:
                    data.append(('cog_list', None))
                # try:
                #     data.append(('nog_list', line[5]))
                # except:
                #     data.append(('nog_list', None))
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_cog_detail']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入cog注释信息：%s出错!" % (cog_path))
            else:
                self.bind_object.logger.info("导入cog注释信息：%s成功!" % (cog_path))

    @report_check
    def add_annotation_cog_detail_all(self, cog_id, r_cog_path, n_cog_path, seq_type, anno_type):
        '''
        r_cog_path: cog_summary.xls(ref)
        n_cog_path: cog_summary.xls(new)
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                raise Exception('cog_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(r_cog_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(r_cog_path))
        if not os.path.exists(n_cog_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(n_cog_path))
        first = ['INFORMATION STORAGE AND PROCESSING', 'CELLULAR PROCESSES AND SIGNALING', 'METABOLISM', 'POORLY CHARACTERIZED']
        func_type = {
            'INFORMATION STORAGE AND PROCESSING': sorted(['J', 'A', 'K', 'L', 'B']),
            'CELLULAR PROCESSES AND SIGNALING': sorted(['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O']),
            'METABOLISM': sorted(['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q']),
            'POORLY CHARACTERIZED': sorted(['R', 'S']),
        }
        func_decs = {
            'J': 'Translation, ribosomal structure and biogenesis',
            'A': 'RNA processing and modification', 'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'B': 'Chromatin structure and dynamics',
            'D': 'Cell cycle control, cell division, chromosome partitioning',
            'Y': 'Nuclear structure', 'V': 'Defense mechanisms', 'T': 'Signal transduction mechanisms',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility', 'Z': 'Cytoskeleton', 'W': 'Extracellular structures',
            'U': 'Intracellular trafficking, secretion, and vesicular transport',
            'O': 'Posttranslational modification, protein turnover, chaperones',
            'C': 'Energy production and conversion', 'G': 'Carbohydrate transport and metabolism',
            'E': 'Amino acid transport and metabolism', 'F': 'Nucleotide transport and metabolism',
            'H': 'Coenzyme transport and metabolism', 'I': 'Lipid transport and metabolism',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
            'R': 'General function prediction only', 'S': 'Function unknown'
        }
        data_list = list()
        cate = {}
        funlist = {'COG': {}, 'NOG': {}}
        cog_fun, nog_fun = {}, {}
        with open(r_cog_path, 'r') as f, open(n_cog_path, 'r') as n:
            lines = f.readlines()
            items = n.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                line[1] = "[" + line[2] + "]" + " " + line[1]
                m = re.match(r"\[(.+)\].+$", line[1])
                if m:
                    fun1 = m.group(1)
                    try:
                        funlist['COG'][fun1] = line[4].split(";")
                    except:
                        funlist['COG'][fun1] = []
                    try:
                        funlist['NOG'][fun1] = line[5].split(";")
                    except:
                        funlist['NOG'][fun1] = []
            for item in items[1:]:
                item = item.strip().split('\t')
                item[1] = "[" + item[2] + "]" + " " + item[1]
                m = re.match(r"\[(.+)\].+$", item[1])
                if m:
                    fun1 = m.group(1)
                    if fun1 in funlist['COG']:
                        try:
                            cog_ids = item[4].split(";")
                            if cog_ids:
                                for cog in cog_ids:
                                    if cog not in funlist['COG'][fun1]:
                                        funlist['COG'][fun1].append(cog)
                        except:
                            pass
                    else:
                        try:
                            funlist['COG'][fun1] = item[4].split(";")
                        except:
                            funlist['COG'][fun1] = []
                    if fun1 in funlist['NOG']:
                        try:
                            nog_ids = item[5].split(";")
                            if nog_ids:
                                for nog in nog_ids:
                                    if nog not in funlist['NOG'][fun1]:
                                        funlist['NOG'][fun1].append(nog)
                        except:
                            pass
                    else:
                        try:
                            funlist['NOG'][fun1] = item[5].split(";")
                        except:
                            funlist['NOG'][fun1] = []
        for g in funlist['COG'].keys():
            detail = func_decs[g]
            category = '[' + g + ']' + ' ' + detail
            try:
                cog_list = list(set(funlist['COG'][g]))
            except:
                cog_list = []
            try:
                nog_list = list(set(funlist['NOG'][g]))
            except:
                nog_list = []
            for cog_class in func_type.keys():
                if g in func_type[cog_class]:
                    cog_type = cog_class
            data = [
                ('cog_id', cog_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('type', cog_type),
                ('function_categories', category),
                ('cog', len([x for x in cog_list if x ])),
                # ('nog', len([x for x in nog_list if x])),
                ('cog_list', ';'.join(cog_list)),
                # ('nog_list', ';'.join(nog_list)),
            ]
            data = SON(data)
            if len(cog_list) + len(nog_list) != 0:
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_cog_detail']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入cog注释all出错：%s， %s" % (r_cog_path, n_cog_path))
            else:
                self.bind_object.logger.info("导入cog注释all成功：%s, %s" % (r_cog_path, n_cog_path))

    @report_check
    def add_annotation_cog_table(self, cog_id, table_path, seq_type, anno_type):
        '''
        table_path:cog_table.xls
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                raise Exception('cog_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(table_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(table_path))
        data_list = list()
        with open(table_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('cog_id', cog_id),
                    # ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('query_name', line[0]),
                    ('query_length', line[1]),
                    ('hsp_start', line[2]),
                    ('hsp_end', line[3]),
                    ('hsp_strand', line[4]),
                    ('hit_name', line[5]),
                    ('hit_description', line[6]),
                    ('hit_length', line[7]),
                    ('hit_start', line[8]),
                    ('hit_end', line[9]),
                    ('group', line[10]),
                    ('group_description', line[11]),
                    ('group_categories', line[12]),
                    ('region_start', line[13]),
                    ('region_end', line[14]),
                    ('region_coverage', line[15]),
                    ('region_identities', line[16]),
                    ('region_positives', line[17])
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_cog_table']
            collection.insert_many(data_list)
        except Exception, e:
            raise Exception("导入cog注释table信息：%s出错!" % (table_path))
        else:
            self.bind_object.logger.info("导入cog注释table信息：%s成功!" % (table_path))

    @report_check
    def add_annotation_go(self, name=None, params=None, result_dir=None):
        """
        go注释导表函数
        """
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationGo_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'go注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        collection = self.db['sg_annotation_go']
        go_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_go!")
        return go_id

    def add_annotation_go_detail(self, go_id, seq_type, anno_type, level, go_path):
        """
        go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                raise Exception('go_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(go_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(go_path))
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
                    raise Exception("导入go注释信息：%s出错!" % (go_path))
                else:
                    self.bind_object.logger.info("导入go注释信息：%s成功!" % (go_path))

    @report_check
    def add_annotation_go_graph(self, go_id, seq_type, anno_type, level, go_path):
        """
        go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                raise Exception('go_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(go_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(go_path))
        data_list = list()
        self.bind_object.logger.info("开始导入go注释画图信息：%s" % (go_path))
        with open(go_path, 'r') as f:
            lines = f.readlines()
            term = {}
            term_list = []
            for i in range(1, len(lines)):
                line = lines[i].strip().split('\t')
                if level == 2:
                    term_type = line[0]
                    go_term = line[1]
                    if go_term not in term:
                        term[go_term] = []
                        term[go_term].append(i)
                    else:
                        term[go_term].append(i)
                if level == 3:
                    term_type = line[0]
                    go_term = line[3]
                    if go_term not in term:
                        term[go_term] = []
                        term[go_term].append(i)
                    else:
                        term[go_term].append(i)
                if level == 4:
                    term_type = line[0]
                    go_term = line[5]
                    if go_term not in term:
                        term[go_term] = []
                        term[go_term].append(i)
                    else:
                        term[go_term].append(i)
        with open(go_path, 'r') as f:
            lines = f.readlines()
            for item in term:
                seq_list = []
                for j in term[item]:
                    line = lines[j].strip().split("\t")
                    term_type = line[0]
                    '''
                    for seq in line[-1].split(";"):
                        if seq not in seq_list:
                            seq_list.append(seq)
                    '''
                data = [
                    ('go_id', go_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('level', level),
                    ('term_type', term_type),
                    ('go_term', item),
                    ('seq_number', line[-3]),
                    ('percent', line[-2]),
                    #('seq_list', seq_list)
                ]
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_go_graph']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入go注释画图信息出错：%s" % (go_path))
            else:
                self.bind_object.logger.info("导入go注释画图信息成功：%s" % (go_path))

    @report_check
    def add_annotation_go_level(self, go_id, seq_type, anno_type, level, level_path):
        """
        level_path: go2level.xls
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                raise Exception('go_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(level_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(level_path))
        data_list = list()
        with open(level_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('go_id', go_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('level', level),
                    ('term_type', line[1]),
                    ('parent_name', line[0]),
                    ('num', int(line[3])),
                    ('percent', round(float(line[4]), 4)),
                    ('go', line[2]),
                    #('seq_list', line[5]),
                ]
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_go_level']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入go注释第二层级信息：%s出错!" % (level_path))
            else:
                self.bind_object.logger.info("导入go注释第二层级信息：%s成功!" % (level_path))

    @report_check
    def add_annotation_go_list(self, go_id, seq_type, anno_type, gos_path):
        """
        gos_path: go.list
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                raise Exception('go_id须为ObjectId对象或其他对应的字符串！')
        if not os.path.exists(gos_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(gos_path))
        data_list = []
        with open(gos_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                data = [
                    ('go_id', go_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('gene_id', line[0]),
                    ('gos_list', line[1]),
                ]
                data = SON(data)
                data_list.append(data)
            if data_list:
                try:
                    collection = self.db['sg_annotation_go_list']
                    collection.insert_many(data_list)
                except Exception, e:
                    raise Exception("导入gos_list注释信息：%s出错:%s" % (gos_path))
                else:
                    self.bind_object.logger.info("导入gos_list注释信息：%s成功!" % (gos_path))

    def add_annotation_go_all(self, go_id, seq_type, anno_type, level, r_go_path):
        """
        r_go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls(ref)
        n_go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls(new)
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                raise Exception('go_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(r_go_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(r_go_path))
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
                self.bind_object.set_error("导入go注释画图all信息出错：%s" % (r_go_path))
                print "导入go注释画图all出错：%s" % (r_go_path)
            else:
                self.bind_object.logger.info("导入go注释画图all信息成功：%s" % (r_go_path))
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
                raise Exception("导入go注释all出错：%s" % (r_go_path))
            else:
                self.bind_object.logger.info("导入go注释all信息成功：%s" % (r_go_path))

    @report_check
    def add_annotation_kegg(self, name=None, params=None, result_dir=None):
        """
        kegg注释导表函数
        """
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationKegg_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'kegg注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'categories': ["M","CP","EIP","GIP","OS","HD","DD"]
        }
        collection = self.db['sg_annotation_kegg']
        kegg_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_kegg!")
        return kegg_id

    @report_check
    def add_annotation_kegg_categories(self, kegg_id, seq_type, anno_type, categories_path):
        """
        categories_path:kegg_layer.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(categories_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(categories_path))
        data_list = list()
        first_type = []
        with open(categories_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                type_abr = ''.join([x[0] for x in line[0].split(' ')])
                first_type.append(type_abr)
                data = [
                    ('kegg_id', kegg_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('first_category', line[0]),
                    ('second_category', line[1]),
                    ('num', int(line[2])),
                    ('seq_list', line[3]),
                ]
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_categories']
                collection.insert_many(data_list)
                self.update_db_record('sg_annotation_kegg', kegg_id, categories=list(set(first_type)))
            except Exception, e:
                raise Exception("导入kegg注释分类信息：%s出错!" % (categories_path))
            else:
                self.bind_object.logger.info("导入kegg注释分类信息：%s 成功!" % categories_path)

    @report_check
    def add_annotation_kegg_level(self, kegg_id, seq_type, anno_type, level_path, png_dir):
        """
        level_path: pathway_table.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(level_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(level_path))
        if not os.path.exists(png_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(png_dir))
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                # fs = gridfs.GridFS(self.db)
                pid = re.sub('path:', '', line[0])
                # pdfid = fs.put(open(png_dir + '/' + pid + '.pdf', 'rb'))
                # graph_png_id = fs.put(open(png_dir + '/' + pid + '.png', 'rb'))
                insert_data = {
                    'kegg_id': kegg_id,
                    'seq_type': seq_type,
                    'anno_type': anno_type,
                    'pathway_id': line[0],
                    'first_category': line[1],
                    'second_category': line[2],
                    'pathway_definition': line[3],
                    'number_of_seqs': int(line[4]),
                    'seq_list': line[5],
                    # 'graph_id': pdfid,
                    # 'graph_png_id': graph_png_id,
                    'hyperlink': line[-1]
                }
                data_list.append(insert_data)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_level']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入kegg注释层级信息：%s、%s出错!" % (level_path, png_dir))
            else:
                self.bind_object.logger.info("导入kegg注释层级信息：%s、%s 成功!" % (level_path, png_dir))

    @report_check
    def add_annotation_kegg_pic(self, kegg_id, seq_type, anno_type, level_path, png_dir):
        """
        level_path: pathway_table.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(level_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(level_path))
        if not os.path.exists(png_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(png_dir))
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                pid = re.sub('path:', '', line[0])
                if os.path.exists(png_dir + '/' + line[0] + '.html.mark'):
                    with open(png_dir + '/' + line[0] + '.html.mark', 'r') as mark_f:
                        for line_mark in mark_f.readlines():
                            # print len(line_mark.strip("\n").split("\t"))
                            if len(line_mark.strip("\n").split("\t")) == 8:
                                [png, shape, bg_color, fg_color, coords, title, kos, href] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                            else:
                                continue

                            insert_data = {
                                'kegg_id': kegg_id,
                                'seq_type': seq_type,
                                'anno_type': anno_type,
                                'pathway_id': line[0],
                                'shape': shape,
                                'bg_colors': bg_color,
                                'fg_colors': fg_color,
                                'coords': coords,
                                'href': href,
                                'kos': kos,
                                'title': title
                            }

                            if bg_color != "" and len(bg_color.split(",")) > 0:
                                insert_data.update({'bg_type': len(bg_color.split(","))})
                            if fg_color != "" and len(fg_color.split(",")) > 0:
                                insert_data.update({'fg_type': len(fg_color.split(","))})
                            data_list.append(insert_data)
                else:
                    self.bind_object.logger.info("kegg 图片{} 不存在html标记!".format(line[0]))

        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_pic']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入kegg注释图片信息：%s、%s出错!" % (level_path, png_dir))
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))


    @report_check
    def add_annotation_kegg_table(self, kegg_id, seq_type, anno_type, table_path):
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(table_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(table_path))
        with open(table_path, 'rb') as r:
            data_list = []
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                insert_data = {
                    'kegg_id': kegg_id,
                    'seq_type': seq_type,
                    'anno_type': anno_type,
                    'transcript_id': line[0],
                    'ko_id': line[1],
                    'ko_name': line[2],
                    'hyperlink': line[3],
                    'paths': line[4],
                }
                data_list.append(insert_data)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_table']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入kegg注释table信息：%s出错!" % (table_path))
            else:
                self.bind_object.logger.info("导入kegg注释table信息：%s成功!" % (table_path))

    @report_check
    def add_annotation_kegg_categories_all(self, kegg_id, seq_type, anno_type, r_cate_path, n_cate_path):
        """
        r_cate_path:kegg_layer.xls(ref)
        n_cate_path:kegg_layer.xls(new)
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(r_cate_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(r_cate_path))
        if not os.path.exists(n_cate_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(n_cate_path))
        data_list = list()
        cate = {}
        with open(r_cate_path, 'r') as f, open(n_cate_path, "r") as n:
            lines = f.readlines()
            items = n.readlines()
            for line in lines:
                line = line.strip().split('\t')
                if line[0] == "Metabolism" and line[1] == "Global and overview maps":
                    pass
                else:
                    cate_d = line[0] + "|||" + line[1]
                    cate[cate_d] = line[3].split(";")
            for item in items:
                item = item.strip().split('\t')
                if item[0] == "Metabolism" and item[1] == "Global and overview maps":
                    pass
                else:
                    cate_d = item[0] + "|||" + item[1]
                    if cate_d not in cate:
                        cate[cate_d] = item[3].split(";")
                    else:
                        ids = item[3].split(";")
                        for q in ids:
                            if q not in cate[cate_d]:
                                cate[cate_d].append(q)
        for f in ["Metabolism", "Genetic Information Processing", "Environmental Information Processing", "Cellular Processes", "Organismal Systems", "Human Diseases", "Drug Development"]:
            for c in cate:
                ca = c.split("|||")
                first_category = ca[0]
                second_category = ca[1]
                if f == first_category:
                    num = len(cate[c])
                    seq_list = ";".join(cate[c])
                    data = [
                        ('kegg_id', kegg_id),
                        ('seq_type', seq_type),
                        ('anno_type', anno_type),
                        ('first_category', first_category),
                        ('second_category', second_category),
                        ('num', num),
                        ('seq_list', seq_list)
                    ]
                    data = SON(data)
                    data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_categories']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入kegg注释分类alls出错：%s, %s" % (r_cate_path, n_cate_path))
            else:
                self.bind_object.logger.info("导入kegg注释分类all成功：%s, %s" % (r_cate_path, n_cate_path))

    def get_pic(self, path, kos_path, png_path):
        """
        画通路图
        """
        fs = gridfs.GridFS(self.mongodb)
        pid = re.sub("map", "ko", path)
        with open("pathway.kgml", "w+") as k, open("pathway.png", "w+") as p:
            result = self.png_coll.find_one({"pathway_id": pid})
            if result:
                kgml_id = result['pathway_ko_kgml']
                png_id = result['pathway_map_png']
                k.write(fs.get(kgml_id).read())
                p.write(fs.get(png_id).read())
        cmd = "{} {} {} {} {} {} {}".format(self.r_path, self.map_path, path, kos_path, png_path, "pathway.kgml", "pathway.png")
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            print "{}画图出错".format(path)
            os.system("cp {} {}".format("pathway.png", png_path))

    @report_check
    def add_annotation_kegg_level_all(self, kegg_id, seq_type, anno_type, r_level_path, n_level_path):
        """
        r_level_path: pathway_table.xls(ref)
        n_level_path: pathway_table.xls(new)
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(r_level_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(r_level_path))
        if not os.path.exists(n_level_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(n_level_path))
        data_list = []
        path_def = {}
        r_path_list = {}
        n_path_list = {}
        fs = gridfs.GridFS(self.mongodb)
        with open(r_level_path, "rb") as r, open(n_level_path, "rb") as n:
            lines = r.readlines()
            items = n.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[1] == "Metabolism" and line[2] == "Global and overview maps":
                    pass
                else:
                    path = line[0] + "|||" + line[1] + "|||" + line[2] + "|||" + line[3]
                    sp = 'http://www.genome.jp/dbget-bin/show_pathway?' + line[0]
                    k_cols = line[7].split(sp)
                    k_ids = k_cols[1].split("%09yellow")
                    k_list = []
                    for k in k_ids:
                        if k.startswith('/'):
                            k_id = k.split('/')[1]
                            k_list.append(k_id)
                    seqlist = line[5].split(";")
                    path_def[line[0]] = path
                    r_path_list[line[0]] = []
                    r_path_list[line[0]].append(seqlist)
                    r_path_list[line[0]].append(k_list)
            for item in items[1:]:
                item = item.strip().split("\t")
                if item[1] == "Metabolism" and item[2] == "Global and overview maps":
                    pass
                else:
                    path = item[0] + "|||" + item[1] + "|||" + item[2] + "|||" + item[3]
                    sp = 'http://www.genome.jp/dbget-bin/show_pathway?' + item[0]
                    k_cols = item[7].split(sp)
                    k_ids = k_cols[1].split("%09green")
                    k_list = []
                    for k in k_ids:
                        if k.startswith('/'):
                            k_id = k.split('/')[1]
                            k_list.append(k_id)
                    seqlist = item[5].split(";")
                    n_path_list[item[0]] = []
                    n_path_list[item[0]].append(seqlist)
                    n_path_list[item[0]].append(k_list)
                    if item[0] not in path_def:
                        path_def[item[0]] = path
        for map_id in path_def:
            link = []
            r_kos, n_kos, b_kos = [], [], []
            ref, new, both = [], [], []
            paths = path_def[map_id].split("|||")
            first_category = paths[1]
            second_category = paths[2]
            pathway_definition = paths[3]
            try:
                seq_list = r_path_list[map_id][0]
                r_ko = r_path_list[map_id][1]
            except:
                seq_list = []
                r_ko = []
            try:
                for q in n_path_list[map_id][0]:
                    seq_list.append(q)
                n_ko = n_path_list[map_id][1]
            except:
                n_ko = []
            seq_list = list(set(seq_list))
            for r_c in r_ko:
                r = r_c.split("%09tomato")[0]
                if r_c not in n_ko:
                    r_id = r + '%09' + 'yellow'
                    link.append(r_id)
                    r_kos.append(r)
                else:
                    b_id = r + '%09' + 'tomato'
                    link.append(b_id)
                    b_kos.append(r)
            for n_c in n_ko:
                n = n_c.split("%09tomato")[0]
                if n_c not in r_ko:
                    n_id = n + '%09' + 'green'
                    link.append(n_id)
                    n_kos.append(n)
                else:
                    b_id = n + '%09' + 'tomato'
                    link.append(b_id)
                    b_kos.append(n)
            link = list(set(link))
            b_kos = list(set(b_kos))
            link = 'http://www.genome.jp/kegg-bin/show_pathway?' + map_id + '/' + '/'.join(link)
            png_path = os.getcwd() + '/' + map_id + ".png"
            pdf_path = os.getcwd() + '/' + map_id + ".pdf"
            kos_path = os.path.join(os.getcwd(), "KOs.txt")
            with open(kos_path, "w") as w:
                w.write("#KO\tbg\tfg\n")
                for k in n_kos:
                    w.write(k + "\t" + "#00CD00" + "\t" + "NA" + "\n")
                for k in r_kos:
                    w.write(k + "\t" + "#FFFF00" + "\t" + "NA" + "\n")
                for k in b_kos:
                    w.write(k + "\t" + "#FFFF00,#00CD00" + "\t" + "NA" + "\n")
            self.get_pic(map_id, kos_path, png_path)
            cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + png_path + ' ' + pdf_path
            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError:
                print '图片格式pdf转png出错'
            # pdfid = fs.put(open(pdf_path, 'rb'))
            # graph_png_id = fs.put(open(png_path, 'rb'))
            insert_data = {
                'kegg_id': kegg_id,
                'seq_type': seq_type,
                'anno_type': anno_type,
                'pathway_id': map_id,
                'first_category': first_category,
                'second_category': second_category,
                'pathway_definition': pathway_definition,
                'number_of_seqs': len(seq_list),
                'seq_list': ";".join(seq_list),
                # 'graph_id': pdfid,
                "hyperlink": link,
                # 'graph_png_id': graph_png_id
            }
            data_list.append(insert_data)
            # os.remove(pdf)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_level']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入kegg注释层级all信息出错：%s、%s" % (r_level_path, n_level_path))
            else:
                self.bind_object.logger.info("导入kegg注释层级all信息成功：%s、%s" % (level_path, png_dir))

    @report_check
    def add_annotation_query(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            print stat_id
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        if self.target_name:
            query_name = self.target_name
        else:
            query_name = 'AnnotationQuery' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': query_name,
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': '注释查询主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_query']
        query_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_query!")
        return query_id

    def get_ko2des(self):
        self.kegg_des = Config().SOFTWARE_DIR + "/database/Annotation/other/ko.des.txt"
        with open(self.kegg_des, 'r') as kegg_des_f:
            ko2des =  [(line.strip().split("\t")[0][3:], line.strip().split("\t")[1].split(";")[-1].strip()) for line in kegg_des_f.readlines()]
        return dict(ko2des)

    @report_check
    def add_annotation_query_denovo_detail2(self, query_id, query_path, seq_type, anno_type, exp_level):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                raise Exception('query_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(query_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(query_path))
        map2class = self.get_kegg_map2class()

        data_list = []
        ko2des= self.get_ko2des()

        tran2gene_dict = dict()
        gene_dict = dict()

        with open(query_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                is_gene = True
                if self.trans_isgene.has_key(line[0]) and self.trans_isgene[line[0]] == False:
                    is_gene = False
                if exp_level.lower() == "gene" and is_gene == False:
                    continue
                if line[0] in self.target_smallrna:
                    small_rnas = self.target_smallrna[line[0]]
                else:
                    small_rnas = []
                data = [
                    ('query_id', query_id),
                    #('anno_type', anno_type),
                    ('seq_type', seq_type),
                    ('transcript_id', line[0]),
                    ('small_rnas', ";".join(small_rnas)),
                    ('small_rnas_list', small_rnas),
                    ('gene_id', line[1]),
                    ('is_gene', is_gene),
                    ('gene_name', ""),
                ]

                try:
                    data.append(('length', line[2]))
                except:
                    data.append(('length', None))
                try:
                    des = line[12].split("[")[0]
                    des = line[12].split("(")[1]
                    data.append(('description', des))
                except:
                    data.append(('description', None))


                if is_gene:
                    data_annotation_add = list()
                    try:
                        data_annotation_add.append(('cog', line[3]))
                        data_annotation_add.append(('cog_description', line[5]))
                    except:
                        data_annotation_add.append(('cog_description', None))

                    # try:
                    #     data_annotation_add.append(('nog', line[4]))
                    #     data_annotation_add.append(('nog_description', line[6]))
                    # except:
                    #     data_annotation_add.append(('nog', None))
                    #     data_annotation_add.append(('nog_description', None))
                    try:
                        data_annotation_add.append(('ko_id', line[7]))
                        ko_des = ";".join([ko2des[ko] if ko in ko2des else "" for ko in line[8].split(";")])
                        data_annotation_add.append(('ko_description', ko_des))
                    except:
                        data_annotation_add.append(('ko_id', None))
                        data_annotation_add.append(('ko_description', None))
                    try:
                        data_annotation_add.append(('ko_name', line[8]))
                    except:
                        data_annotation_add.append(('ko_name', None))

                    try:
                        data_annotation_add.append(('pathways', line[9]))
                    except:
                        data_annotation_add.append(('pathways', None))

                    try:
                        pathways = line[9].split("; ")
                        pathways_class1 = [map2class[path.split("(")[0]][0] for path in pathways]
                        data_annotation_add.append(('pathways_class1', ";".join(pathways_class1)))
                    except:
                        data_annotation_add.append(('pathways_class1', ""))

                    try:
                        pathways = line[9].split("; ")
                        pathways_class2 = [map2class[path.split("(")[0]][1] for path in pathways]
                        data_annotation_add.append(('pathways_class2', ";".join(pathways_class2)))
                    except:
                        data_annotation_add.append(('pathways_class2', ""))

                    try:
                        data_annotation_add.append(('pfam', line[10]))
                    except:
                        data_annotation_add.append(('pfam', None))
                    try:
                        data_annotation_add.append(('go', line[11]))
                    except:
                        data_annotation_add.append(('go', None))
                    try:
                        data_annotation_add.append(('nr', line[12]))
                    except:
                        data_annotation_add.append(('nr', None))
                    try:
                        data_annotation_add.append(('swissprot', line[13]))
                    except:
                        data_annotation_add.append(('swissprot', None))
                    '''
                    try:
                        data.append(('enterz', line[14]))
                    except:
                        data.append(('enterz', None))
                    '''

                    gene_dict[line[1]] = data_annotation_add
                data_list.append(data)

        data_list_gene = list()
        for data in data_list:
            data_dict = dict(data)
            annotation_add = gene_dict[data_dict['gene_id']]
            data_list_gene.append(SON(data + annotation_add))

        try:
            collection = self.db['sg_annotation_query_detail']
            collection.insert_many(data_list_gene)
        except Exception, e:
            raise Exception("导入转录本注释统计信息：%s出错!" % (query_path))
        else:
            self.bind_object.logger.info("导入转录本注释统计信息：%s成功!" % (query_path))


    @report_check
    def add_annotation_query_denovo_detail(self, query_id, query_path, seq_type, anno_type, exp_level):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                raise Exception('query_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(query_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(query_path))
        map2class = self.get_kegg_map2class()

        data_list = []
        ko2des= self.get_ko2des()

        tran2gene_dict = dict()
        gene_dict = dict()
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                is_gene = True
                if self.trans_isgene.has_key(line[1]) and self.trans_isgene[line[1]] == False:
                    is_gene = False
                if exp_level.lower() == "gene" and is_gene == False:
                    continue
                if line[1] in self.target_smallrna:
                    small_rnas = self.target_smallrna[line[1]]
                else:
                    small_rnas = []
                data = [
                    ('query_id', query_id),
                    #('anno_type', anno_type),
                    ('seq_type', seq_type),
                    ('transcript_id', line[1]),
                    ('small_rnas', ";".join(small_rnas)),
                    ('small_rnas_list', small_rnas),
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

                if is_gene:
                    data_annotation_add = list()
                    try:
                        data_annotation_add.append(('cog', line[6]))
                        data_annotation_add.append(('cog_description', line[7]))
                    except:
                        data_annotation_add.append(('cog_description', None))
                    # try:
                    #     data_annotation_add.append(('nog', line[4]))
                    #     data_annotation_add.append(('nog_description', line[6]))
                    # except:
                    #     data_annotation_add.append(('nog', None))
                    #     data_annotation_add.append(('nog_description', None))
                    try:
                        data_annotation_add.append(('ko_id', line[8]))
                        ko_des = ";".join([ko2des[ko] if ko in ko2des else "" for ko in line[8].split(";")])
                        data_annotation_add.append(('ko_description', ko_des))
                    except:
                        data_annotation_add.append(('ko_id', None))
                        data_annotation_add.append(('ko_description', None))
                    try:
                        data_annotation_add.append(('ko_name', line[9]))
                    except:
                        data_annotation_add.append(('ko_name', None))

                    try:
                        data_annotation_add.append(('pathways', line[10]))
                    except:
                        data_annotation_add.append(('pathways', None))

                    try:
                        pathways = line[10].split("; ")
                        pathways_class1 = [map2class[path.split("(")[0]][0] for path in pathways]
                        data_annotation_add.append(('pathways_class1', ";".join(pathways_class1)))
                    except:
                        data_annotation_add.append(('pathways_class1', ""))

                    try:
                        pathways = line[10].split("; ")
                        pathways_class2 = [map2class[path.split("(")[0]][1] for path in pathways]
                        data_annotation_add.append(('pathways_class2', ";".join(pathways_class2)))
                    except:
                        data_annotation_add.append(('pathways_class2', ""))

                    try:
                        data_annotation_add.append(('pfam', line[11]))
                    except:
                        data_annotation_add.append(('pfam', None))
                    try:
                        data_annotation_add.append(('go', line[12]))
                    except:
                        data_annotation_add.append(('go', None))
                    try:
                        data_annotation_add.append(('nr', line[13]))
                    except:
                        data_annotation_add.append(('nr', None))
                    try:
                        data_annotation_add.append(('swissprot', line[14]))
                    except:
                        data_annotation_add.append(('swissprot', None))
                    try:
                        data_annotation_add.append(('enterz', line[15]))
                    except:
                        data_annotation_add.append(('enterz', None))
                    gene_dict[line[0]] = data_annotation_add
                data_list.append(data)
        data_list_gene = list()
        for data in data_list:
            data_dict = dict(data)
            annotation_add = gene_dict[data_dict['gene_id']]
            data_list_gene.append(SON(data + annotation_add))

        try:
            collection = self.db['sg_annotation_query_detail']
            collection.insert_many(data_list_gene)
        except Exception, e:
            raise Exception("导入转录本注释统计信息：%s出错!" % (query_path))
        else:
            self.bind_object.logger.info("导入转录本注释统计信息：%s成功!" % (query_path))


    @report_check
    def add_annotation_query_detail(self, query_id, query_path, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                raise Exception('query_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(query_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(query_path))
        data_list = []

        
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('query_id', query_id),
                    ('anno_type', anno_type),
                    ('transcript_id', line[0]),
                    ('gene_id', line[1]),
                ]
                try:
                    data.append(('gene_name', line[2]))
                except:
                    data.append(('gene_name', None))
                try:
                    data.append(('length', line[3]))
                except:
                    data.append(('length', None))
                try:
                    data.append(('cog', line[4]))
                    data.append(('cog_description', line[6]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                try:
                    data.append(('nog', line[5]))
                    data.append(('nog_description', line[7]))
                except:
                    data.append(('nog', None))
                    data.append(('nog_description', None))
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
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            collection.insert_many(data_list)
        except Exception, e:
            raise Exception("导入转录本注释统计信息：%s出错!" % (query_path))
        else:
            self.bind_object.logger.info("导入转录本注释统计信息：%s成功!" % (query_path))

    @report_check
    def add_annotation_gene_query_detail(self, query_id, query_path, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                raise Exception('query_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(query_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(query_path))
        data_list = []
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for j in range(1, len(lines)):
                line = lines[j].strip().split('\t')
                data = [
                    ('query_id', query_id),
                    ('gene_id', line[0]),
                    ('anno_type', anno_type),
                ]
                try:
                    data.append(('gene_name', line[1]))
                except:
                    data.append(('gene_name', None))
                try:
                    data.append(('cog', line[2]))
                    data.append(('cog_description', line[4]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                try:
                    data.append(('nog', line[3]))
                    data.append(('nog_description', line[5]))
                except:
                    data.append(('nog', None))
                    data.append(('nog_description', None))
                try:
                    data.append(('ko_id', line[6]))
                except:
                    data.append(('ko_id', None))
                try:
                    data.append(('ko_name', line[7]))
                except:
                    data.append(('ko_name', None))
                try:
                    data.append(('pathways', line[8]))
                except:
                    data.append(('pathways', None))
                try:
                    data.append(('pfam', line[9]))
                except:
                    data.append(('pfam', None))
                try:
                    data.append(('go', line[10]))
                except:
                    data.append(('go', None))
                try:
                    data.append(('nr', line[11]))
                except:
                    data.append(('nr', None))
                try:
                    data.append(('swissprot', line[12]))
                except:
                    data.append(('swissprot', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            collection.insert_many(data_list)
        except Exception, e:
            raise Exception("导入基因注释查询信息：%s出错!" % (query_path))
        else:
            self.bind_object.logger.info("导入基因注释统计信息：%s成功!" % (query_path))

    @report_check
    def add_annotation_gene_query_denovo_detail(self, query_id, query_path, seq_type, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                raise Exception('query_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(query_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(query_path))
        data_list = []
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for j in range(1, len(lines)):
                line = lines[j].strip().split('\t')
                if line[1] in self.target_trans:
                    mirna_list = self.target_smallrna[line[1]]
                else:
                    continue

                data = [
                    ('query_id', query_id),
                    ('mirna', ";".join(mirna_list)),
                    ('mirna_list', mirna_list),
                    ('gene_id', line[0]),
                    ('tran_id', line[1]),
                   # if $line[2] == "yes":
                   #     ('is_gene', True),
                    ('anno_type', anno_type),
                ]

                try:
                    data.append(('cog', line[1]))
                    data.append(('cog_description', line[3]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                try:
                    data.append(('nog', line[2]))
                    data.append(('nog_description', line[4]))
                except:
                    data.append(('nog', None))
                    data.append(('nog_description', None))
                try:
                    data.append(('ko_id', line[5]))
                except:
                    data.append(('ko_id', None))
                try:
                    data.append(('ko_name', line[6]))
                except:
                    data.append(('ko_name', None))
                try:
                    data.append(('pathways', line[7]))
                except:
                    data.append(('pathways', None))
                try:
                    data.append(('pfam', line[8]))
                except:
                    data.append(('pfam', None))
                try:
                    data.append(('go', line[9]))
                except:
                    data.append(('go', None))
                try:
                    data.append(('nr', line[10]))
                except:
                    data.append(('nr', None))
                try:
                    data.append(('swissprot', line[11]))
                except:
                    data.append(('swissprot', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            collection.insert_many(data_list)
        except Exception, e:
            raise Exception("导入基因注释查询信息：%s出错!" % (query_path))
        else:
            self.bind_object.logger.info("导入基因注释统计信息：%s成功!" % (query_path))

class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """
    def test_mongo(test):
        from mbio.workflows.small_rna.small_rna_test_api import SmallRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            "id": "majorbio_306492",
            #+ str(random.randint(1,10000)),
            #"id": "denovo_rna_v2",
            "project_sn": "35910_5fd96fe4d0589",
            #+ str(random.randint(1,10000)),
            "type": "workflow",
            "name": "AnnotationStat_origin_weihu_20201222_042038",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = SmallRnaTestApiWorkflow(wsheet)

        result_dir = "/mnt/lustre/users/sanger/app/database/Genome_DB_finish/plants/Vaccinium_corymbosum/v1.0_v1.0/Annotation_v2/annot_class"
        g2t2p = "/mnt/lustre/users/sanger/app/database/Genome_DB_finish/plants/Vaccinium_corymbosum/v1.0_v1.0/Annotation_v2/annot_class/tran2gene.txt"
        target = "/mnt/ilustre/users/sanger-dev/workspace/20200817/Smallrna_tsg_158369/TargetPredict/output/"

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("small_rna.target_annotation_weihu")
        params = {
            "nr_evalue": 1e-5,
            "nr_similarity": 0,
            "nr_identity": 0,
            "swissprot_evalue":1e-5,
            "swissprot_similarity": 0,
            "swissprot_identity": 0,
            "cog_evalue": 1e-5,
            "cog_similarity": 0,
            "cog_identity": 0,
            "kegg_evalue": 1e-5,
            "kegg_similarity": 0,
            "kegg_identity": 0,
            "pfam_evalue": 1e-5,
        }
        wf.test_api.species_name = "Homo_sapiens"

        params_target = {
            "miranda": "yes",
            "targetscan": "yes",
            "rnahybrid": "yes",
            "min_support": 2,
        }

        params_target = {
                "psrobot": "yes",
                "targetfinder": "no",
                "rnahybrid": "no",
                "min_support": str(2),
                'miranda_score': "160.0",
                'miranda_energy': "-20",
                'miranda_strict': "on",
                'rnahybrid_num': "100",
                'rnahybrid_energy': "-20",
                'rnahybrid_pvalue': "0.01",
                'ps_robot_score': "2.5",
                'targetfinder_score': "4"

        }
        target = "/mnt/ilustre/users/sanger-dev/workspace/20200824/Smallrna_tsg_218443/TargetPredict/output/"
        new_target_file = target + 'novol_target.xls'
        known_target_file = "/mnt/lustre/users/sanger/sg-users/liubinxu/weihu/target_predict_detail.xls"
        novol_target_file = None
        new_seq = '/mnt/ilustre/users/sanger-dev/workspace/20200824/Smallrna_tsg_218443/Srna/output/novel_mirna/novel_mature_seq.fa'
        known_seq = '/mnt/ilustre/users/sanger-dev/workspace/20200824/Smallrna_tsg_218443/Srna/KnownMirna/output/mature.fa'
        diff_summary = '/mnt/ilustre/users/sanger-dev/workspace/20200817/Smallrna_tsg_158369/Diffexp/output/DESeq2_diff_summary.xls'
        target_dir = "/mnt/ilustre/users/sanger-dev/workspace/20200817/Smallrna_tsg_158369/TargetPredict/output"

        wf.test_api.run(known_target_file, None, result_dir, g2t2p, params)
        # wf.test_api.import_target_detail(new_target_file, known_target_file, params_target, new_seq, known_seq, anno_type="origin", species_name="Homo_sapiens", target_dir=target_dir, version="v2")
        # wf.test_api.add_target_geneset(new_target_file, known_target_file, diff_summary)
        # wf.test_api.anno_type = 'latest'
        # wf.test_api.run(test_dir, trans2gene,  params)

if __name__ == '__main__':
    unittest.main()
