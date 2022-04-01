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
from mbio.api.database.denovo_rna_v2.api_assembly import ApiBase

class DenovoassAnnotation(ApiBase):
    def __init__(self, bind_object):
        super(DenovoassAnnotation, self).__init__(bind_object)
        self.result_dir = ''
        self.result_file = {}
        self.trans_gene = {}
        self.trans_isgene = {}
        self.task_id = self.bind_object.sheet.id
        self.anno_type = 'origin'
        #self._db_name = Config().MONGODB + '_ref_rna'

    def set_result_dir(self, annotation_mudule_dir):
        '''
        根据注释模块结果导入结果路径
        '''
        self.result_dir = annotation_mudule_dir
        self.bind_object.logger.info("导入**** {}".format(annotation_mudule_dir))
        self.result_file['new_stat_path'] = os.path.join(annotation_mudule_dir, "anno_stat/all_annotation_statistics.xls")
        self.result_file['new_venn_path'] = os.path.join(annotation_mudule_dir, "anno_stat/venn")
        blast_path = annotation_mudule_dir + "/anno_stat/blast"
        # for db in ["nr", "swissprot", "string", "kegg"]:
        for db in ["nr", "swissprot"]:
            db_gene_name = db + "gene_blast_path"
            db_trans_name = db + "trans_blast_path"
            self.result_file[db_gene_name] = blast_path + '/gene_' + db + '.xls'
            self.result_file[db_trans_name] = blast_path + '/' + db + '.xls'

        nr_path = annotation_mudule_dir + "/anno_stat/blast_nr_statistics"
        self.result_file['evalue_nr_path'] = nr_path + "/nr_evalue.xls"
        self.result_file['similar_nr_path'] = nr_path + "/nr_similar.xls"
        self.result_file['gene_evalue_nr_path'] = nr_path + "/gene_nr_evalue.xls"
        self.result_file['gene_similar_nr_path'] = nr_path + "/gene_nr_similar.xls"
        self.result_file['gene_species_nr_path'] = nr_path + "/gene_nr_species.xls"


        for key, value in self.result_file.items():
            if os.path.exists(value):
                pass
            else:
                self.bind_object.set_error('%s对应的结果文件%s 不存在，请检查', variables=(key, value), code="52001501")
        self.bind_object.logger.info("数据路径正确，文件完整 {}")

    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('assemble_id必须为ObjectId对象或其对应的字符串！', code="52001502")
        return object_id

    def get_trans2gene(self, trans2gene):
        trans2gene = str(trans2gene)
        if os.path.exists(trans2gene):
            pass
        else:
            self.bind_object.set_error('转录本基因对应的结果文件%s不存在，请检查', variables=(trans2gene), code="52001503")
        self.bind_object.logger.info("读入基因转录本对应关系文件 {}".format(trans2gene))
        with open(trans2gene, 'rb') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split("\t")
                self.trans_gene[line[0]] = line[1]
                if line[2] == "yes":
                    self.trans_isgene[line[0]] = True
                else:
                    self.trans_isgene[line[0]] = False
        self.bind_object.logger.info("读入基因转录本对应关系文件结束")


    def run(self, result_dir, trans2gene, params_dict, gene_exp=None, trans_exp = None, gene_species_count=None, trans_species_count=None, samples=None, taxon='Animals'):
        """
        new_anno_path: 新序列注释的结果文件夹
        pfam_path:转录本的pfam_domain
        merge_tran_output: 转录本的merge_annot tool输出结果路径
        merge_gene_output: 基因的merge_annot tool"
        """
        self.bind_object.logger.info("开始到表情数据路径为 {}".format(result_dir))
        self.set_result_dir(result_dir)
        self.get_trans2gene(trans2gene)

        task_id = self.task_id

        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        self.remove_table_by_main_record(main_table='sg_annotation_stat', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_stat_detail'], detail_table_key='stat_id')
        stat_id = self.add_annotation_stat(name=None, params=params, seq_type="new", database="nr,swissprot,pfam,cog,go,kegg", taxon=taxon, result_dir=self.result_dir)
        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['new_stat_path'], venn_path=self.result_file['new_venn_path'], gene_exp = gene_exp, trans_exp = trans_exp)
        self.update_db_record('sg_annotation_stat', stat_id, status="end", main_id=stat_id)


        params_select_nr = dict([(k,params_dict.get(k,None)) for k in ('nr_evalue', 'nr_similarity', 'nr_identity')])
        params_select_nr = json.dumps(params_select_nr, sort_keys=True, separators=(',', ':'))
        self.remove_table_by_main_record(main_table='sg_annotation_nr', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_nr_detail', 'sg_annotation_nr_pie'], detail_table_key='nr_id')
        nr_id = self.add_annotation_nr(name=None, params=params_select_nr, stat_id=stat_id, result_dir=self.result_dir, samples=samples)
        self.add_annotation_nr_pie(nr_id=nr_id, species_path=self.result_file['gene_species_nr_path'] ,  evalue_path=self.result_file['evalue_nr_path'], similar_path=self.result_file['similar_nr_path'], seq_type="new", anno_type="T")
        self.add_annotation_nr_pie(nr_id=nr_id, species_path=self.result_file['gene_species_nr_path'],  evalue_path=self.result_file['gene_evalue_nr_path'], similar_path=self.result_file['gene_similar_nr_path'], seq_type="new", anno_type="G")
        self.add_annotation_nr_species_detail(nr_id=nr_id, gene_species_count=gene_species_count, trans_species_count=trans_species_count, seq_type="new", anno_type="G")


        self.update_db_record('sg_annotation_nr', nr_id, status="end",  main_id=nr_id)


    @report_check
    def add_annotation_stat(self, name=None, params=None, seq_type=None, database=None, result_dir=None, taxon='Animals' ):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationStat_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            "taxonomy": taxon,
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
    def add_annotation_stat_detail(self, stat_id, stat_path, venn_path, gene_exp=None, trans_exp=None):
        """
        database: 进行统计的数据库
        stat_path: all_annotation_statistics.xls
        venn_path: venn图目录
        """
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52001504")
        if not os.path.exists(stat_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(stat_path), code="52001505")
        if not os.path.exists(venn_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(venn_path), code="52001506")


        trans_list = self.trans_gene.keys()
        genes_list = list(set(self.trans_gene.values()))
        trans2enum = {t:str(n) for n, t in enumerate(trans_list)}
        genes2enum = {g:str(n) for n, g in enumerate(genes_list)}
        exp_g_set = set()
        with open(gene_exp, 'r') as gene_exp_f:
            gene_exp_f.readline()
            for line in gene_exp_f:
                exps = line.strip().split("\t")[1:]
                if sum(map(float, exps)) != 0:
                    exp_g_set.add(genes2enum[line.strip().split("\t")[0]])
                else:
                    pass
        exp_t_set = set()
        with open(trans_exp, 'r') as trans_exp_f:
            trans_exp_f.readline()
            for line in trans_exp_f:
                exps = line.strip().split("\t")[1:]
                if sum(map(float, exps)) != 0:
                    exp_t_set.add(trans2enum[line.strip().split("\t")[0]])
                else:
                    pass

        gene_dict = dict()
        trans_dict = dict()

        data_list = []
        database_venn = {
            'NR': 'nr',
            'Swiss-Prot': 'swissprot',
            'Swiss-prot': 'swissprot',
            'Pfam': 'pfam',
            'KEGG': 'kegg',
            'GO': 'go',
            'COG': 'cog',
        }
        with open(stat_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                if line[0] == "Total":
                    line[3] = "1.0"
                    line[4] = "1.0"
                data = [
                    ('stat_id', stat_id),
                    ('type', line[0]),
                    ('transcript', int(line[1])),
                    ('gene', int(line[2])),
                    ('transcript_percent', round(float(line[3]), 4)),
                    ('gene_percent', round(float(line[4]), 4)),
                    ('anno_type', 'all')
                ]
                venn_list, gene_venn_list = None, None
                database = ["nr", "swissprot", "pfam", "kegg", "go", "string", "cog"]
                if database_venn.has_key(line[0]) and database_venn[line[0]] in database:
                    db = line[0]
                    venn = venn_path + "/" + database_venn[line[0]] + "_venn.txt"
                    gene_venn = venn_path + "/gene_" + database_venn[line[0]] + "_venn.txt"
                    if os.path.exists(venn) and os.path.exists(gene_venn):
                        with open(venn, "rb") as f:
                            venn_list = f.readline().strip('\n')
                            for line in f:
                                venn_list += ',{}'.format(trans2enum[line.strip('\n')])
                            trans_dict[db] = set(venn_list.split(","))
                        with open(gene_venn, "rb") as f:
                            gene_venn_list = f.readline().strip('\n')
                            for line in f:
                                gene_venn_list += ',{}'.format(genes2enum[line.strip('\n')])
                            gene_dict[db] = set(gene_venn_list.split(","))
                        data.append(("gene_list", gene_venn_list))
                        if len(gene_venn_list) + len(venn_list) < 15000000:
                            data.append(("transcript_list", venn_list))
                    else:
                        self.bind_object.set_error("%s对应的venn.txt文件不存在", variables=(line[0]), code="52001507")
                data = SON(data)
                data_list.append(data)

            anno_gene_list = set()
            anno_trans_list = set()
            for db in ['GO', 'KEGG', 'COG', 'NR', 'Swiss-Prot', 'Pfam']:
                trans_set = exp_t_set & set(trans_dict[db])
                gene_set = exp_g_set & set(gene_dict[db])
                data = [
                    ('stat_id', stat_id),
                    ('type', db),
                    ('transcript', len(trans_set)),
                    ('gene', len(gene_set)),
                    ('transcript_percent', round(float(len(trans_set))/float(len(exp_t_set)), 4)),
                    ('gene_percent', round(float(len(gene_set))/float(len(exp_g_set)), 4)),
                    ('anno_type', 'exp'),
                    ('gene_list', ",".join(list(gene_set))),
                    ('transcript_list', ",".join(list(trans_set)))
                ]
                anno_gene_list |= gene_set
                anno_trans_list |= trans_set
                data = SON(data)
                data_list.append(data)
            data = [
                ('stat_id', stat_id),
                ('type', 'Total_anno'),
                ('transcript', len(anno_trans_list)),
                ('gene', len(anno_gene_list)),
                ('transcript_percent', round(float(len(anno_trans_list))/float(len(exp_t_set)), 4)),
                ('gene_percent', round(float(len(anno_gene_list))/float(len(exp_g_set)), 4)),
                ('anno_type', 'exp')
            ]
            data = SON(data)
            data_list.append(data)
            data = [
                ('stat_id', stat_id),
                ('type', 'Total'),
                ('transcript', len(exp_t_set)),
                ('gene', len(exp_g_set)),
                ('transcript_percent', "1.0"),
                ('gene_percent', "1.0"),
                ('anno_type', 'exp')
            ]
            data = SON(data)
            data_list.append(data)

        try:
            collection = self.db['sg_annotation_stat_detail']
            for data in data_list:
                collection.insert_one(data)
            # collection.insert_many(data_list)
            # collection.insert_one(data)
        except Exception, e:
            self.bind_object.set_error("导入注释统计信息失败", code="52001508")
        '''
        try:
            collection = self.db['sg_annotation_stat_detail']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert(tmp_list)
            else:
                collection.insert(data_list)
        except Exception, e:
            raise Exception("导入注释统计信息失败:%s" % stat_path)
        else:
            self.bind_object.logger.info("导入注释统计信息成功：%s" % (stat_path))
        '''

    @report_check
    def add_stat_detail(self, old_stat_id, stat_id, nr_evalue, gene_nr_evalue, sw_evalue, gene_sw_evalue):
        """
        注释重运行时注释统计导表sg_annotation_stat_detail
        """
        if not isinstance(old_stat_id, ObjectId):
            if isinstance(old_stat_id, types.StringTypes):
                old_stat_id = ObjectId(old_stat_id)
            else:
                self.bind_object.set_error('old_stat_id必须为ObjectId对象或其对应的字符串！', code="52001509")
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52001510")
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
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                    collection.insert_many(data_list)
        except:
            self.bind_object.set_error("导入注释统计信息出错", code="52001511")
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
    def add_annotation_nr(self, name=None, params=None, stat_id=None, result_dir=None, samples=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52001512")
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
            'stat_id': stat_id,
            'samples': samples
        }
        collection = self.db['sg_annotation_nr']
        nr_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_nr!")
        return nr_id

    @report_check
    def add_annotation_nr_species_detail(self, nr_id, gene_species_count, trans_species_count, seq_type="new", anno_type="G"):
        if not isinstance(nr_id, ObjectId):
            if isinstance(nr_id, types.StringTypes):
                nr_id = ObjectId(nr_id)
            else:
                self.bind_object.set_error('nr_id必须为ObjectId对象或其对应的字符串！', code="52001513")
        if not os.path.exists(gene_species_count):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(evalue_path), code="52001514")
        if not os.path.exists(trans_species_count):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(similar_path), code="52001515")
        data_list = []
        with open(gene_species_count, 'r') as gene_spe_f:
            header = gene_spe_f.readline().strip().split("\t")
            for line in gene_spe_f:
                cols = line.strip().split("\t")
                data = zip(header, cols)
                data.extend([("anno_type", "G"), ("nr_id", nr_id)])
                data = SON(data)
                data_list.append(data)
        with open(trans_species_count, 'r') as trans_spe_f:
            header = trans_spe_f.readline().strip().split("\t")
            for line in trans_spe_f:
                cols = line.strip().split("\t")
                data = zip(header, cols)
                data.extend([("anno_type", "T"), ("nr_id", nr_id)])

                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_annotation_nr_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % ("sg_annotation_nr_detail", e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % ("sg_annotation_nr_detail"))

    @report_check
    def add_annotation_nr_pie(self, nr_id, species_path, evalue_path, similar_path, seq_type, anno_type):
        if not isinstance(nr_id, ObjectId):
            if isinstance(nr_id, types.StringTypes):
                nr_id = ObjectId(nr_id)
            else:
                self.bind_object.set_error('nr_id必须为ObjectId对象或其对应的字符串！', code="52001516")
        if not os.path.exists(evalue_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(evalue_path), code="52001517")
        if not os.path.exists(similar_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(similar_path), code="52001518")
        evalue, evalue_list,species_list_a, species_a, similar, similar_list = [], [], [], [], [], []
        data_list = []

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
            species_other = []
            species_other_list = []
            for line3 in lines3[16:]:
                line3 = line3.strip().split("\t")
                if anno_type == "T":
                    species = {"key": line3[0], "value": int(line3[1])}
                    try:
                        species_list = {"key": line3[0], "value": line3[5]}
                    except:
                        species_list = {"key": line3[0], "value": None}
                    other += int(line3[1])
                elif anno_type == "G":
                    try:
                        species_list = {"key": line3[0], "value": line3[6]}
                    except:
                        species_list = {"key": line3[0], "value": None}
                    other += int(line3[2])
                species_other.append(species)
                species_other_list.append(species_list)
            species = {"key": 'other', "value": other}
            species_list = {"key": 'other', "value": None}
            species_a.append(species)
            species_list_a.append(species_list)

            data = [
                ('nr_id', nr_id),
                #('seq_type', seq_type),
                ('anno_type', anno_type),
                ('species', species_a),
                ('species_list', species_list_a),
                ('species_other', species_other),
                ('species_other_list', species_other_list),
            ]

        data = SON(data)
        data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_nr_pie']
                if len(data_list) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入nr库注释作图信息evalue,similar：%s、%s出错!" , variables=(evalue_path, similar_path), code="52001519")
            else:
                self.bind_object.logger.info("导入nr库注释作图信息evalue,similar：%s、%s成功!" % (evalue_path, similar_path))


class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """
    def test_mongo(test):
        from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            "id": "denovo_assemble",
            #+ str(random.randint(1,10000)),
            #"id": "denovo_rna_v2",
            "project_sn": "denovo_assemble",
            #+ str(random.randint(1,10000)),
            "type": "workflow",
            "name": "denovo_rna_v2.denovo_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = DenovoTestApiWorkflow(wsheet)

        test_dir = '/mnt/ilustre/users/sanger-dev/workspace/20190618/DenovoAssemble_denovo_ass6741/output/annotation'
        trans2gene = '/mnt/ilustre/users/sanger-dev/workspace/20190618/DenovoAssemble_denovo_ass6741/output/assemble/Trinity.filter.gene_trans_map'
        samples = ['CL1', 'CL2', 'CL5', 'HFL3', 'HFL4', 'HFL6', 'HGL1', 'HGL3', 'HGL4']

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("denovo_rna_v2.denovoass_annotation")
        params = {
            "submit_location": "annotationstat",
            "task_id": "i-sanger_96683",
            "task_type": 2,
            "nr_evalue": 1e-3,
            "nr_similarity": 0,
            "nr_identity": 0,
            "swissprot_evalue":1e-3,
            "swissprot_similarity": 0,
            "swissprot_identity": 0,
            "cog_evalue": 1e-3,
            "cog_similarity": 0,
            "cog_identity": 0,
            "kegg_evalue": 1e-3,
            "kegg_similarity": 0,
            "kegg_identity": 0,
            "pfam_evalue": 1e-3,
        }
        wf.test_api.anno_type = 'origin'
        wf.test_api.run(test_dir, trans2gene,  params, gene_exp="/mnt/ilustre/users/sanger-dev/workspace/20190618/DenovoAssemble_denovo_ass6741/Quant/gene.tpm.matrix", trans_exp = "/mnt/ilustre/users/sanger-dev/workspace/20190618/DenovoAssemble_denovo_ass6741/Quant/transcript.tpm.matrix", gene_species_count= "/mnt/ilustre/users/sanger-dev/workspace/20190618/DenovoAssemble_denovo_ass6741/Quant/outpath", trans_species_count="/mnt/ilustre/users/sanger-dev/workspace/20190618/DenovoAssemble_denovo_ass6741/Quant/outpath", samples=samples, taxon="Plant")

if __name__ == '__main__':
    unittest.main()
