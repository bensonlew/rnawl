# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# from mainapp.config.db import get_mongo_client
from bson.objectid import ObjectId
import types
from bson import SON
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.meta import Meta
import random
import json
import re
import os


class Bacgenome(Meta):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(Bacgenome, self).__init__(self._bind_object)
        self._project_type = "bacgenome"

    def get_sample_genefile(self, task_id, sample, type="fnn"):
        if type not in ["fnn", "gff", "faa"]:
            raise Exception("输入type参数必须为fnn, gff, faa!")
        collection = self.db['gene_predict']
        result = collection.find_one({"task_id": task_id})
        file = result['file_path']
        path_list = []
        for i in sample.split(","):
            file_path = ''.join([file[0],'/', i,  file[1], i, file[2], type])
            path_list.append(file_path)
        path = ','.join(path_list)
        return path

    def get_genefile_bysample(self, task_id, sample, type="faa", predict = "gene"):  # add by zhujuan 20180409
        if type not in ["fnn", "gff", "faa"]:
            raise Exception("输入type参数必须为fnn, gff, faa!")
        if predict not in ["gene", "trna", "rrna"]:
            raise Exception("输入predict必须为gene, trna, rrna!")
        collection = self.db[predict + '_predict']
        result = collection.find_one({"task_id": task_id})
        result_path = result['file_path']
        sample_gene_path = {}
        for i in sample.split(","):
            file_path = ''.join([result_path[0], i, result_path[1], i, result_path[2], type])
            sample_gene_path[i] = file_path
        return sample_gene_path

    def get_assemble_bysample(self, task_id, sample, type, seq_type):
        if type not in ["uncomplete", "complete"]:
            raise Exception("输入type参数必须为uncomplete,complete!")
        collection = self.db['gene_predict']
        result = collection.find_one({"task_id": task_id})
        result_path = result['file_path']
        sample_gene_path = {}
        if type in ['uncomplete']:
            if seq_type in ['scaffold']:
                for i in sample.split(";"):
                    i = i.split(":")[0]
                    file_path = ''.join([result_path[0], i, "/assembly_predict/assembly/", i, '_scaf.fna'])
                    sample_gene_path[i] = file_path
            elif seq_type in ['contig']:
                for i in sample.split(";"):
                    i = i.split(":")[0]
                    file_path = ''.join([result_path[0], i, "/assembly_predict/assembly/", i, '_ctg.fna'])
                    sample_gene_path[i] = file_path
        elif type in ['complete']:
            for i in sample.split(";"):
                i = i.split(":")[0]
                if re.search(r'://', result_path[0]):
                    file_path = ''.join([result_path[0], i, "/assembly_predict/assembly/seq_dir/"])
                else:
                    file_path = ''.join([result_path[0], i, "/assembly_predict/assembly/"])
                sample_gene_path[i] = file_path
        return sample_gene_path



    def get_anno_summary(self, task_id, sample): # add by ysh 20190423
        collection = self.db['gene_predict']
        result = collection.find_one({"task_id": task_id})
        result_path = result['file_path']
        result_dir = result_path[0]
        sample_anno = {}
        for i in sample.split(","):
            file_path = ''.join([result_path[0], i, "/annotation/Summary/", i, "_anno_summary.xls"])
            sample_anno[i] = file_path
        return sample_anno

    def get_orthmcl_file(self, task_id, homology_id):
        if isinstance(homology_id, types.StringTypes):
            homology_id = ObjectId(homology_id)
        elif isinstance(homology_id, ObjectId):
            homology_id = homology_id
        else:
            raise Exception("输入homology_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['compare_homology']
        result = collection.find_one({"task_id": task_id, "_id": homology_id})
        if not result:
            raise Exception('compare_homology没有相应homology_id:{}对应的信息!'.format(homology_id))
        stat_file = result['stat_file']
        return stat_file

    def get_projectsn(self, task_id):
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        project_sn = result["project_sn"]
        return project_sn

    def get_specimen_name(self, task_id, specimen_id):
        main_collection = self.db['datastat']
        main = main_collection.find_one({"task_id": task_id})
        main_id = main["_id"]
        spe_collection = self.db['datastat_specimen']
        result = spe_collection.find_one({"datastat_id": main_id, "init_name": specimen_id})
        specimen_name = result['anay_name']
        return specimen_name

    def get_pre_file(self, task_id, specimen_id, location=None, return_color_version=False):
        main_collection = self.db['circos_table']
        main = main_collection.find_one({"task_id": task_id})
        main_id = main["_id"]
        detail_collection = self.db['circos_table_detail']
        if location:
            result = detail_collection.find_one(
                {"circos_table_id": main_id, "specimen_id": specimen_id, "location": location})
        else:
            result = detail_collection.find_one({"circos_table_id": main_id, "specimen_id": specimen_id})
        pre_file = result['file_path']
        if return_color_version:
            if 'color_version' in main:
                return pre_file, main['color_version']
            else:
                return pre_file, ''
        else:
            return pre_file


    def get_xml_file(self,task_id,specimen_id,genome_type):
        main_collection = self.db['cgview']
        #main = main_collection.find_one({"task_id": task_id})
        #main_id = main["_id"]
        if genome_type == 'scaffold':
            result = main_collection.find_one({"task_id": task_id, "specimen_id": specimen_id,"genome_type":"scaffold"})
        else:
            result = main_collection.find_one({"task_id": task_id, "specimen_id": specimen_id, "genome_type": genome_type})
        xml_file = result['xml_path']
        return xml_file

    def get_gbk(self,gbk_detail_id,gbk_id):
        main_collection = self.db['gbk_detail']
        result = main_collection.find_one({"_id": gbk_detail_id,"gbk_id":gbk_id})
        gbk_file = result['params']
        return gbk_file

    def get_gbk_path(self,gbk_detail_id):
        main_collection = self.db['gbk_detail']
        result = main_collection.find_one({"_id": gbk_detail_id})
        gbk_path = result['path']
        return gbk_path

    def get_one_common_info(self,table,search_dic,ref=False):
        if ref:
            collection = self.ref_db[table]
        else:
            collection = self.db[table]
        result = collection.find_one(search_dic)
        if result:
            return result
        else:
            return False

    def get_all_common_info(self,table,search_dic,ref=False):
        if ref:
            collection = self.ref_db[table]
        else:
            collection = self.db[table]
        result = collection.find(search_dic)
        if result:
            return result
        else:
            return False

    def get_genome_info(self, task_id, genome, sample):
        fasta_path = self.db["complete_stat"].find_one({"task_id": task_id, "samp": sample, "genome": genome, "status":"end"})["seq_path"]
        return fasta_path

    def get_genome_seq(self, task_id, soft_type, sample):
        main_collection = self.db["draft"]
        main_id = main_collection.find_one({"task_id": task_id})['_id']
        fasta_path = self.db["draft_stat_detail"].find_one({"draft_id": ObjectId(main_id), "sof_type": soft_type, "samp": sample})["seq_path"]
        return fasta_path

    def get_pe_reads(self, task_id, sample):
        main_collection = self.db["sample_info"]
        main_id = main_collection.find_one({"task_id": task_id})["_id"]
        detail_collection = self.db["qc_stat_detail"]
        result = detail_collection.find_one({"info_id": main_id, "lib_type": "PE", "samp": sample})
        if result:
            return result["fq1"],result["fq2"]
        else:
            raise Exception("get_pe_reads error, check db:sample_info, qc_stat_detail, main_id:%s" % main_id)

    def get_gap_fill(self, task_id, sample, sof_type):
        main_coll = self.db["draft"]
        main_id = main_coll.find_one({"task_id": task_id})['main_id']
        fasta_path = self.db["draft_stat_detail"].find_one({"draft_id": ObjectId(main_id),"samp": sample, "sof_type": sof_type})['seq_path']
        return fasta_path

    def get_gap_fill2(self, main_id):
        main_coll = self.db["gap_fill"]
        result = main_coll.find_one({"main_id": ObjectId(main_id)})
        seq_path = result["seq_path"]
        return seq_path

    def get_genome_fill(self, task_id, sample, genome):
        main_coll = self.db["genome"]
        resul = main_coll.find_one({"task_id": task_id, "sample": sample, "genome": genome})
        genome_id = ObjectId(resul['_id'])
        gap_id = resul['gap_id']
        main_coll2 = self.db["gap_fill"]
        result = main_coll2.find_one({"main_id": ObjectId(gap_id)})
        seq_path = result["seq_path"]
        return seq_path, genome_id

    def get_gap_fill_path(self, genome_id):
        genome_id = ObjectId(genome_id)
        main_collection = self.db["genome"]
        result = main_collection.find_one({"_id": genome_id})
        ids = result["gap_detail_ids"]
        new_ids = []
        for id in ids:
            new_ids.append(str(id))
        return ",".join(new_ids)

    def update_mongo(self,db_name,search, change):
        db = self.db[db_name]
        ret = db.find_one(search)
        if ret:
            db.update({"_id":ret["_id"]},{"$set":change})
            return str(ret["_id"])
        else:
            return 'not find'

















