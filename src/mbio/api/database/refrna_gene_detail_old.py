#!/usr/bin/python
# -*- coding: utf-8 -*-
import sqlite3
import subprocess
import os,re,pickle,json
from bs4 import BeautifulSoup
from pymongo import MongoClient
from bson.objectid import ObjectId
import types
from types import StringTypes
import re
import json, time
import pandas as pd
import numpy as np
import datetime, os
from bson.son import SON
from collections import Counter
import glob
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import math
import subprocess
from bson.json_util import loads


class RefrnaGeneDetail(Base):
    def __init__(self, bind_object):
        super(RefrnaGeneDetail, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_ref_rna'
        db = Config().mongo_client[Config().MONGODB + "_ref_rna"]

    # def __init__(self):
    #     self._db_name = Config().MONGODB + '_ref_rna'
    #     db = Config().mongo_client[Config().MONGODB + "_ref_rna"]

    # @report_check
    def add_express_diff_class_code_detail(self, class_code, class_code_id, species=None):
        if not isinstance(class_code_id, ObjectId):
            if isinstance(class_code_id, types.StringTypes):
                express_id = ObjectId(class_code_id)
            else:
                raise Exception('class_code_id必须为ObjectId对象或其对应的字符串！')
        db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        data_list = []
        with open(class_code, 'r+') as f1:
            f1.readline()
            for lines in f1:
                line = lines.strip().split("\t")
                data = [
                    ('assembly_trans_id', line[0]),
                    ('assembly_gene_id', line[1]),
                    ('class_code', line[2]),
                    ('gene_name', line[3]),
                    ('class_code_id', ObjectId(class_code_id)),
                    ("type", "express_diff")
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = db["sg_express_class_code_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            print("导入%s表出错:%s" % (class_code, e))
        else:
            print("导入%s表成功！" % (class_code))

    def biomart(self,biomart_path, species_name=None):
        start = time.time()
        ss = 0
        biomart = dict()
        def check(_id):
            if not _id:
                return '-'
            else:
                return _id
        with open(biomart_path, 'r+') as f1:
            for lines in f1:
                ss += 1
                line = lines.strip().split("\t")
                gene_id = check(line[0])
                trans_id = check(line[1])
                gene_name = check(line[2])
                chromosome = check(line[9])
                gene_type = check(line[17])
                desc = check(line[8])
                strand_tmp = check(line[12])
                if strand_tmp == "1":
                    strand = "+"
                elif strand_tmp == "-1":
                    strand = "-"
                else:
                    strand = "."
                start = check(str(line[10]))
                end = check(str(line[11]))
                # location = "{}:{}-{}".format(check(line[12]), check(str(line[10])), check(str(line[11])))
                pep_id = check(line[6])
                if gene_id:
                    if gene_id not in biomart.keys():
                        biomart[gene_id] = {}
                        biomart[gene_id]= {"trans_id":[trans_id],"gene_name":[gene_name],"chromosome":[chromosome],
                                          "gene_type":[gene_type],"description":[desc],"strand":[strand],"pep_id":[pep_id],
                                           "start":[start],"end":[end]}
                    else:
                        biomart[gene_id]['trans_id'].append(trans_id)
                        biomart[gene_id]['gene_name'].append(gene_name)
                        biomart[gene_id]['chromosome'].append(chromosome)
                        biomart[gene_id]['gene_type'].append(gene_type)
                        biomart[gene_id]['description'].append(desc)
                        #biomart[gene_id]['location'].append(location)
                        biomart[gene_id]['pep_id'].append(pep_id)
                        biomart[gene_id]['strand'].append(strand)
                        biomart[gene_id]['start'].append(start)
                        biomart[gene_id]['end'].append(end)
                        # print ("同一个gene_id对应的trans_id{} 含有的pep_id为{}".format(",".join(biomart[gene_id]['trans_id']),",".join(biomart[gene_id]['pep_id'])))
        if not biomart:
            print "没有生成biomart信息"
        print "biomart共统计出{}行信息".format(str(ss))
        return biomart

    def entrez_id(self,entrez_path,entrez_db_path):
        """此函数暂时没有用到"""
        start = time.time()
        c=sqlite3.connect(entrez_db_path)
        conn = c.cursor()
        drop = lambda cur:cur.execute('DROP TABLE FILE')
        def create(cur):
            try:
                conn.execute("""CREATE TABLE gene2ensembl(id INT PRIMARY KEY NOT NULL, tax_id CHAR, geneid CHAR, ensem_gene_id CHAR, entrez_gene_id CHAR,
                ensem_trans_id CHAR, entrez_pep_id CHAR, ensem_pep_id CHAR);""")
            except sqlite3.OperationalError,e:
                drop(cur)
                create(cur)
        create(conn)
        ll=0
        with open(entrez_path, 'r+') as f1:
            f1.readline()
            for lines in f1:
                ll += 1
                line = lines.strip().split("\t")
                conn.execute("INSERT INTO gene2ensembl (id,tax_id,geneid,ensem_gene_id,entrez_gene_id,ensem_trans_id,entrez_pep_id,ensem_pep_id) VALUES (?,?,?,?,?,?,?,?)",
                             (ll,line[0],line[1],line[2],line[3],line[4],line[5],line[6]))
                if ll%1000 ==0:
                    print ll
        c.commit()
        c.close()
        end = time.time()
        duration = end - start
        m, s = divmod(duration, 60)
        h, m = divmod(m, 60)
        print('整个entrez_id建db的运行的时间为{}h:{}m:{}s'.format(h, m, s))
        print 'end'

    def insert_entrez_id(self,entrez_path,task_id,project_sn):
        db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        count = 0
        with open(entrez_path,'r+') as f1:
            f1.readline()
            data_list=[]
            #tax_id,geneid,ensem_gene_id,entrez_gene_id,ensem_trans_id,entrez_pep_id,ensem_pep_id
            for lines in f1:
                count+=1
                line = lines.strip().split("\t")
                data=[
                    ('taxid',line[0]),
                    ('geneid',line[1]),
                    ('ensem_gene_id',line[2]),
                    ('entrez_gene_id',line[3]),
                    ('ensem_trans_id',line[4]),
                    ('entrez_pep_id',line[5]),
                    ('ensem_pep_id',line[6])
                ]
                data=SON(data)
                data_list.append(data)
                if count%1000==0:
                    print count
            try:
                collection = db['sg_entrez_id']
                collection.insert_many(data_list)
            except Exception, e:
                print ("导入entrez_id表出错:%s" % e)
            else:
                print ("导入entrez_id表成功！")

    def get_cds_seq(self,cds_path):
        """转录本的cds信息已经提取出来了"""
        start = time.time()
        trans = dict()
        j = 0
        with open(cds_path, 'r+') as f1:
            for lines in f1:
                line = lines.strip()
                if re.search(r'>', line):
                    j += 1
                    cds_m = re.search(r'\>(\w+\.\w+).\w+.+gene:(\w+)', line)
                    trans_m  = re.search(r'\>(\w+).\w+',line)
                    if cds_m:
                        if j > 1:
                            if trans_id not in trans.keys():
                                trans[trans_id]={"name":cds_id,"sequence":cds_sequence,"sequence_length":len(cds_sequence)}
                        cds_id = cds_m.group(1)
                        trans_id = trans_m.group(1)
                        cds_sequence = ''
                else:
                    cds_sequence += line
            if trans_id not in trans.keys():
                trans[trans_id] = {"name":cds_id,"sequence":cds_sequence,"sequence_length":len(cds_sequence)}
        if not trans:
            print '提取cds序列信息为空'
        print "共统计出{}行信息".format(str(j))
        end = time.time()
        duration = end - start
        m, s = divmod(duration, 60)
        h, m = divmod(m, 60)
        print('cds提取运行的时间为{}h:{}m:{}s'.format(h, m, s))
        return trans

    def get_pep_seq(self,pep_path):
        start = time.time()
        trans = dict()
        j = 0
        with open(pep_path, 'r+') as f1:
            for lines in f1:
                line = lines.strip()
                if re.search(r'>', line):
                    j += 1
                    if j > 1:
                        if trans_id not in trans.keys():
                            trans[trans_id] = {"name":pep_id, "sequence":pep_sequence, "sequence_length":len(pep_sequence)}
                    pep_m = re.search(r'\>(\w+\.\w+)', line)
                    if pep_m:
                        pep_id = pep_m.group(1)
                    trans_m = re.search(r'transcript:(\w+)',line)
                    if trans_m:
                        trans_id = trans_m.group(1)
                    pep_sequence = ''
                else:
                    pep_sequence += line
            if trans_id not in trans.keys():
                trans[trans_id] = {"name":pep_id, "sequence":pep_sequence, "sequence_length":len(pep_sequence)}
        if not trans:
            print '提取cds序列信息为空'
        print "共统计出{}行信息".format(str(j))
        end = time.time()
        duration = end - start
        m, s = divmod(duration, 60)
        h, m = divmod(m, 60)
        print('pep提取运行的时间为{}h:{}m:{}s'.format(h, m, s))
        return trans

    def query_entrezid(self,ensem_gene_id=None):
        """只能对参考基因组上的ensembl id查询"""
        db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        collection = db["sg_entrez_id"]
        data = collection.find_one({"ensem_gene_id":ensem_gene_id})
        if data:
            return int(data['geneid'])
        else:
            return '-'

    def new_gene_location(self,new_gene_path):
        new_gene_info = {}
        with open(new_gene_path,'r+') as f1:
            f1.readline()
            for lines in f1:
                line = lines.strip().split("\t")
                if line[0] not in new_gene_info.keys():
                    new_gene_info[line[0]]={"chr":line[1], "strand":line[2], "start":line[3], "end":line[4]}
                else:
                    raise Exception("在提取新基因的位置信息中,基因{}有重复的位置信息".format(line[0]))
        return new_gene_info

    def extract_class_code(self,class_code):
        start = time.time()
        with open(class_code,'r+') as f1:
            f1.readline()
            extract_class_code = {}
            for lines in f1:
                line =lines.strip().split("\t")
                if line[1] not in extract_class_code.keys():
                    extract_class_code[line[1]] = [line[0]]
                else:
                    extract_class_code[line[1]].append(line[0])
            if extract_class_code:
                end = time.time()
                duration = end - start
                m, s = divmod(duration, 60)
                h, m = divmod(m, 60)
                print('class_code提取运行的时间为{}h:{}m:{}s'.format(h, m, s))
                return extract_class_code
            else:
                print '没有提取出基因对应的转录本信息'
                raise Exception("error")

    def add_gene_detail_class_code_detail(self,class_code, class_code_id=None, biomart_path=None,new_gene_location_path=None,cds_path=None,pep_path=None,transcript_path=None,gene_path=None,species=None):
        # if not isinstance(class_code_id, ObjectId):
        #     if isinstance(class_code_id, types.StringTypes):
        #         express_id = ObjectId(class_code_id)
        #     else:
        #         raise Exception('class_code_id必须为ObjectId对象或其对应的字符串！')
        db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        collection = db["sg_express_class_code_detail"]
        start = time.time()
        db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        data_list = []
        class_code_gene_trans = self.extract_class_code(class_code)
        new_gene_location_info = self.new_gene_location(new_gene_location_path)
        trans_cds_info = self.get_cds_seq(cds_path)  #转录本和cds对应关系
        trans_pep_info = self.get_pep_seq(pep_path)  #转录本和pep对应关系
        biomart_data = self.biomart(biomart_path)
        trans_sequence = self.get_transcript_seq(transcript_path=transcript_path)
        gene_sequence = self.get_transcript_seq(transcript_path=gene_path,_type="gene")
        gene_id_list = []
        line_id = 0
        with open(class_code, 'r+') as f1:
            f1.readline()
            for lines in f1:
                line_id +=1
                # if line_id <=2:
                if line_id%1000 ==0:
                    print line_id
                line = lines.strip().split("\t")
                if line[1] in gene_id_list:
                    continue
                else:
                    gene_id_list.append(line[1])
                    transcript_info = {}
                    entrez_id = self.query_entrezid(line[1])
                    trans_list = class_code_gene_trans[line[1]]
                    transcript_num = int(len(trans_list))
                    if not trans_list:
                        print '{}基因没有对应的转录本信息'.format(line[1])
                        raise Exception("error!")
                    else:
                        transcript = ",".join(trans_list)
                    new_trans_list = []
                    for i in trans_list:
                        m_ = re.search(r'\.|\_',i)
                        if m_:
                            continue
                        else:
                            new_trans_list.append(i)
                    for trans_ll in new_trans_list:
                        if trans_ll not in transcript_info.keys():
                            transcript_info[trans_ll]={}
                            count_none =0
                            if trans_ll in trans_cds_info.keys():
                                transcript_info[trans_ll]["cds"]=trans_cds_info[trans_ll]
                            else:
                                count_none+=1
                            if trans_ll in trans_pep_info.keys():
                                transcript_info[trans_ll]["pep"] = trans_pep_info[trans_ll]
                            else:
                                count_none+=1
                            if trans_ll in trans_sequence.keys():
                                transcript_info[trans_ll]["length"] = trans_sequence[trans_ll]
                            else:
                                # transcript_info.pop(trans_ll)
                                transcript_info[trans_ll]["length"] = 0
                            if count_none ==2:
                                transcript_info.pop(trans_ll)
                        else:
                            print "{}没有找到对应的transcript_info信息".format(trans_ll)
                            raise Exception("error")
                    if line[1] in biomart_data.keys():
                        desc = biomart_data[line[1]]["description"][0]
                        # location = biomart_data[line[1]]["location"][0].split(":")
                        strand = biomart_data[line[1]]["strand"][0]
                        start = int(biomart_data[line[1]]["start"][0])
                        end = int(biomart_data[line[1]]["end"][0])
                        # start,end =location[1].split('-')
                        # start = int(start)
                        # end = int(end)
                        gene_name = biomart_data[line[1]]['gene_name'][0]
                        gene_type = biomart_data[line[1]]['gene_type'][0]
                        chrom = biomart_data[line[1]]['chromosome'][0]
                    elif line[1] in new_gene_location_info.keys():
                        desc='-'
                        strand = new_gene_location_info[line[1]]["strand"]
                        start = new_gene_location_info[line[1]]["start"]
                        end = new_gene_location_info[line[1]]["end"]
                        gene_name = '-'
                        chrom = new_gene_location_info[line[1]]["chr"]
                        gene_type = "-"
                    else:
                        continue
                        # raise Exception("{}既不属于新基因也不在biomart中".format(line[1]))
                    if line[1] in gene_sequence.keys():
                        gene_seq = gene_sequence[line[1]]
                        gene_len = int(end) - int(start)
                    elif line[3] in gene_sequence.keys():
                        gene_seq = gene_sequence[line[3]]
                        gene_len = int(end) - int(start)
                    else:
                        gene_seq = '-'
                        gene_len = int(end) - int(start)
                    data =[
                        ("type","gene_detail"),
                        ("class_code_id",class_code_id),
                        ("gene_id",line[1]),
                        ("entrez_id",entrez_id),
                        ("description",desc),
                        ("strand",strand),
                        ("start",start),
                        ("location","{}-{}".format(str(start),str(end))),
                        ("end",end),
                        ("chrom",chrom),
                        ("gene_name",gene_name),
                        ("transcript",transcript),
                        ("transcript_number",transcript_num),
                        ("gene_type",gene_type),
                        ("gene_sequence",gene_seq),
                        ("gene_length",gene_len),
                        ("gene_ncbi", "https://www.ncbi.nlm.nih.gov/gquery/?term={}".format(line[1])),
                        ("gene_ensembl", "http://www.ensembl.org/{}/Gene/Summary?g={}".format(species, line[1])),
                    ]
                    if not new_trans_list:
                        data.append(("trans_info",None))
                    elif not transcript_info:
                        data.append(("trans_info",None))
                    else:
                        data.append(("trans_info",transcript_info))
                    # print data
                    data=SON(data)
                    data_list.append(data)
                # else:
                #         break

        print "共统计出{}基因".format(str(len(gene_id_list)))
        # end = time.time()
        # duration = end - start
        # m, s = divmod(duration, 60)
        # h, m = divmod(m, 60)
        # print('整个程序运行的时间为{}h:{}m:{}s'.format(h, m, s))
        # print data_list
        try:
            collection = db["sg_express_class_code_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            print("导入%s表出错:%s" % (class_code, e))
        else:
            print("导入%s表成功！" % (class_code))

    def add_class_code(self,assembly_method, name=None, major_gene_detail=False, major_express_diff=False,class_code_path=None, new_gene_location_path=None,biomart_path=None,species=None,cds_path=None,pep_path=None,transcript_path=None,gene_path=None,):
        db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        task_id = self.bind_obj.sheet.task_id
        project_sn = self.bind_obj.sheet.project_sn
        # task_id = 'demo2'
        # project_sn = '10001368'
        data = [
            ('task_id', task_id),
            ('project_sn', project_sn),
            ('assembly_method', assembly_method),
            ('desc', 'class_code信息'),
            ('created_ts', datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')),
            ('status', 'end'),
            ('name', name if name else 'Classcode_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")))
        ]
        if not os.path.exists(transcript_path):
            raise Exception("{}文件不存在!".format(transcript_path))
        else:
            transcript_length =0
            with open(transcript_path,'r+') as f1:
                for lines in f1:
                    if re.search(r'\>',lines):
                        continue
                    else:
                        transcript_length += int(len(lines.strip()))
            data.append(('transcripts_total_length',transcript_length))
        collection = db["sg_express_class_code"]
        print data
        try:
            class_code_id = collection.insert_one(SON(data)).inserted_id
        except Exception, e:
            print("导入class_code主表信息错误" )
        else:
            print("导入class_code主表信息成功!")
        if major_gene_detail:
            self.add_gene_detail_class_code_detail(class_code=class_code_path, class_code_id=class_code_id, new_gene_location_path=new_gene_location_path,biomart_path=biomart_path,cds_path=cds_path, pep_path=pep_path,transcript_path=transcript_path,
                                                   gene_path=gene_path, species=species)
        if major_express_diff:
            self.add_express_diff_class_code_detail(class_code_path, class_code_id, species=species)


    def insert_sequence_id(self,cds_info,species_name,sequence_type=None,task_id=None,project_sn=None):
        db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        count = 0
        data_list = []
        for keys,values in cds_info.items():
            data = [
                ("species_name",species_name),
                ("gene_id",keys),
                ("type",sequence_type)
            ]
            cds_id = []
            cds_sequence=[]
            for keys1,values1 in values.items():
                if keys1 not in cds_id:
                    cds_id.append(keys1)
                else:
                    raise Exception("gene_id{}有重复的cds_id{}".format(keys,keys1))
                if values1 not in cds_sequence:
                    cds_sequence.append(values1)
            length = [str(len(i))for i in cds_sequence]
            if sequence_type == 'cds':
                data +=[
                    ("cds_id",",".join(cds_id)),
                    ("cds_sequence",",".join(cds_sequence)),
                    ("cds_length",",".join(length))
                ]
            else:
                data+=[
                    ("pep_id",",".join(cds_id)),
                    ("pep_sequence",",".join(cds_sequence)),
                    ("pep_length",length)
                ]
            data=SON(data)
            data_list.append(data)
        try:
            collection = db["sg_cds_pep_sequence"]
            collection.insert_many(data_list)
        except Exception, e:
            print ("导入%s信息出错:%s" % (sequence_type,e))
        else:
            print ("导入{}信息成功！".format(sequence_type))

    def get_transcript_seq(self,transcript_path,_type="transcript",assembly_method=None):
        """提取gene和transcript序列信息是同一个函数"""
        start = time.time()
        seq = dict()
        j = 0
        count = 0
        with open(transcript_path, 'r+') as f1:
            for lines in f1:
                line = lines.strip()
                if re.search(r'>', line):
                    j += 1
                    if line.startswith('>MSTRG'):
                        m_ = re.search(r'\>(\w+\.\w+).+gene=(\w+\.\w+)',line)
                        if not m_:
                            m_=re.search(r'\>(\w+\.\w+).+gene=(\w+)',line)
                    # elif line.startswith('>TCONS'):
                    #     m_ = re.search(r'\>(\w+\_\w+).+gene=(\w+\_\w+)',line)
                    else:
                        m_ = re.search(r'\>(\w+).+gene=(\w+)', line)
                    if m_:
                        if _type == "transcript":
                            if j > 1:
                                seq[trans_id] = len(sequence)
                            trans_id = m_.group(1)
                            gene_id = m_.group(2)
                            # if trans_id not in seq.keys():
                            #     seq[trans_id] = {}
                            sequence = ''
                        else:
                            if j>1:
                                seq[gene_id]=sequence
                            gene_id = m_.group(2)
                            # if gene_id not in seq.keys():
                            #     seq[gene_id]= {}
                            sequence = ''
                else:
                     sequence += line
                if _type == 'transcript':
                    seq[trans_id] = len(sequence)
                else:
                    seq[gene_id] = sequence
        if not seq:
            print '提取{}序列信息为空'.format(_type)
        print "共统计出{}行信息".format(str(j))
        end = time.time()
        duration = end - start
        m, s = divmod(duration, 60)
        h, m = divmod(m, 60)
        print('{}序列提取运行的时间为{}h:{}m:{}s'.format(_type,h, m, s))
        return seq

if __name__ == "__main__":
    # biomart_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/biomart/finish/Ensembl_Genes_89/oniloticus_gene_ensembl_gene.txt"
    # entrez_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/ncbi_gene2ensembl/gene2ensembl"
    # entrez_db_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/ncbi_gene2ensembl/gene2ensembl.db"
    # cds_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/Vertebrates/Fish/Oreochromis_niloticus/Oreochromis_niloticus.Orenil1.0.cds.all.fa"
    # pep_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/Vertebrates/Fish/Oreochromis_niloticus/Oreochromis_niloticus.Orenil1.0.pep.all.fa"
    # pep_info = get_pep_seq(pep_path)
    # # # insert_sequence_id(cds_info=cds_info, species_name="oreochromis_niloticus", sequence_type="pep", task_id=None, project_sn=None)
    # cds_info = get_cds_seq(cds_path)
    # insert_sequence_id(cds_info=cds_info, species_name="oreochromis_niloticus", sequence_type="cds", task_id=None,project_sn=None)
    # transcript_path = "/mnt/ilustre/users/sanger-dev/workspace/20170616/Single_transcript_abstract_1/TranscriptAbstract/output/exons.fa"
    # gene_path = "/mnt/ilustre/users/sanger-dev/workspace/20170616/Single_transcript_abstract_1/TranscriptAbstract/output/the_longest_exons.fa"
    # class_code = "/mnt/ilustre/users/sanger-dev/workspace/20170523/Single_merge_rsem_fpkm_11/MergeRsem/new_class_code"
    # transcript_path = "/mnt/ilustre/users/sanger-dev/workspace/20170619/Single_transcript_abstract_2/TranscriptAbstract/output/exons.fa"
    # gene_path = "/mnt/ilustre/users/sanger-dev/workspace/20170619/Single_transcript_abstract_2/TranscriptAbstract/output/the_longest_exons.fa"
    # class_code = "/mnt/ilustre/users/sanger-dev/workspace/20170523/Single_merge_rsem_fpkm_11/MergeRsem/new_class_code"
    # new_gene_path = "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/class_code/new_gene_location"
    # add_class_code_detail(class_code, class_code_id=None, biomart_path=biomart_path, cds_path=cds_path, pep_path=pep_path,
    #                           transcript_path=transcript_path, gene_path=gene_path, species='oreochromis_niloticus')


    # add_class_code(assembly_method="stringtie", name=None, major=True, biomart_path=biomart_path,new_gene_location_path=new_gene_path,
    #                class_code_path=class_code, cds_path = cds_path,pep_path=pep_path, species='oreochromis_niloticus',
    #                transcript_path=transcript_path,gene_path=gene_path)

    ###################################################################################################################################
    ########################---------------------------------导ore_test_for_api的测试数据
    ###################################################################################################################################
    # biomart_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/biomart/finish/Ensembl_Genes_89/oniloticus_gene_ensembl_gene.txt"
    # entrez_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/ncbi_gene2ensembl/gene2ensembl"
    # entrez_db_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/ncbi_gene2ensembl/gene2ensembl.db"
    # cds_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/Vertebrates/Fish/Oreochromis_niloticus/Oreochromis_niloticus.Orenil1.0.cds.all.fa"
    # pep_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/Vertebrates/Fish/Oreochromis_niloticus/Oreochromis_niloticus.Orenil1.0.pep.all.fa"
    #
    # transcript_path = "/mnt/ilustre/users/sanger-dev/workspace/20170613/Refrna_ore_test_for_api/TranscriptAbstract/output/exons.fa"
    # gene_path = "/mnt/ilustre/users/sanger-dev/workspace/20170613/Refrna_ore_test_for_api/TranscriptAbstract/output/the_longest_exons.fa"
    # class_code = "/mnt/ilustre/users/sanger-dev/workspace/20170613/Refrna_ore_test_for_api/Express/MergeRsem/new_class_code"
    # new_gene_path = '/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/class_code/ore_test_for_api.new_gene_location'

    biomart_path = '/mnt/ilustre/users/sanger-dev/app/database/refGenome/Vertebrates/Mammalia/Mus_musculus/ref/mmusculus_gene_ensembl_gene.txt'
    cds_path = '/mnt/ilustre/users/sanger-dev/app/database/refGenome/Vertebrates/Mammalia/Mus_musculus/Mus_musculus.GRCm38.cds.all.fa'
    pep_path = '/mnt/ilustre/users/sanger-dev/app/database/refGenome/Vertebrates/Mammalia/Mus_musculus/ref/Mus_musculus.GRCm38.pep.all.fa'
    transcript_path = '/mnt/ilustre/users/sanger-dev/workspace/20170629/Single_rsem_stringtie_mouse_total_1/Express/TranscriptAbstract/output/exons.fa'
    gene_path = '/mnt/ilustre/users/sanger-dev/workspace/20170629/Single_rsem_stringtie_mouse_total_1/Express/TranscriptAbstract/output/the_longest_exons.fa'
    class_code = '/mnt/ilustre/users/sanger-dev/workspace/20170630/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/MergeRsem/class_code'
    new_gene_path = '/mnt/ilustre/users/sanger-dev/workspace/20170630/Single_rsem_stringtie_mouse_fpkm_2_diff_stat_2/Expresstest3/MergeRsem/new_gene_location/new_gene_location'
    # a= RefrnaGeneDetail()
    # a.add_class_code(assembly_method="stringtie", name=None, major_gene_detail=True, major_express_diff=True, biomart_path=biomart_path,new_gene_location_path=new_gene_path,
    #                class_code_path=class_code, cds_path = cds_path,pep_path=pep_path, species='mus_musculus',
    #                transcript_path=transcript_path,gene_path=gene_path)
    # # a.add_express_diff_class_code_detail(class_code,'59573737a4e1af4b9a3ace22','mus_musculus')
    # gene_sequence_tmp = a.get_transcript_seq(gene_path,_type='gene',assembly_method='stringtie')
    # print "共获得{}个基因".format(str(len(data.keys())))
    # if 'ENSMUSG00000022517' in data.keys():
    #     print data['ENSMUSG00000022517']
    # else:
    #     print "{}没有获得对应的gene sequence".format('ENSMUSG00000022517')
    # print '导入差异分析的详情表成功！'
    # trans_seq = get_transcript_seq(transcript_path, _type="transcript")
    # print trans_seq
    # gene_seq = get_transcript_seq(gene_path, _type="gene")
    # print gene_seq
    ##########################################################################################################################################################################
    ##############################-----------对统计出的错误的基因序列进行更新
    ##########################################################################################################################################################################
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    # collection = db["sg_express_class_code_detail"]
    # data = collection.find({"class_code_id" : ObjectId("59573737a4e1af4b9a3ace22"),"type":"gene_detail"})
    # count = 0
    # for d in data:
    #     _id = d["_id"]
    #     gene_id = d["gene_id"]
    #     gene_sequence = d['gene_sequence']
    #     if gene_id in gene_sequence_tmp.keys():
    #         count +=1
    #         collection.update({"_id":_id},{"$unset":{"gene_sequence":gene_sequence}})
    #         collection.update({"_id": _id}, {"$set": {"gene_sequence": gene_sequence_tmp[gene_id]}})
    # print '共更新基因序列{}个'.format(str(count))
    db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    collection = db["sg_express_class_code_detail"]
    data = collection.find({"class_code_id" : ObjectId("5957539ba4e1af530d2373c0"),"type":"gene_detail"})
    count = 0
    import copy
    for d in data:
        _id = d["_id"]
        gene_id = d["gene_id"]
        if 'trans_info' in d.keys():
            trans_info = d['trans_info']
            new_transc_info = copy.deepcopy(trans_info)
            try:
                if isinstance(trans_info[trans_info.keys()[0]]['length'],int):
                    print trans_info[trans_info.keys()[0]]['length']
                    continue
                else:
                    for keys,values in trans_info.items():
                        new_transc_info[keys]['length'] = int(len(values['length']))
                    collection.update({"_id":_id},{"$unset":{"trans_info":trans_info}})
                    collection.update({"_id": _id}, {"$set": {"trans_info": new_transc_info}})
            except Exception:
                print d

    # biomart(biomart_path)
    # entrez_id_info = entrez_id(entrez_path,entrez_db_path)
    # insert_entrez_id(entrez_path, tasks_id='tsg_2000', project_sn='20000')
