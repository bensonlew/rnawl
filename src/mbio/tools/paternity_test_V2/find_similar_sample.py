#!/bin/env python
#coding=utf-8

"""
	这个程序是用来做个体识别的，主要是检查两个样品之间是不是相似。
	主要是利用存放起来的FATHER的tab文件。以后可能还会用到mother的tab文件
	把所有的结果计算一遍之后，按照相似性去前5个进行导表。
	但是本地上会存放所有结果的文档。
"""
from __future__ import print_function
import codecs
import os
import pymongo
import datetime
import re
from biocluster.config import Config
from operator import itemgetter, attrgetter
import time

def make_dict(input):
	dic = {}
	with codecs.open(input,"r","utf-8") as FA:
		for line in FA:
			if len(line) < 8:               #过滤掉列数小于8的位点
				continue
			else:
				line = line.strip()
				tmp = line.split("\t")
				sample_name, chrom, pos, ref, alt, dp, ref_dp, alt_dp = tmp
				if alt == ".":
					alt = ref
				if int(ref_dp) + int(alt_dp) < 30:    # 过滤掉深度小于30的位点
					continue
				alt1 = alt[0]       #取alt的第一个碱基作为突变碱基
				YeChun = "%s/%s" % (ref,ref)
				TuChun = "%s/%s" % (alt1,alt1)
				ZaHe  =  "%s/%s" % (ref,alt1)
				rf = int(ref_dp)/(int(ref_dp) + int(alt_dp))
				if len(ref) > 1:
					continue           #去除ref碱基数超过1的位点
				if rf >= 0.9:
					tmp.append(YeChun)
				elif rf <= 0.1:
					tmp.append(TuChun)
				else:
					tmp.append(ZaHe)
				sample_name, chrom, pos, ref, alt, dp, ref_dp, alt_dp, genetype = tmp
				content = sample_name+'\t'+chrom+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+dp+'\t'+ref_dp+'\t'+alt_dp+'\t'+genetype+'\n'
				dic[chrom+"_"+pos] = genetype              #以染色体及位置作为键,基因型作为值,生成字典
	return dic

def get_new_sample(sample_id):
	father_list = []
	mother_list = []
	db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
	collection = db['sg_sample_similarity']
	library_path = "/mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/tab_data_bak"
	dir_result = os.walk(library_path)
	for root, dirs, files in dir_result:
		for each in files:
			if re.search(r"-F",each) != None:
				father_list.append(each)
			elif re.search(r"-M",each) != None:
				mother_list.append(each)
	####获取新增样品
	sample_id_mongo = collection.find_one({'sample_id':sample_id})
	if sample_id_mongo != None:
		new_sample = list(set(father_list).difference(set(sample_id_mongo["dec"]))) ##只要参考文件夹内，与monggo数据库中不同的样品，都算是新增样品
		old_similarity = sample_id_mongo["result"]
	else:
		new_sample = father_list
		old_similarity = "empty"
	abs_sample = []
	for each in new_sample:
		if re.search(r"{}".format(sample_id),each) == None:
			abs_name = "{}/{}".format(library_path,each)
			abs_sample.append(abs_name)
	return([abs_sample,old_similarity])

def Identification(current_sample,sample_list,sample_name,old_similarity):
	begin_time = time.time()
	overall_num = len(sample_list)
	run_num = 1
	f1 = codecs.open("{}.tmpfile.xls".format(sample_name),'w',"utf-8")
	current_dict = make_dict(current_sample)	         #生成当前样本的字典
	###写入中间文件###
	if old_similarity != "empty":
		for each in old_similarity:
			content = "\t".join(each)
			print(content,file=f1)
	for sn in sample_list:
		sn =sn.strip()
		other_dict = make_dict(sn)                              #生成所需要对比的字典
		same_gene = 0
		all_gene = 0
		for cur_k,cur_v  in current_dict.iteritems():
			if other_dict.has_key(cur_k):
				all_gene  += 1                              #当对比的两个样本的位置一致时,累加,作为分母
				if cur_v == other_dict[cur_k]:
					same_gene += 1                          #当位置一致的时候,基因型一致时,累加,作为分子
				else:
					continue
			else:
				continue
		if all_gene > 0 :
			percent = '{:.4f}'.format(100*float(same_gene)/float(all_gene))+"%"
		else:
			percent = 0
		content = [str(sn),str(all_gene),str(same_gene),str(percent)]
		content = "\t".join(content)
		print(content,file=f1)
		#print("{}/{}".format(run_num,overall_num))
		run_num += 1                                       #将计算出的每一行结果写入临时文件中,后续对匹配度进行排序工作
	f1.close()
	cmd = "less {}.tmpfile.xls|sort -nrk2 > {}.outfile.xls".format(sample_name,sample_name)
	os.system(cmd)
	finish_time = time.time()
	overall_time  = finish_time - begin_time
	print(overall_time)

def dict_test(current_sample,sample_list):
	begin_time = time.time()
	overall_num = len(sample_list)
	run_num = 1
	current_dict = make_dict(current_sample)	         #生成当前样本的字典     
	for sn in sample_list:
		other_dict = make_dict(sn)
		print("{}/{}".format(run_num,overall_num))
		run_num += 1   
	finish_time = time.time()
	overall_time  = finish_time - begin_time
	print(overall_time)

def dump_similarity(sample_name,old_similarity):
	sample_sort_file = "{}.outfile.xls".format(sample_name)
	count = 0;
	similarity_data = []
	similarity_sample = []
	with codecs.open(sample_sort_file,"r","utf-8") as IN:
		for line in IN:
			line = line.strip()
			linetemp = line.split("\t")
			linetemp[0] = linetemp[0].split("/")[-1].replace('.tab','')
			similarity_sample.append(linetemp[0])
			count += 1
			if count <= 10:
				similarity_data.append(linetemp)
	insert_table = {
		"sample_id": sample_name,
		"upload_time":datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
		"result":similarity_data,
		"dec":similarity_sample
	}
	db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
	collection = db['sg_sample_similarity']
	collection.update({"sample_id":sample_name},{'$set':insert_table},upsert = True)

###test
group_return = get_new_sample("WQ17083258-F1.tab")	
father_list = group_return[0]
old_similarity = group_return[1]
Identification("WQ17083258-F1.tab",father_list,"WQ17083258-F1",old_similarity)
#dict_test("WQ17083258-F1.tab",father_list)
dump_similarity("WQ17083258-F1",old_similarity)
