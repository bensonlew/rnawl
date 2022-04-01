#!/bin/env python
# coding=utf-8

import os
import sys
# import xlrd
from collections import OrderedDict
import subprocess
import argparse
import json
import time	
import pymysql

def add_sourse_db(db):
	cwd = os.getcwd()
	hash_conn = 0
	n = 0
	while True:
		n = n + 1
		if n == 5:
			break
		
		try:
			conn = pymysql.connect(host="192.168.12.61", user="changliang.qiu", passwd="changliang.qiu123",port=4000,db="hdna",charset="utf8")
			break
		except:
			print("connect 192.168.12.61 fail, try again")
			time.sleep(3)
	if hash_conn == 1:
		with open(db,'r') as fin:
			for line in fin:
				sample_sourse_dict = OrderedDict()
				line = line.strip()
				sample_name = line.split()[0]
				clade = line.split()[1]
				snps = line.split()[2]
				analysis_path = os.path.join(cwd,sample_name)
				sample_sourse_dict["haplotype"] = clade + "-" + snps

				conn = pymysql.connect(host="192.168.12.46", user="changliang.qiu", passwd="changliang.qiu123",port=3306,db="hdna",charset="utf8")

				if conn:
					print ("database has already connected")
				else:
					print("database has already not connected")
			# with conn:
			# 	#conn.select_db('hdna')
				cur = conn.cursor()
				cur.execute("select analysis_path from Y_source_gene where analysis_path='{}'".format(analysis_path))
				result = cur.fetchone()
				if result:
					for key, value in sample_sourse_dict.items():
						cur.execute("update Y_source_gene set {}='{}' where analysis_path='{}'".format(key,value,analysis_path))
						conn.commit()
				else:
					cur.execute("insert into Y_source_gene set analysis_path='{}'".format(analysis_path))
					conn.commit()
					for key, value in sample_sourse_dict.items():
						cur.execute("update Y_source_gene set {}='{}' where analysis_path='{}'".format(key,value,analysis_path))
						conn.commit()
				cur.close()
				conn.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser()  # 建立解析器
	parser.add_argument("-i", help="单倍群结果表")
	args = parser.parse_args()  # 从外部传递参数到此
	add_sourse_db(args.i)


