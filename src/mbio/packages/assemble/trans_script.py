# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/25 18:34

import re, os, Bio, argparse, sys, fileinput, urllib2,pickle


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-input_gtf", "--gtf", help="file_path", required=True)
	parser.add_argument("-out", "--out", help="out_file_path", required=True)
	# parser.add_argument("-n", "--number", help="number th file", required=True)
	parser.add_argument("-ref_dic", "--ref_dic_pk", help="ref_gtf_pickle", required=True)
	args = vars(parser.parse_args())
	gtf = args['gtf']
	out = args['out']
	ref_dic_pk = args['ref_dic_pk']
	# number = args['number']
	
	with open(ref_dic_pk, 'rb') as handle:
		gname_gid_dic = pickle.load(handle)
	# print 'ref gene name 和gene id 对应dic 已经装载完毕，开始修饰第{}个chunk文件：{},结果文件为{}'.format(number,gtf,out)
	
	fw = open(out, 'w')
	for line in open(gtf):
		old_gid = ""
		old_tid = ""
		cls = ""
		gid = ""
		tid = ""
		oId = ""
		old_gid_m = re.search(r'gene_id\s+\"(\S+)\";', line)
		old_tid_m = re.search(r'transcript_id\s+\"(\S+)\";', line)
		cls_m = re.search(r'class_code\s+\"(\S+)\";', line)
		gname_m = re.search(r'gene_name\s+\"(\S+)\";', line)
		oId_m = re.search(r'oId\s+\"(\S+)\";', line)
		if old_gid_m:
			old_gid = old_gid_m.group(1)
		if old_tid_m:
			old_tid = old_tid_m.group(1)
		if cls_m:
			cls = cls_m.group(1)
		if gname_m:
			gname = gname_m.group(1)
		if oId_m:
			oId = oId_m.group(1)
		
		if re.search(r'class_code\s+\"=\";', line) and oId:
			if oId:
				tid = oId
			else:
				tid = old_tid
			if gname and (gname in gname_gid_dic.keys()):
				gid = gname_gid_dic[gname]
			else:
				gid = old_gid
				print '{}在ref.gtf中没有记录'
			newline = re.sub(r'(\s*)(gene_id)(\s+)\"(\S+)\";', '\g<1>\g<2>\g<3>\"{},{}\"'.format(gid, gname), line)
			newline = re.sub(r'(\s*)(transcript_id)(\s+)\"(\S+)\";', '\g<1>\g<2>\g<3>\"{},{}\"'.format(tid, gname),
			                 newline)
			fw.write(newline)
			continue
		if not re.search(r'class_code\s+\"[=u]\";', line):
			tid = old_tid
			if gname and (gname in gname_gid_dic.keys()):
				gid = gname_gid_dic[gname]
			else:
				gid = old_gid
			newline = re.sub(r'(\s*)(gene_id)(\s+)\"(\S+)\";', '\g<1>\g<2>\g<3>\"{},{}\"'.format(gid, gname), line)
			newline = re.sub(r'(\s*)(transcript_id)(\s+)\"(\S+)\";', '\g<1>\g<2>\g<3>\"{},{}\"'.format(tid, gname),
			                 newline)
			fw.write(newline)
			continue
		if re.search(r'class_code\s+\"u\";', line):
			fw.write(line)
			continue