# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/24 14:42


import argparse
import subprocess
import re
import time
import os,pickle

'''
此脚本主要用作cufflinks的merged.gtf 转换ID使用。
原则如下：
1. class code 为=时，认为此转录本是ref.gtf中的已知转录本，则gene/transcript name/ID 都要替换为ref.gtf中的ID和name。其中gene id 来源于 记录中gene name 在ref.gtf 对应的gene ID，
txpt 的ID为 记录中的oId
2. class code 为u时,不改变记录
3. 其他情况只替换其gene name
使用说明：
因merged.gtf较大，因此将其劈成指定行数的小文件集合，同步进行修饰，增加速度
参数：-s 处理脚本 必须在一起合用 是trans_script.py 的路径
  -tmp  name 暂存文件的文件夹名字  绝对路径为merged。gtf 的文件夹路径下的name文件夹
   -rm 是否结束程序后删除 暂存文件夹
   -lines merged gtf 被split后，每个小文件里有多少行内容
   
'''

# file_path = "F:\\code_lib\\merged.gtf"
# new_file_path = "F:\\code_lib\\new_merged.gtf"
# ref_gtf = "F:\\code_lib\\ref.gtf"

# =================函数区==================
def get_dic_from_list(**kwargs):
	d = {}
	lst = kwargs['lst']
	sep = kwargs['primary_sep']
	for record in lst:
		record_items = record.split(sep)
		gene_id = ''
		gname = ''
		for ele in record_items:
			ele_m = re.match(r'\s*(\S+)\s+"(\S+)"', ele)
			if ele_m:
				attr = ele_m.group(1)
				val = ele_m.group(2)
				if attr == 'gene_id':
					gene_id = val
					continue
				if attr == 'gene_name':
					gname = val
					continue
			else:
				print '{} 的{} 不正常'.format(record, ele)
		if gname and gene_id:
			d[gname] = gene_id
	return d


def split_merged_gtf(f, line_num, tmp_dir_name,name_prefix):
	tmp_dir = os.path.join(os.path.dirname(f), tmp_dir_name)
	abs_prefix = os.path.join(tmp_dir, name_prefix)
	if not os.path.isdir(tmp_dir):
		os.mkdir(tmp_dir)
	subprocess.call('split -l {} {} {}'.format(line_num, f, abs_prefix),shell=True)
	chunk_files = [os.path.join(tmp_dir, chunk) for chunk in os.listdir(tmp_dir)]
	return chunk_files, tmp_dir


# =============================主函数=====================
if __name__ == '__main__':
	time1 = time.time()
	parser = argparse.ArgumentParser(description="file")
	parser.add_argument("-merge", "--merged_gtf", help="input merged gtf ", required=True)
	parser.add_argument("-out_merge", "--new_merged_gtf", help="out new merged gtf", required=True)
	parser.add_argument("-ref_gtf", "--ref_gtf", help="ref_gtf", required=True)
	parser.add_argument("-method", "--method", help="cufflinks/stringtie", required=True)
	parser.add_argument("-s", "--trans_script", help="trans content script path", required=True)
	parser.add_argument("-tmp", "--tmp_dir", help="temp dir name", required=True)
	parser.add_argument("-rm", "--rm_tmp", help="rm temp dir or not,default is true(1),false is 0  ", default=1)
	parser.add_argument("-lines", "--split_line_number", help="split big file to  how many lines per file default is 40000 ", default=40000)
	args = vars(parser.parse_args())
	
	file_path = args["file_path"]
	new_file_path = args["new_file_path"]
	ref_gtf = args["ref_gtf"]
	method = args["method"]
	trans_script = args['trans_script']
	tmp_dir = args['tmp_dir']
	rm = args['rm_tmp']
	line_number = int(args['split_line_number'])
	
	ref_tmp_cmd = """ awk -F '\t|;' 'NF>=11{print $10";"$11}' %s  | uniq """ % (ref_gtf)
	print '开始awk ref.gtf'
	ref_tmp_content = subprocess.check_output(ref_tmp_cmd, shell=True).split('\n')
	print 'ref gtf 暂时文件读取完毕'
	gname_gid_dic = get_dic_from_list(lst=ref_tmp_content, primary_sep=';')
	print 'ref gtf 信息装载完毕'
	
	chunk_files_lst, tmp_folder = split_merged_gtf(file_path, line_number,tmp_dir, 'chunk_file_')
	ref_dic_pickle = os.path.join(tmp_folder,'ref_gname_gid_dic.pk')
	with open(ref_dic_pickle, 'wb') as handle:
		pickle.dump(gname_gid_dic, handle, protocol=pickle.HIGHEST_PROTOCOL)
	print  'ref 的gene name gene id 字典对象已保存，保存地址为{}'.format(ref_dic_pickle)
	pro_lst = []
	out_file_lst = []
	for chunk in chunk_files_lst:
		number = chunk_files_lst.index(chunk)+1
		out_file = os.path.join(tmp_folder, 'modified_' + os.path.basename(chunk))
		mod_gtf_cmd = 'python {}  -input_gtf  {} -out {} -ref_dic {} '.format(trans_script, chunk, out_file,ref_dic_pickle)
		child_pro = subprocess.Popen(mod_gtf_cmd, shell=True)
		print '开始修饰第{}个文件'.format(number)
		# print '命令为：{}'.format(number,mod_gtf_cmd)
		pro_lst.append(child_pro)
		out_file_lst.append(out_file)
	for pro in pro_lst:
		pro.communicate()
	out_file_str = "  ".join(out_file_lst)
	subprocess.call('cat {} > {}'.format(out_file_str,new_file_path),shell=True)
	print '已将最后结果汇聚为{}'.format(new_file_path)
	if rm:
		subprocess.call('rm -rf {}'.format(tmp_folder),shell=True)
	time2 = time.time()
	duration = time2 - time1
	m, s = divmod(duration, 60)
	h, m = divmod(m, 60)
	print '整个程序运行的时间为{}h:{}m:{}s'.format(h,m,s)
