# coding=utf-8
import xml.etree.ElementTree as ET
from biocluster.config import Config
import re
import collections
import json
from itertools import islice
import subprocess
import gridfs
import os
import sys
import pandas as pd
import json
import math
import numpy as np
from scipy import stats
import xml.etree.ElementTree as ET
import lxml.html
from collections import OrderedDict
from scipy.stats import pearsonr
import random
from bson.son import SON
from mbio.packages.rna.dendrogram2newick import convert2
import csv

class Chart(object):
	def __init__(self):
		"""
		设置数据库，连接到mongod数据库，kegg_ko,kegg_gene,kegg_pathway_png三个collections2
		本itraq_tmt静态图开发所使用的测试项目为：http://report.nsg.com/itraq/proteinsetgo_two/task_id/s5jn_s0ieqm5rf9dfelc1t88ajh.html
		开发的静态图都尽量与该项目的交互页面的图片一致。
		"""

		"""
		class: sample[A1, B1], group[A, B], classify: [samples, groups], level[G, T], kind[ref, new, all], category[mRNA, lncRNA]
		"""
		self.sample_list = list()
		self.command_list = list()
		# chart_dir = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/work/sg_chart/test2"
		chart_dir = Config().SOFTWARE_DIR + "/bioinfo/sg_chart"
		self.phantomjs_dir = Config().SOFTWARE_DIR + "/program/phantomjs/phantomjs-2.1.1-linux-x86_64/bin/"
		self.wkhtml_dir = Config().SOFTWARE_DIR + "/install_packages/wkhtmltox-0.12.6-1/usr/local/bin/"
		self.node_path = Config().SOFTWARE_DIR + "/bioinfo/sg_chart/node-v14.16.0-linux-x64/bin/node"
		self.puppeteer_path = Config().SOFTWARE_DIR + "/bioinfo/sg_chart/model/sg_chart_puppeteer.js"
		self.puppeteer_dir = Config().SOFTWARE_DIR + "/install_packages/wkhtmltox-0.12.6-1/usr/local/bin/"
		self.mode_dir = chart_dir + "/dia_v3"
		self.mode_mode = chart_dir + "/model"
		self.mode_mode2 = chart_dir + "/dia_v3/highchart_model"
		self.js_list = list()
		self.work_dir = os.environ.get("PWD", "") + "/"
		self.output_dir = os.environ.get("PWD", "") + "/output/"
		self.rcolor_dict = self.get_rcolor_dict()

	def chart_protein_mw(self, Protein_mw_distribution_xls):
		mw_distribution_df = pd.read_table(Protein_mw_distribution_xls)
		category_list = mw_distribution_df.iloc[:,0].values.tolist()
		num_list = mw_distribution_df.iloc[:,1].values.tolist()
		self.chart_protein_mw_tojs("mw", category_list, num_list, "protein_mw_distribution.json")

	def chart_protein_mw_tojs(self, name, category_list, num_list, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["categories"] = category_list
			a["data"] = [num_list]
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":800,"height":500}])

	def chart_protein_seq_cover(self, Protein_seq_cover_distribution_xls):
		mw_distribution_df = pd.read_table(Protein_seq_cover_distribution_xls)
		coverage_list = mw_distribution_df.iloc[:,0].values.tolist()
		num_list = mw_distribution_df.iloc[:,1].values.tolist()
		dic_list = [{'y':num,'name':cov} for num,cov in zip(num_list, coverage_list)]
		self.chart_protein_seq_cover_tojs("seq_cov", dic_list, "protein_seq_cover_distribution.json")

	def chart_protein_seq_cover_tojs(self, name, dic_list, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = dic_list
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showPie","width":600,"height":400}])

	def chart_protein_infomation(self, Protein_infomation_xls):
		infomation_df = pd.read_csv(Protein_infomation_xls)
		num_list = infomation_df.iloc[0,:].values.tolist()
		self.chart_protein_infomation_tojs("info", num_list, "protein_infomation.json")

	def chart_protein_infomation_tojs(self, name, num_list, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [num_list]
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":600,"height":400}])

	def chart_peptite_number_distribution(self, Peptite_number_distribution_xls):
		peptite_number_distribution_df = pd.read_table(Peptite_number_distribution_xls)
		category_list = peptite_number_distribution_df.iloc[:,0].values.tolist()
		num_list = peptite_number_distribution_df.iloc[:,1].values.tolist()
		self.chart_peptite_number_distribution_tojs("pep_num", category_list, num_list, "peptite_number_distribution.json")

	def chart_peptite_number_distribution_tojs(self, name, category_list, num_list, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["categories"] = category_list
			a["data"] = [num_list]
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":1000,"height":500}])

	def chart_peptite_length_distribution(self, Peptite_length_distribution_xls):
		peptite_length_distribution_df = pd.read_table(Peptite_length_distribution_xls)
		category_list = peptite_length_distribution_df.iloc[:,0].values.tolist()
		num_list = peptite_length_distribution_df.iloc[:,1].values.tolist()
		self.chart_peptite_length_distribution_tojs("pep_len", category_list, num_list, "peptite_length_distribution.json")

	def chart_peptite_length_distribution_tojs(self, name, category_list, num_list, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["categories"] = category_list
			a["data"] = [num_list]
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":1000,"height":600}])

	# def chart_peptide_dmass(self, Peptide_dmass_xls):
	# 	peptide_dmass_df = pd.read_table(Peptide_dmass_xls)
	# 	da_list = peptide_dmass_df.iloc[:,0].values.tolist()
	# 	ppm_list = peptide_dmass_df.iloc[:,1].values.tolist()
	# 	list_list = [[da,ppm] for da,ppm in zip(da_list, ppm_list)]
	# 	self.chart_peptide_dmass_tojs("pep_error", list_list, "peptide_dmass.json")

	# def chart_peptide_dmass_tojs(self, name, list_list, json_mode):
	# 	json_mode = self.mode_dir + "/" + json_mode
	# 	with open(json_mode, 'r') as f, open(self.output_dir + name + ".scatter.js", 'w') as fo:
	# 		a = json.loads(f.read())
	# 		a["data"] = [list_list]
	# 		fo.write("var options = ")
	# 		fo.write(json.dumps(a, indent=4))
	# 		self.js_list.append([self.output_dir + name + ".scatter.js", {"model":"highchart","highchart_type":"plot_big_scatter","width":1000,"height":500}])

	def chart_all_annotation_stat(self, all_annotation_statistics_xls):
		all_annotation_stat_df = pd.read_table(all_annotation_statistics_xls,nrows=5,index_col=False)
		categories_list = all_annotation_stat_df.iloc[:,0].values.tolist()
		num_list = all_annotation_stat_df.iloc[:,1].values.tolist()
		self.chart_all_annotation_stat_tojs("anno_stat", categories_list, num_list, "all_annotation_statistics.json")

	def chart_all_annotation_stat_tojs(self, name, categories_list, num_list, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [num_list]
			a["categories"] = categories_list
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":600,"height":400}])

	def chart_go12level_statistics(self, golevel_statistics_xls):
		golevel_statistics_df = pd.read_table(golevel_statistics_xls,dtype={'Percent':object})
		def get_ten(x):
		  df = x.drop_duplicates(subset=None,keep='first',inplace=False).sort_values(by = [x.columns[2],x.columns[1]], ascending=[False, True])
		  return df.iloc[:10,:]
		list_list_ = golevel_statistics_df.iloc[:,[0,1,3,4]].groupby(golevel_statistics_df.columns[0],as_index=False).apply(get_ten)
		list_list = list_list_.values.tolist()
		self.chart_golevel_statistics_tojs("go_lev2_bar", list_list, "golevel_statistics.json")
		
		def t_l(x):
		  return [{"name":str(f+1)+" : "+i[1],"value":i[2]} for f,i in enumerate(x.values.tolist())]
		list_list_for_pie = list_list_.groupby(list_list_.columns[0],as_index=False).apply(t_l).values.tolist()
		self.chart_golevel_statistics_pie_tojs("go_lev2_pie", list_list_for_pie, "golevel_statistics_pie.json")

	def chart_go123level_statistics(self, golevel_statistics_xls):
		golevel_statistics_df = pd.read_table(golevel_statistics_xls,dtype={'Percent':object})
		def get_ten(x):
		  df = x.drop_duplicates(subset=None,keep='first',inplace=False).sort_values(by = [x.columns[2],x.columns[1]], ascending=[False, True])
		  return df.iloc[:10,:]
		list_list_ = golevel_statistics_df.iloc[:,[0,3,5,6]].groupby(golevel_statistics_df.columns[0],as_index=False).apply(get_ten)
		list_list = list_list_.values.tolist()
		self.chart_golevel_statistics_tojs("go_lev3_bar", list_list, "golevel_statistics.json")
		
		def t_l(x):
		  return [{"name":str(f+1)+" : "+i[1],"value":i[2]} for f,i in enumerate(x.values.tolist())]
		list_list_for_pie = list_list_.groupby(list_list_.columns[0],as_index=False).apply(t_l).values.tolist()
		self.chart_golevel_statistics_pie_tojs("go_lev3_pie", list_list_for_pie, "golevel_statistics_pie.json")

	def chart_go1234level_statistics(self, golevel_statistics_xls):
		golevel_statistics_df = pd.read_table(golevel_statistics_xls,dtype={'Percent':object})
		def get_ten(x):
		  df = x.drop_duplicates(subset=None,keep='first',inplace=False).sort_values(by = [x.columns[2],x.columns[1]], ascending=[False, True])
		  return df.iloc[:10,:]
		list_list_ = golevel_statistics_df.iloc[:,[0,5,7,8]].groupby(golevel_statistics_df.columns[0],as_index=False).apply(get_ten)
		list_list = list_list_.values.tolist()
		self.chart_golevel_statistics_tojs("go_lev4_bar", list_list, "golevel_statistics.json")

		def t_l(x):
		  return [{"name":str(f+1)+" : "+i[1],"value":i[2]} for f,i in enumerate(x.values.tolist())]
		list_list_for_pie = list_list_.groupby(list_list_.columns[0],as_index=False).apply(t_l).values.tolist()
		self.chart_golevel_statistics_pie_tojs("go_lev4_pie", list_list_for_pie, "golevel_statistics_pie.json")

	def chart_golevel_statistics_tojs(self, name, list_list, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".go_bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = list_list
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4, sort_keys=True))
			self.js_list.append([self.output_dir + name + ".go_bar.js", {"model":"highchart","highchart_type":"go_bar","width":1100,"height":600}])
	
	def chart_golevel_statistics_pie_tojs(self, name, list_list, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".go_pie.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = list_list
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4, sort_keys=True))
			self.js_list.append([self.output_dir + name + ".go_pie.js", {"model":"highchart","highchart_type":"go_bar","width":1500,"height":550}])

	def chart_kegg_layer(self, kegg_layer_xls):
		kegg_layer_df = pd.read_table(kegg_layer_xls, header=None)
		kegg_type = kegg_layer_df.iloc[:,0].values.tolist()
		categories = kegg_layer_df.iloc[:,1].values.tolist()
		data = kegg_layer_df.iloc[:,2].values.tolist()
		color_dic = {"Metabolism":"#F44336", "Genetic Information Processing":"#388E3C", "Environmental Information Processing":"#FF00FF", "Cellular Processes":"#0288D1", "Organismal Systems":"#ff0", "Human Diseases":"#FF9800"}
		taxon_data = []
		data_colors = []
		for i in sorted(set(kegg_type), key=kegg_type.index):
			tmp_index = [idx for idx, e in enumerate(kegg_type) if e==i]
			taxon_data.append({
				"name":i,
				"color":color_dic[i],
				"start":tmp_index[0],
				"end":tmp_index[-1]
			})
			data_colors.extend([color_dic[i]]*len(tmp_index))

		self.chart_kegg_layer_tojs("path_class", data, categories, taxon_data, data_colors, "kegg_layer.json")

	def chart_kegg_layer_tojs(self, name, data, categories, taxon_data, data_colors, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".barline.js", 'w') as fo:
			a = json.loads(f.read())
			a["data_colors"] = data_colors
			a["taxon_data"] = taxon_data
			a["data"] = [data]
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".barline.js", {"model":"highchart","highchart_type":"showBarLine","width":1000,"height":600}])

	def chart_kegg_pathway(self,kegg_pathway_xls):
		kegg_pathway_df = pd.read_table(kegg_pathway_xls)
		categories = kegg_pathway_df.sort_values(by = kegg_pathway_df.columns[4],ascending=False).iloc[:20,3].values.tolist()
		num = kegg_pathway_df.sort_values(by = kegg_pathway_df.columns[4],ascending=False).iloc[:20,4].values.tolist()
		self.chart_kegg_pathway_tojs("key_path", categories, num, "kegg_pathway.json")

	def chart_kegg_pathway_tojs(self, name, categories, num, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [num]
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":1000,"height":600}])

	def chart_cog_summary(self,cog_summary_xls):
		cog_summary_df = pd.read_table(cog_summary_xls, skiprows=1)
		categories = [i[1] for i in cog_summary_df.iloc[:,1].values.tolist()]
		num = cog_summary_df.iloc[:,2].values.tolist()
		self.chart_cog_summary_tojs("cog_bar", categories, num, "cog_summary.json")

	def chart_cog_summary_tojs(self, name, categories, num, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [num]
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":800,"height":600}])

	def chart_pfam_domain(self,pfam_domain_xls):
		pfam_domain_df = pd.read_table(pfam_domain_xls)
		pfam_series =  pd.value_counts(pfam_domain_df.iloc[:,3]).head(20)
		categories = pfam_series.index.tolist()
		num = pfam_series.tolist()
		self.chart_pfam_domain_tojs("pfam_bar", categories, num, "pfam_domain.json")

	def chart_pfam_domain_tojs(self, name, categories, num, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [num]
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":800,"height":500}])

	def chart_subloc_stat(self,subloc_stat_xls):
		subloc_stat_df = pd.read_table(subloc_stat_xls)
		categories =  subloc_stat_df.iloc[:,0].values.tolist()
		num = subloc_stat_df.iloc[:,1].values.tolist()
		self.chart_subloc_stat_tojs("subloc_bar", categories, num, "subloc_stat.json")

	def chart_subloc_stat_tojs(self, name, categories, num, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [num]
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":800,"height":500}])

	def chart_diff(self,group_list,num_summary_file,allsummary_file,group_file_list):
		num_summary_df = pd.read_table(num_summary_file)
		allsummary_df = pd.read_table(allsummary_file)
		for i,this_group_name in enumerate(group_list):
			tmp_group_name = this_group_name.replace("_vs_", "|")
			categories = num_summary_df[tmp_group_name].tolist()
			this_group_file = group_file_list[i]
			diff_df = pd.read_table(this_group_file)
			diff_df = diff_df.merge(allsummary_df, how='left', on=diff_df.columns[0])
			data_nosig = diff_df[diff_df[tmp_group_name] == 'no'][[diff_df.columns[2], diff_df.columns[1], diff_df.columns[0], diff_df.columns[5], diff_df.columns[9]]].values.tolist()
			data_down = diff_df[diff_df[tmp_group_name+'_down'] == 'yes'][[diff_df.columns[2], diff_df.columns[1], diff_df.columns[0], diff_df.columns[5], diff_df.columns[9]]].values.tolist()
			data_up = diff_df[diff_df[tmp_group_name+'_up'] == 'yes'][[diff_df.columns[2], diff_df.columns[1], diff_df.columns[0], diff_df.columns[5], diff_df.columns[9]]].values.tolist()
			data = [data_nosig,data_down,data_up]
			self.chart_diff_tojs(this_group_name+'_scatter', categories, data, this_group_name, "diff.json")
			#
			diff_df['-log2fc'] = -diff_df.iloc[:,[9]]
			data_nosig_volcano = diff_df[diff_df[tmp_group_name] == 'no'][[diff_df.columns[5], diff_df.columns[-1], diff_df.columns[0], diff_df.columns[5], diff_df.columns[9]]].values.tolist()
			data_down_volcano = diff_df[diff_df[tmp_group_name+'_down'] == 'yes'][[diff_df.columns[5], diff_df.columns[-1], diff_df.columns[0], diff_df.columns[5], diff_df.columns[9]]].values.tolist()
			data_up_volcano = diff_df[diff_df[tmp_group_name+'_up'] == 'yes'][[diff_df.columns[5], diff_df.columns[-1], diff_df.columns[0], diff_df.columns[5], diff_df.columns[9]]].values.tolist()
			data_volcano = [data_nosig_volcano,data_down_volcano,data_up_volcano]
			self.chart_diff_volcano_tojs(this_group_name+'_volcano', categories, data_volcano, this_group_name, "diff_volcano.json")

	def chart_diff_tojs(self, name, categories, data, this_group_name, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".scatter.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["categories"] = categories
			a["params"]["title"] = this_group_name + ".Scatter"
			a["params"]["x_label"] = this_group_name.split('_vs_')[0]
			a["params"]["y_label"] = this_group_name.split('_vs_')[1]
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".scatter.js", {"model":"highchart","highchart_type":"plot_big_scatter","width":600,"height":480}])

	def chart_diff_volcano_tojs(self, name, categories, data, this_group_name, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".volcano.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["categories"] = categories
			a["params"]["title"] = this_group_name + ".volcano"
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".volcano.js", {"model":"highchart","highchart_type":"plot_big_scatter","width":600,"height":480}])

	def chart_proteinsetcluster(self, proteinsetcluster_group_list, expression_matrix_files_list, seq_tree_files_list, sample_tree_files_list):
		for i,this_group_name in enumerate(proteinsetcluster_group_list):
			this_expression_matrix_file = expression_matrix_files_list[i]
			this_seq_tree_file = seq_tree_files_list[i]
			this_sample_tree_file = sample_tree_files_list[i]
			this_seq_tree_file_str = open(this_seq_tree_file, 'r').read()
			this_tree_1 = this_seq_tree_file_str.split('\n')[0]
			this_tree_rows = this_seq_tree_file_str.split('\n')[1].strip().split(';')
			this_sample_tree_file_str = open(this_sample_tree_file, 'r').read()
			this_tree_2 = this_sample_tree_file_str.split('\n')[0]
			this_tree_columns = this_sample_tree_file_str.split('\n')[1].strip().split(';')

			this_expression_matrix_df = pd.read_table(this_expression_matrix_file, index_col=0)
			this_expression_matrix_df = this_expression_matrix_df[this_tree_columns]		# sort columns by sample tree
			heatmap_data = this_expression_matrix_df[this_expression_matrix_df.index.isin(this_tree_rows)].reindex(this_tree_rows).values.tolist()

			self.chart_proteinsetcluster_tojs(this_group_name+'____'+this_group_name, this_tree_1, this_tree_2, this_tree_rows, this_tree_columns, heatmap_data, "proteinsetcluster_heatmap.json")

			#______子聚类趋势图________
			this_group_dir = os.path.dirname(this_expression_matrix_file)
			tmp_file_prefix_list = ['seq.subcluster_'+str(i+1)+'_' for i in range(int(this_group_name.split('_')[-1]))]
			sub_files_path = []
			for file in os.listdir(this_group_dir):
				for prefix in tmp_file_prefix_list:
					if file.startswith(prefix):
						sub_files_path.append(os.path.join(this_group_dir,file))
			for sub_file in sub_files_path:
				this_sub_file_df = pd.read_table(sub_file, index_col=0)
				this_sub_file_df.loc['mean'] = this_sub_file_df.mean()
				data = this_sub_file_df.values.tolist()
				tmp_n = os.path.basename(sub_file)[4:-4].split('_')
				title = tmp_n[0]+'_'+tmp_n[1]+'('+tmp_n[2]+')'
				symbols = [""]*int(tmp_n[2])+['circle']
				colors = ["#F0F0F0"]*int(tmp_n[2])+['#000099']
				legend = this_sub_file_df.index.tolist()
				categories = this_sub_file_df.columns.tolist()
				self.chart_proteinsetcluster_sub_tojs(this_group_name+'____'+tmp_n[0][:3]+tmp_n[1], data, title, symbols, colors, legend, categories, "proteinsetcluster_SubLine.json")
				

	def chart_proteinsetcluster_tojs(self, name, this_tree_1, this_tree_2, this_tree_rows, this_tree_columns, heatmap_data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".heatmap.js", 'w') as fo:
			a = json.loads(f.read())
			a["tree_2"] = this_tree_2
			a["tree_1"] = this_tree_1
			a["rows"] = this_tree_rows
			a["columns"] = this_tree_columns
			a["heatmap_data"] = heatmap_data
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".heatmap.js", {"model":"highchart","highchart_type":"tree_heatmap","width":950,"height":5940,"use_puppeteer":'yes'}])

	def chart_proteinsetcluster_sub_tojs(self, name, data, title, symbols, colors, legend, categories, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".line.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["params"]['title'] = title
			a["params"]["symbols"] = symbols
			a["params"]["colors"] = colors
			a["legend"] = legend
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".line.js", {"model":"highchart","highchart_type":"showCurve","width":600,"height":400}])

	def chart_proteinsetcluster_goclass(self, diff_group_list, proteinsetgo_class_all_files_list, proteinsetgo_class_updown_files_list):
		for i,this_group_name in enumerate(diff_group_list):
			proteinsetgo_class_all_file = proteinsetgo_class_all_files_list[i]
			proteinsetgo_class_updown_file = proteinsetgo_class_updown_files_list[i]
			# proteinsetgo_class_all_df = pd.read_table(golevel_statistics_xls,dtype={'Percent':object})
			proteinsetgo_class_all_df = pd.read_table(proteinsetgo_class_all_file)
			proteinsetgo_class_updown_df = pd.read_table(proteinsetgo_class_updown_file)
			
			#____________单向直方图-all__________
			proteinsetgo_class_all_df_head20 = proteinsetgo_class_all_df.iloc[:,[0,1,3,4]].drop_duplicates(subset=None,keep='first',inplace=False).sort_values(by = [proteinsetgo_class_all_df.columns[3],proteinsetgo_class_all_df.columns[1]], ascending=[False,True]).iloc[:20,:]
			proteinsetgo_class_all_df_head20[proteinsetgo_class_all_df_head20.columns[3]] = proteinsetgo_class_all_df_head20[proteinsetgo_class_all_df_head20.columns[3]].map(lambda x:x.split('(')[0]).astype(float)
			list_list = proteinsetgo_class_all_df_head20.values.tolist()
			self.chart_proteinsetcluster_goclass_tojs(this_group_name+"____"+'_all_go', list_list, "proteinsetgo_class.json")
			self.chart_proteinsetcluster_goclass_tojs(this_group_name + "____" + this_group_name + '_all_go', list_list, "proteinsetgo_class.json", y_right=True)

			#____________单向直方图-up_down__________
			# proteinsetgo_class_updown_df['tmp_sum'] = proteinsetgo_class_updown_df.iloc[:,3] + proteinsetgo_class_updown_df.iloc[:,6]
			proteinsetgo_class_updown_df_head20 = proteinsetgo_class_updown_df.sort_values(by = [proteinsetgo_class_updown_df.columns[3], proteinsetgo_class_updown_df.columns[1]], ascending=[False,True]).iloc[:,[0,1,3,4,6,7]].drop_duplicates(subset=None,keep='first',inplace=False).iloc[:20,:]
			proteinsetgo_class_updown_df_head20[proteinsetgo_class_updown_df_head20.columns[3]] = proteinsetgo_class_updown_df_head20[proteinsetgo_class_updown_df_head20.columns[3]].map(lambda x:x.split('(')[0]).astype(float)
			proteinsetgo_class_updown_df_head20[proteinsetgo_class_updown_df_head20.columns[5]] = proteinsetgo_class_updown_df_head20[proteinsetgo_class_updown_df_head20.columns[5]].map(lambda x:x.split('(')[0]).astype(float)
			list_list = proteinsetgo_class_updown_df_head20.values.tolist()
			sub_cat = [this_group_name + '_up', this_group_name + '_down']
			self.chart_proteinsetcluster_goclass_tojs(this_group_name + "____" + this_group_name + '_up_down_go', list_list, "proteinsetgo_class_updown.json", sub_cat=sub_cat)

			proteinsetgo_class_up_df_head20 = proteinsetgo_class_updown_df.iloc[:, [0, 1, 3, 4]].drop_duplicates(subset=None, keep='first', inplace=False).sort_values(by=[proteinsetgo_class_updown_df.columns[3], proteinsetgo_class_updown_df.columns[1]], ascending=[False, True]).iloc[:20, :]
			proteinsetgo_class_up_df_head20[proteinsetgo_class_up_df_head20.columns[3]] = proteinsetgo_class_up_df_head20[proteinsetgo_class_up_df_head20.columns[3]].map(lambda x: x.split('(')[0]).astype(float)
			up_list = proteinsetgo_class_up_df_head20.values.tolist()
			self.chart_proteinsetcluster_goclass_tojs(this_group_name + "____" + this_group_name + '_up_go', up_list,"proteinsetgo_class_updown.json", y_right=True)
			proteinsetgo_class_down_df_head20 = proteinsetgo_class_updown_df.iloc[:, [0, 1, 6, 7]].drop_duplicates(subset=None, keep='first', inplace=False).sort_values(by=[proteinsetgo_class_updown_df.columns[6], proteinsetgo_class_updown_df.columns[1]], ascending=[False, True]).iloc[:20, :]
			proteinsetgo_class_down_df_head20[proteinsetgo_class_down_df_head20.columns[3]] = proteinsetgo_class_down_df_head20[proteinsetgo_class_down_df_head20.columns[3]].map(lambda x: x.split('(')[0]).astype(float)
			down_list = proteinsetgo_class_down_df_head20.values.tolist()
			self.chart_proteinsetcluster_goclass_tojs(this_group_name + "____" + this_group_name + '_down_go', down_list, "proteinsetgo_class_updown.json", y_right=True)

			#___________双向直方图-up_down___________
			proteinsetgo_class_updown_td_df_head20_ = proteinsetgo_class_updown_df.iloc[:,[0,1,3,6]].drop_duplicates(subset=None,keep='first',inplace=False).sort_values(by = [proteinsetgo_class_updown_df.columns[3], proteinsetgo_class_updown_df.columns[6], proteinsetgo_class_updown_df.columns[1]], ascending=[False, False, True]).iloc[:20,:]
			def get_sort(x):
			  df = x.sort_values(by=[x.columns[2], x.columns[3], x.columns[1]], ascending=[False, False, True])
			  return df
			proteinsetgo_class_updown_td_df_head20 = proteinsetgo_class_updown_td_df_head20_.groupby(proteinsetgo_class_updown_td_df_head20_.columns[0],as_index=False).apply(get_sort)
			color_dic = {'biological_process':'#FF00FF', 'cellular_component':'#388E3C', 'molecular_function':'#F44336'}
			abbr_dic = {'biological_process':'BP', 'cellular_component':'CC', 'molecular_function':'MF'}
			categories = [i[1]+':'+abbr_dic[i[0]] for i in proteinsetgo_class_updown_td_df_head20.iloc[:,[0,1]].values.tolist()]
			categories_color = [{'color':color_dic[i[0]], 'group':abbr_dic[i[0]], 'value':i[1]+':'+abbr_dic[i[0]]} for i in proteinsetgo_class_updown_td_df_head20.iloc[:,[0,1]].values.tolist()]
			colors_range_ = [i['color'] for i in categories_color]
			colors_range = []
			colors_range_1 = [i for i,val in enumerate(colors_range_) if val=='#FF00FF']
			# colors_range_11 = {'start':min(colors_range_1),'end':max(colors_range_1),'top_color':'#FF00FF','bottom_color':'#FF00FF'}
			if colors_range_1:
				colors_range_11 = {'start':min(colors_range_1),'end':max(colors_range_1),'top_color':'#FF00FF','bottom_color':'#FF00FF'}
				colors_range.append(colors_range_11)
			colors_range_2 = [i for i,val in enumerate(colors_range_) if val=='#388E3C']
			# colors_range_22 = {'start':min(colors_range_2),'end':max(colors_range_2),'top_color':'#388E3C','bottom_color':'#388E3C'}
			if colors_range_2:
				colors_range_22 = {'start':min(colors_range_2),'end':max(colors_range_2),'top_color':'#388E3C','bottom_color':'#388E3C'}
				colors_range.append(colors_range_22)
			colors_range_3 = [i for i,val in enumerate(colors_range_) if val=='#F44336']
			# colors_range_33 = {'start':min(colors_range_3),'end':max(colors_range_3),'top_color':'#F44336','bottom_color':'#F44336'}
			# colors_range = [colors_range_11,colors_range_22,colors_range_33]
			if colors_range_3:
				colors_range_33 = {'start':min(colors_range_3),'end':max(colors_range_3),'top_color':'#F44336','bottom_color':'#F44336'}
				colors_range.append(colors_range_33)

			proteinsetgo_class_updown_td_df_head20.iloc[:,[3]] = -proteinsetgo_class_updown_td_df_head20.iloc[:,[3]]
			data = proteinsetgo_class_updown_td_df_head20.iloc[:,[2,3]].T.values.tolist()
			self.chart_proteinsetgo_class_updown_td_tojs(this_group_name+"____"+this_group_name+'_up_down_go_neg', data, colors_range, categories, categories_color, "proteinsetgo_class_updown_td.json")

	def chart_proteinsetcluster_goclass_web(self, proteinsetgo_class_updown_file):
		proteinsetgo_class_updown_df = pd.read_table(proteinsetgo_class_updown_file, index_col=False)

		#____________单向直方图-up_down__________
		colnum = int(proteinsetgo_class_updown_df.shape[1])/3
		need_list_1 = [(i+1)*3-3 for i in range(colnum)][1:]   #[3,6]
		need_list_2 = list(set(list(range(colnum*3)))-set([(i+1)*3-1 for i in range(colnum)]))   #[0,1,3,4,6,7]
		need_list_3 = [(i+1)*2+1 for i in range(colnum-1)]  #[3,5]
		print need_list_1,need_list_2,need_list_3

		sub_cat = [list(proteinsetgo_class_updown_df.columns)[i].split(' ')[0] for i in need_list_1]
		for i in range(len(sub_cat)):
			sort_cols = [0, 1, (i+1)*3, (i+1)*3+1]
			proteinsetgo_class_ps_df_head20 = proteinsetgo_class_updown_df.sort_values(by=[proteinsetgo_class_updown_df.columns[need_list_1[i]], proteinsetgo_class_updown_df.columns[1]], ascending=[False, True]).iloc[:, sort_cols].drop_duplicates(subset=None, keep='first', inplace=False).iloc[:20, :]
			proteinsetgo_class_ps_df_head20[proteinsetgo_class_ps_df_head20.columns[3]] = proteinsetgo_class_ps_df_head20[proteinsetgo_class_ps_df_head20.columns[3]].map(lambda x: x.split('(')[0]).astype(float)
			list_list = proteinsetgo_class_ps_df_head20.values.tolist()
			self.chart_proteinsetcluster_goclass_tojs(sub_cat[i] + '.webmulti_go', list_list, "proteinsetgo_class_updown.json", y_right=True)

		if len(sub_cat) > 1:
			# proteinsetgo_class_updown_df['tmp_sum'] = 0
			# for ii in need_list_1:
			# 	proteinsetgo_class_updown_df['tmp_sum'] += proteinsetgo_class_updown_df.iloc[:,ii]
			proteinsetgo_class_updown_df_head20 = proteinsetgo_class_updown_df.sort_values(by = [proteinsetgo_class_updown_df.columns[3],proteinsetgo_class_updown_df.columns[1]], ascending=[False,True]).iloc[:,need_list_2].drop_duplicates(subset=None,keep='first',inplace=False).iloc[:20,:]
			for ii in need_list_3:
				proteinsetgo_class_updown_df_head20[proteinsetgo_class_updown_df_head20.columns[ii]] = proteinsetgo_class_updown_df_head20[proteinsetgo_class_updown_df_head20.columns[ii]].map(lambda x:x.split('(')[0]).astype(float)
			list_list = proteinsetgo_class_updown_df_head20.values.tolist()
			self.chart_proteinsetcluster_goclass_tojs('webmulti_go', list_list, "proteinsetgo_class_updown.json", sub_cat=sub_cat)

	def chart_proteinsetcluster_goclass_bothdir_web(self, proteinsetgo_class_updown_file):
		proteinsetgo_class_updown_df = pd.read_table(proteinsetgo_class_updown_file, index_col=False)
		if proteinsetgo_class_updown_df.columns[3].split(' ')[0].endswith('down') and proteinsetgo_class_updown_df.columns[6].split(' ')[0].endswith('up'):
			proteinsetgo_class_updown_df = proteinsetgo_class_updown_df.iloc[:, [0,1,2,6,7,8,3,4,5]]
		#___________双向直方图-up_down___________
		proteinsetgo_class_updown_td_df_head20_ = proteinsetgo_class_updown_df.iloc[:,[0,1,3,6]].drop_duplicates(subset=None,keep='first',inplace=False).sort_values(by = [proteinsetgo_class_updown_df.columns[3], proteinsetgo_class_updown_df.columns[6], proteinsetgo_class_updown_df.columns[1]], ascending=[False, False, True]).iloc[:20,:]
		def get_sort(x):
		  df = x.sort_values(by=[x.columns[2], x.columns[3], x.columns[1]], ascending=[False, False, True])
		  return df
		proteinsetgo_class_updown_td_df_head20 = proteinsetgo_class_updown_td_df_head20_.groupby(proteinsetgo_class_updown_td_df_head20_.columns[0],as_index=False).apply(get_sort)
		color_dic = {'biological_process':'#FF00FF', 'cellular_component':'#388E3C', 'molecular_function':'#F44336'}
		abbr_dic = {'biological_process':'BP', 'cellular_component':'CC', 'molecular_function':'MF'}
		categories = [i[1]+':'+abbr_dic[i[0]] for i in proteinsetgo_class_updown_td_df_head20.iloc[:,[0,1]].values.tolist()]
		categories_color = [{'color':color_dic[i[0]], 'group':abbr_dic[i[0]], 'value':i[1]+':'+abbr_dic[i[0]]} for i in proteinsetgo_class_updown_td_df_head20.iloc[:,[0,1]].values.tolist()]
		colors_range_ = [i['color'] for i in categories_color]
		colors_range = []
		colors_range_1 = [i for i,val in enumerate(colors_range_) if val=='#FF00FF']
		colors_range_11 = {'start':min(colors_range_1),'end':max(colors_range_1),'top_color':'#FF00FF','bottom_color':'#FF00FF'}
		if colors_range_1:
			colors_range_11 = {'start':min(colors_range_1),'end':max(colors_range_1),'top_color':'#FF00FF','bottom_color':'#FF00FF'}
			colors_range.append(colors_range_11)
		colors_range_2 = [i for i,val in enumerate(colors_range_) if val=='#388E3C']
		colors_range_22 = {'start':min(colors_range_2),'end':max(colors_range_2),'top_color':'#388E3C','bottom_color':'#388E3C'}
		if colors_range_2:
			colors_range_22 = {'start':min(colors_range_2),'end':max(colors_range_2),'top_color':'#388E3C','bottom_color':'#388E3C'}
			colors_range.append(colors_range_22)
		colors_range_3 = [i for i,val in enumerate(colors_range_) if val=='#F44336']
		colors_range_33 = {'start':min(colors_range_3),'end':max(colors_range_3),'top_color':'#F44336','bottom_color':'#F44336'}
		colors_range = [colors_range_11,colors_range_22,colors_range_33]
		if colors_range_3:
			colors_range_33 = {'start':min(colors_range_3),'end':max(colors_range_3),'top_color':'#F44336','bottom_color':'#F44336'}
			colors_range.append(colors_range_33)
		proteinsetgo_class_updown_td_df_head20.iloc[:,[3]] = -proteinsetgo_class_updown_td_df_head20.iloc[:,[3]]
		data = proteinsetgo_class_updown_td_df_head20.iloc[:,[2,3]].T.values.tolist()
		self.chart_proteinsetgo_class_updown_td_tojs('webmulti_go_neg', data, colors_range, categories, categories_color, "proteinsetgo_class_updown_td.json")

	def chart_proteinsetcluster_goclass_tojs(self, name, list_list, json_mode, sub_cat=None, y_right=False):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".go_bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = list_list
			if sub_cat:
				a["sub_categories"] = sub_cat
			if y_right:
				a['params']["yAxis2Display"] = True
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4, sort_keys=True))
			self.js_list.append([self.output_dir + name + ".go_bar.js", {"model":"highchart","highchart_type":"go_bar","width":1100,"height":600}])

	def chart_proteinsetgo_class_updown_td_tojs(self, name, data, colors_range, categories, categories_color, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar_neg.js", 'w') as fo:
			a = json.loads(f.read())
			a["categories_color"] = categories_color
			a["params"]["colors_range"] = colors_range
			a["data"] = data
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4, sort_keys=True))
			self.js_list.append([self.output_dir + name + ".bar_neg.js", {"model":"highchart","highchart_type":"showBarNegative","width":800,"height":549}])

	def chart_proteinsetcluster_gorich(self, diff_group_list, go_enrich_all_files, go_enrich_down_files, go_enrich_up_files, pset_name=None):
		for i,this_group_name in enumerate(diff_group_list):
			go_enrich_all_file = go_enrich_all_files[i]
			go_enrich_down_file = go_enrich_down_files[i]
			go_enrich_up_file = go_enrich_up_files[i]

			# Bar pdf
			go_enrich_all_file_df = pd.read_table(go_enrich_all_file)
			go_enrich_all_file_df = go_enrich_all_file_df[(go_enrich_all_file_df.iloc[:,2] == 'e')]
			go_enrich_all_file_df['new_percent'] = go_enrich_all_file_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/go_enrich_all_file_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
			go_enrich_all_data_ = go_enrich_all_file_df.sort_values(by = [go_enrich_all_file_df.columns[6],go_enrich_all_file_df.columns[3]], ascending=True).iloc[:20,[1,3,11,6]]
			go_enrich_all_data_.iloc[:,1] = go_enrich_all_data_.iloc[:,1]+':'+go_enrich_all_data_.iloc[:,0]
			go_enrich_all_data = go_enrich_all_data_.values.tolist()

			if this_group_name != '':
				self.chart_proteinsetcluster_gorichBar_tojs(this_group_name+"_all____"+this_group_name+'_go_enrich_bar', go_enrich_all_data, "proteinsetgo_richBar.json", this_group_name)

				go_enrich_down_file_df = pd.read_table(go_enrich_down_file)
				go_enrich_down_file_df = go_enrich_down_file_df[(go_enrich_down_file_df.iloc[:,2] == 'e')]
				go_enrich_down_file_df['new_percent'] = go_enrich_down_file_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/go_enrich_down_file_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
				go_enrich_down_data_ = go_enrich_down_file_df.sort_values(by = [go_enrich_down_file_df.columns[6],go_enrich_down_file_df.columns[3]], ascending=True).iloc[:20,[1,3,11,6]]
				go_enrich_down_data_.iloc[:,1] = go_enrich_down_data_.iloc[:,1]+':'+go_enrich_down_data_.iloc[:,0]
				go_enrich_down_data = go_enrich_down_data_.values.tolist()
				self.chart_proteinsetcluster_gorichBar_tojs(this_group_name+"_down____"+this_group_name+'_go_enrich_bar', go_enrich_down_data, "proteinsetgo_richBar.json", this_group_name)

				go_enrich_up_file_df = pd.read_table(go_enrich_up_file)
				go_enrich_up_file_df = go_enrich_up_file_df[(go_enrich_up_file_df.iloc[:,2] == 'e')]
				go_enrich_up_file_df['new_percent'] = go_enrich_up_file_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/go_enrich_up_file_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
				go_enrich_up_data_ = go_enrich_up_file_df.sort_values(by = [go_enrich_up_file_df.columns[6],go_enrich_up_file_df.columns[3]], ascending=True).iloc[:20,[1,3,11,6]]
				go_enrich_up_data_.iloc[:,1] = go_enrich_up_data_.iloc[:,1]+':'+go_enrich_up_data_.iloc[:,0]
				go_enrich_up_data = go_enrich_up_data_.values.tolist()
				self.chart_proteinsetcluster_gorichBar_tojs(this_group_name+"_up____"+this_group_name+'_go_enrich_bar', go_enrich_up_data, "proteinsetgo_richBar.json", this_group_name)
			else:
				self.chart_proteinsetcluster_gorichBar_tojs(this_group_name + "_all____" + this_group_name + '_go_enrich_bar', go_enrich_all_data,"proteinsetgo_richBar.json", pset_name)

			## DenseBubble pdf
			def get_sort(x):
			  df = x.sort_values(by = x.columns[6], ascending=True)
			  return [{'y':i[0],'x':i[1],'desc':i[2],'name':i[3], 'size': (math.log(i[4])+1)*3, 'realsize':i[4]} for i in df.values.tolist()]

			def bubble_each(go_enrich_file_df, data_type, group_name):
				go_enrich_file_df['logUnpvalue'] = -np.log2(go_enrich_file_df.iloc[:,6])
				go_enrich_all_data_ = go_enrich_file_df.sort_values(by = [go_enrich_file_df.columns[6],go_enrich_file_df.columns[1]], ascending=True).iloc[:20,[12,11,3,0,8,1,6]]
				# go_enrich_all_data__ = go_enrich_all_data_.sort_values(by = [go_enrich_all_data_.columns[5],go_enrich_all_data_.columns[6]], ascending=True).groupby(go_enrich_all_data_.columns[5],as_index=False).apply(get_sort)
				go_enrich_all_data__ = go_enrich_all_data_.sort_values(by=[go_enrich_all_data_.columns[5], go_enrich_all_data_.columns[6]], ascending=True)
				bp_data = get_sort(go_enrich_all_data__[go_enrich_all_data__[go_enrich_all_data__.columns[5]] == 'BP'])
				cc_data = get_sort(go_enrich_all_data__[go_enrich_all_data__[go_enrich_all_data__.columns[5]] == 'CC'])
				mf_data = get_sort(go_enrich_all_data__[go_enrich_all_data__[go_enrich_all_data__.columns[5]] == 'MF'])
				go_enrich_all_data = [bp_data, cc_data, mf_data]
				# go_enrich_all_data = go_enrich_all_data__.values.tolist()
				self.chart_proteinsetcluster_gorichDenseBubble_tojs(group_name+"_{}____".format(data_type)+group_name+'_go_enrich_bubble', go_enrich_all_data, False,'densebubble', "proteinsetgo_DenseBubble.json", group_name)
				self.chart_proteinsetcluster_gorichDenseBubble_tojs(group_name+"_{}____".format(data_type)+group_name+'_go_enrich_bubble2', go_enrich_all_data, True,'scatterbubble', "proteinsetgo_DenseBubble.json", group_name)

			if this_group_name != '':
				bubble_each(go_enrich_all_file_df, 'all', this_group_name)
				bubble_each(go_enrich_down_file_df, 'down', this_group_name)
				bubble_each(go_enrich_up_file_df, 'up', this_group_name)
			else:
				bubble_each(go_enrich_all_file_df, 'all', pset_name)

			# go_enrich_down_file_df['logUnpvalue'] = -np.log2(go_enrich_down_file_df.iloc[:,6])
			# go_enrich_down_data_ = go_enrich_down_file_df.sort_values(by = [go_enrich_down_file_df.columns[6],go_enrich_down_file_df.columns[1]], ascending=True).iloc[:20,[12,11,3,0,8,1,6]]
			# go_enrich_down_data__ = go_enrich_down_data_.sort_values(by = [go_enrich_down_data_.columns[5],go_enrich_down_data_.columns[6]], ascending=True).groupby(go_enrich_down_data_.columns[5],as_index=False).apply(get_sort)
			# go_enrich_down_data = go_enrich_down_data__.values.tolist()
			# self.chart_proteinsetcluster_gorichDenseBubble_tojs(this_group_name+"_down____"+this_group_name+'_go_enrich_bubble', go_enrich_down_data, False,'densebubble', "proteinsetgo_DenseBubble.json", this_group_name)
			# self.chart_proteinsetcluster_gorichDenseBubble_tojs(this_group_name+"_down____"+this_group_name+'_go_enrich_bubble2', go_enrich_down_data, True, 'scatterbubble',"proteinsetgo_DenseBubble.json", this_group_name)
			#
			# go_enrich_up_file_df['logUnpvalue'] = -np.log2(go_enrich_up_file_df.iloc[:,6])
			# go_enrich_up_data_ = go_enrich_up_file_df.sort_values(by = [go_enrich_up_file_df.columns[6],go_enrich_up_file_df.columns[1]], ascending=True).iloc[:20,[12,11,3,0,8,1,6]]
			# go_enrich_up_data__ = go_enrich_up_data_.sort_values(by = [go_enrich_up_data_.columns[5],go_enrich_up_data_.columns[6]], ascending=True).groupby(go_enrich_up_data_.columns[5],as_index=False).apply(get_sort)
			# go_enrich_up_data = go_enrich_up_data__.values.tolist()
			# self.chart_proteinsetcluster_gorichDenseBubble_tojs(this_group_name+"_up____"+this_group_name+'_go_enrich_bubble', go_enrich_up_data, False,'densebubble', "proteinsetgo_DenseBubble.json", this_group_name)
			# self.chart_proteinsetcluster_gorichDenseBubble_tojs(this_group_name+"_up____"+this_group_name+'_go_enrich_bubble2', go_enrich_up_data, True,'scatterbubble', "proteinsetgo_DenseBubble.json", this_group_name)


	def chart_proteinsetcluster_gorichBar_tojs(self, name, data, json_mode, group_name):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".shadowbar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["params"]["text"] = a["params"]["text"].format(group_name)
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4, sort_keys=True))
			self.js_list.append([self.output_dir + name + ".shadowbar.js", {"model":"highchart","highchart_type":"shadow_bar","width":850,"height":650}])

	def chart_proteinsetcluster_gorichDenseBubble_tojs(self, name, data, scatter_or_not,bubble_type, json_mode, group_name):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + "."+bubble_type+".js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a['params']['disperse'] = scatter_or_not
			a["params"]["title"] = a["params"]["title"].format(group_name)
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4, sort_keys=True))
			self.js_list.append([self.output_dir + name + "."+bubble_type+".js", {"model":"highchart","highchart_type":"graph.bubble","width":750,"height":550}])

	def chart_proteinsetcluster_keggclass(self, diff_group_list, kegg_class_stat_all_files, kegg_class_stat_updown_files, kegg_class_level_file):
		color_dic = {"Metabolism":"#FF9800", "Genetic Information Processing":"#F44336", "Environmental Information Processing":"#0288D1", "Cellular Processes":"#FF00FF", "Organismal Systems":"#ff0", "Human Diseases":"#388E3C"}
		kegg_level_df = pd.read_table(kegg_class_level_file)
		for index,this_group_name in enumerate(diff_group_list):
			###### all #######
			kegg_stat_all_xls = kegg_class_stat_all_files[index]
			if os.path.exists(kegg_stat_all_xls):
				kegg_stat_all_df = pd.read_table(kegg_stat_all_xls)
				kegg_stat_all_merge_df = pd.merge(kegg_level_df.iloc[:,[0,4,9]],kegg_stat_all_df.iloc[:,[0,3]],on=kegg_level_df.columns[0])
				merge_df_ = kegg_stat_all_merge_df.iloc[:,[2,1,3]].replace(np.nan, '', regex=True).values.tolist()
				tmp_dic = {}
				tmp_dic_fs = {}
				for i in merge_df_:
					tmp_dic_fs[i[0]] = i[1]
					if i[0] not in tmp_dic.keys():
						tmp_dic[i[0]] = i[2]
					else:
						tmp_dic[i[0]] += ';'+i[2]
				tmp_list = [[i,len(set(filter(None,tmp_dic[i].split(';')))),tmp_dic_fs[i]] for i in tmp_dic.keys()]
				merge_df__ = pd.DataFrame(tmp_list, columns = ['second_category','num','first_category'])
				merge_df__ = merge_df__.drop(merge_df__[merge_df__['num']==0].index)
				merge_df__['tmp_category'] = pd.Categorical(
					merge_df__.iloc[:,2], 
					categories=['Metabolism','Genetic Information Processing','Environmental Information Processing','Cellular Processes','Organismal Systems','Human Diseases'], 
					ordered=True
				)
				merge_df___sort_by_category = merge_df__.sort_values(by=['tmp_category','second_category']).iloc[:,[0,1,2]]
				categories = merge_df___sort_by_category.iloc[:,0].values.tolist()
				data = merge_df___sort_by_category.iloc[:,1].values.tolist()
				data_first_category = merge_df___sort_by_category.iloc[:,2].values.tolist()
				taxon_data = []
				data_colors = []
				for i in sorted(set(data_first_category), key=data_first_category.index):
					tmp_index = [idx for idx, e in enumerate(data_first_category) if e==i]
					taxon_data.append({
						"name":i,
						"color":color_dic[i],
						"start":tmp_index[0],
						"end":tmp_index[-1]
					})
					data_colors.extend([color_dic[i]]*len(tmp_index))
				self.chart_proteinsetcluster_keggclass_tojs(this_group_name+"_all____"+this_group_name+'_path_all', data, categories, taxon_data, data_colors, "proteinsetkegg_BarLine.json")

			###### up #######
			kegg_stat_updown_xls = kegg_class_stat_updown_files[index]
			if os.path.exists(kegg_stat_updown_xls):
				kegg_stat_updown_df = pd.read_table(kegg_stat_updown_xls)
				kegg_stat_updown_merge_df = pd.merge(kegg_level_df.iloc[:,[0,4,9]],kegg_stat_updown_df.iloc[:,[0,3]],on=kegg_level_df.columns[0])
				merge_df_ = kegg_stat_updown_merge_df.iloc[:,[2,1,3]].replace(np.nan, '', regex=True).values.tolist()
				tmp_dic = {}
				tmp_dic_fs = {}
				for i in merge_df_:
					i[2] = re.sub(r'\(.*?\)', '', i[2])
					tmp_dic_fs[i[0]] = i[1]
					if i[0] not in tmp_dic.keys():
						tmp_dic[i[0]] = i[2]
					else:
						tmp_dic[i[0]] += ';'+i[2]
				tmp_list = [[i,len(set(filter(None,tmp_dic[i].split(';')))),tmp_dic_fs[i]] for i in tmp_dic.keys()]
				merge_df__ = pd.DataFrame(tmp_list, columns = ['second_category','num','first_category'])
				merge_df__ = merge_df__.drop(merge_df__[merge_df__['num']==0].index)
				merge_df__['tmp_category'] = pd.Categorical(
					merge_df__.iloc[:,2], 
					categories=['Metabolism','Genetic Information Processing','Environmental Information Processing','Cellular Processes','Organismal Systems','Human Diseases'], 
					ordered=True
				)
				merge_df___sort_by_category = merge_df__.sort_values(by=['tmp_category','second_category']).iloc[:,[0,1,2]]
				categories = merge_df___sort_by_category.iloc[:,0].values.tolist()
				data = merge_df___sort_by_category.iloc[:,1].values.tolist()
				data_first_category = merge_df___sort_by_category.iloc[:,2].values.tolist()
				taxon_data = []
				data_colors = []
				for i in sorted(set(data_first_category), key=data_first_category.index):
					tmp_index = [idx for idx, e in enumerate(data_first_category) if e==i]
					taxon_data.append({
						"name":i,
						"color":color_dic[i],
						"start":tmp_index[0],
						"end":tmp_index[-1]
					})
					data_colors.extend([color_dic[i]]*len(tmp_index))
				self.chart_proteinsetcluster_keggclass_tojs(this_group_name+"_up_down____"+this_group_name+'_path_up', data, categories, taxon_data, data_colors, "proteinsetkegg_BarLine.json")

			###### down #######
			kegg_stat_updown_xls = kegg_class_stat_updown_files[index]
			if os.path.exists(kegg_stat_updown_xls):
				kegg_stat_updown_df = pd.read_table(kegg_stat_updown_xls)
				kegg_stat_updown_merge_df = pd.merge(kegg_level_df.iloc[:,[0,4,9]],kegg_stat_updown_df.iloc[:,[0,5]],on=kegg_level_df.columns[0])
				merge_df_ = kegg_stat_updown_merge_df.iloc[:,[2,1,3]].replace(np.nan, '', regex=True).values.tolist()
				tmp_dic = {}
				tmp_dic_fs = {}
				for i in merge_df_:
					i[2] = re.sub(r'\(.*?\)', '', i[2])
					tmp_dic_fs[i[0]] = i[1]
					if i[0] not in tmp_dic.keys():
						tmp_dic[i[0]] = i[2]
					else:
						tmp_dic[i[0]] += ';'+i[2]
				tmp_list = [[i,len(set(filter(None,tmp_dic[i].split(';')))),tmp_dic_fs[i]] for i in tmp_dic.keys()]
				merge_df__ = pd.DataFrame(tmp_list, columns = ['second_category','num','first_category'])
				merge_df__ = merge_df__.drop(merge_df__[merge_df__['num']==0].index)
				merge_df__['tmp_category'] = pd.Categorical(
					merge_df__.iloc[:,2], 
					categories=['Metabolism','Genetic Information Processing','Environmental Information Processing','Cellular Processes','Organismal Systems','Human Diseases'], 
					ordered=True
				)
				merge_df___sort_by_category = merge_df__.sort_values(by=['tmp_category','second_category']).iloc[:,[0,1,2]]
				categories = merge_df___sort_by_category.iloc[:,0].values.tolist()
				data = merge_df___sort_by_category.iloc[:,1].values.tolist()
				data_first_category = merge_df___sort_by_category.iloc[:,2].values.tolist()
				taxon_data = []
				data_colors = []
				for i in sorted(set(data_first_category), key=data_first_category.index):
					tmp_index = [idx for idx, e in enumerate(data_first_category) if e==i]
					taxon_data.append({
						"name":i,
						"color":color_dic[i],
						"start":tmp_index[0],
						"end":tmp_index[-1]
					})
					data_colors.extend([color_dic[i]]*len(tmp_index))
				self.chart_proteinsetcluster_keggclass_tojs(this_group_name+"_up_down____"+this_group_name+'_path_down', data, categories, taxon_data, data_colors, "proteinsetkegg_BarLine.json")

	def chart_proteinsetcluster_keggclass_web(self, kegg_class_stat_updown_file, kegg_class_level_file):
		color_dic = {"Metabolism":"#FF9800", "Genetic Information Processing":"#F44336", "Environmental Information Processing":"#0288D1", "Cellular Processes":"#FF00FF", "Organismal Systems":"#ff0", "Human Diseases":"#388E3C"}
		kegg_level_df = pd.read_table(kegg_class_level_file)
		kegg_stat_updown_df = pd.read_table(kegg_class_stat_updown_file)
		if int(kegg_stat_updown_df.shape[1]) < 6:
			this_group_names = [[list(kegg_stat_updown_df.columns)[2].split('_numbers')[0], 3]]
		else:
			this_group_names = [[list(kegg_stat_updown_df.columns)[2].split('_numbers')[0], 3], [list(kegg_stat_updown_df.columns)[4].split('_numbers')[0], 5]]
		for this_group_name,this_index in this_group_names:
			kegg_stat_updown_merge_df = pd.merge(kegg_level_df.iloc[:,[0,4,9]],kegg_stat_updown_df.iloc[:,[0,this_index]],on=kegg_level_df.columns[0])
			merge_df_ = kegg_stat_updown_merge_df.iloc[:,[2,1,3]].replace(np.nan, '', regex=True).values.tolist()
			tmp_dic = {}
			tmp_dic_fs = {}
			for i in merge_df_:
				i[2] = re.sub(r'\(.*?\)', '', i[2])
				tmp_dic_fs[i[0]] = i[1]
				if i[0] not in tmp_dic.keys():
					tmp_dic[i[0]] = i[2]
				else:
					tmp_dic[i[0]] += ';'+i[2]
			tmp_list = [[i,len(set(filter(None,tmp_dic[i].split(';')))),tmp_dic_fs[i]] for i in tmp_dic.keys()]
			merge_df__ = pd.DataFrame(tmp_list, columns = ['second_category','num','first_category'])
			merge_df__ = merge_df__.drop(merge_df__[merge_df__['num']==0].index)
			merge_df__['tmp_category'] = pd.Categorical(
				merge_df__.iloc[:,2],
				categories=['Metabolism','Genetic Information Processing','Environmental Information Processing','Cellular Processes','Organismal Systems','Human Diseases'],
				ordered=True
			)
			merge_df___sort_by_category = merge_df__.sort_values(by=['tmp_category','second_category']).iloc[:,[0,1,2]]
			categories = merge_df___sort_by_category.iloc[:,0].values.tolist()
			data = merge_df___sort_by_category.iloc[:,1].values.tolist()
			data_first_category = merge_df___sort_by_category.iloc[:,2].values.tolist()
			taxon_data = []
			data_colors = []
			for i in sorted(set(data_first_category), key=data_first_category.index):
				tmp_index = [idx for idx, e in enumerate(data_first_category) if e==i]
				taxon_data.append({
					"name":i,
					"color":color_dic[i],
					"start":tmp_index[0],
					"end":tmp_index[-1]
				})
				data_colors.extend([color_dic[i]]*len(tmp_index))
			self.chart_proteinsetcluster_keggclass_tojs(this_group_name+"_kegg", data, categories, taxon_data, data_colors, "proteinsetkegg_BarLine.json")

	def chart_proteinsetcluster_keggclass_tojs(self, name, data, categories, taxon_data, data_colors, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".barline.js", 'w') as fo:
			a = json.loads(f.read())
			a["data_colors"] = data_colors
			a["taxon_data"] = taxon_data
			a["data"] = [data]
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".barline.js", {"model":"highchart","highchart_type":"showBarLine","width":1000,"height":600}])

	def chart_proteinsetcluster_keggrich(self, diff_group_list, kegg_enrich_all_files, kegg_enrich_down_files, kegg_enrich_up_files):

		abbrev_dic = {'Cellular Processes':'CP','Environmental Information Processing':'EIP','Genetic Information Processing':'GIP','Human Diseases':'HD','Metabolism':'M','Organismal Systems':'OS'}
		for i,this_group_name in enumerate(diff_group_list):
			
			def get_sort(x):
			  df = x.sort_values(by = x.columns[6], ascending=True)
			  return [{'y':i[0],'x':i[1],'desc':i[2],'name':i[3],'size':i[4]} for i in df.values.tolist()]

			kegg_enrich_all_file = kegg_enrich_all_files[i]
			kegg_enrich_down_file = kegg_enrich_down_files[i]
			kegg_enrich_up_file = kegg_enrich_up_files[i]

			if os.path.exists(kegg_enrich_all_file):
				# Bar pdf
				kegg_enrich_all_file_df = pd.read_table(kegg_enrich_all_file)
				kegg_enrich_all_file_df['new_percent'] = kegg_enrich_all_file_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/kegg_enrich_all_file_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
				kegg_enrich_all_data_ = kegg_enrich_all_file_df.sort_values(by = [kegg_enrich_all_file_df.columns[6],kegg_enrich_all_file_df.columns[11]], ascending=[True,False]).iloc[:20,[11,1,12,6]]
				kegg_enrich_all_data_['typeI'] = kegg_enrich_all_data_['typeI'].map(abbrev_dic)
				kegg_enrich_all_data_.iloc[:,1] = kegg_enrich_all_data_.iloc[:,1]+':'+kegg_enrich_all_data_.iloc[:,0]
				kegg_enrich_all_data = kegg_enrich_all_data_.values.tolist()
				self.chart_proteinsetcluster_keggrichBar_tojs(this_group_name+"_all____"+this_group_name+'_kegg_enrich_bar', kegg_enrich_all_data, "proteinsetkegg_richBar.json", this_group_name)
				
				########## DenseBubble pdf

				kegg_enrich_all_file_df['logUnpvalue'] = -np.log10(kegg_enrich_all_file_df.iloc[:,6])
				kegg_enrich_all_data_ = kegg_enrich_all_file_df.sort_values(by = [kegg_enrich_all_file_df.columns[6],kegg_enrich_all_file_df.columns[11]], ascending=True).iloc[:20,[13,12,1,3,0,11,6]]
				#
				suqence_list = ['Environmental Information Processing','Genetic Information Processing','Metabolism','Cellular Processes','Organismal Systems','Human Diseases']
				kegg_enrich_all_data_['tmp_category'] = pd.Categorical(
					kegg_enrich_all_data_.iloc[:,5], 
					categories=suqence_list, 
					ordered=True
				)
				kegg_enrich_all_data____sort_by_category = kegg_enrich_all_data_.sort_values(by=['tmp_category',kegg_enrich_all_data_.columns[6]]).iloc[:,[0,1,2,3,4,5,6]]
				kegg_enrich_all_data = []
				for fc in suqence_list:
					this_fc_dicts_list = []
					for i in kegg_enrich_all_data____sort_by_category.values.tolist():
						if i[5] == fc:
							this_fc_dicts_list.append({'y':i[0],'x':i[1],'desc':i[2],'name':i[3],'size':i[4]})
					kegg_enrich_all_data.append(this_fc_dicts_list)
				#
				# kegg_enrich_all_data__ = kegg_enrich_all_data_.sort_values(by = [kegg_enrich_all_data_.columns[5],kegg_enrich_all_data_.columns[6]], ascending=True).groupby(kegg_enrich_all_data_.columns[5],as_index=False).apply(get_sort)
				# kegg_enrich_all_data = kegg_enrich_all_data__.values.tolist()
				self.chart_proteinsetcluster_keggrichDenseBubble_tojs(this_group_name+"_all____"+this_group_name+'_kegg_enrich_bubble', kegg_enrich_all_data, False,'densebubble', "proteinsetkegg_DenseBubble.json", this_group_name)
				self.chart_proteinsetcluster_keggrichDenseBubble_tojs(this_group_name+"_all____"+this_group_name+'_kegg_enrich_bubble2', kegg_enrich_all_data, True,'scatterbubble', "proteinsetkegg_DenseBubble.json", this_group_name)

			#########################################################################################################

			if os.path.exists(kegg_enrich_down_file):
				kegg_enrich_down_file_df = pd.read_table(kegg_enrich_down_file)
				kegg_enrich_down_file_df['new_percent'] = kegg_enrich_down_file_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/kegg_enrich_down_file_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
				kegg_enrich_down_data_ = kegg_enrich_down_file_df.sort_values(by = [kegg_enrich_down_file_df.columns[6],kegg_enrich_down_file_df.columns[11]], ascending=[True,False]).iloc[:20,[11,1,12,6]]
				kegg_enrich_down_data_['typeI'] = kegg_enrich_down_data_['typeI'].map(abbrev_dic)
				kegg_enrich_down_data_.iloc[:,1] = kegg_enrich_down_data_.iloc[:,1]+':'+kegg_enrich_down_data_.iloc[:,0]
				kegg_enrich_down_data = kegg_enrich_down_data_.values.tolist()
				self.chart_proteinsetcluster_keggrichBar_tojs(this_group_name+"_down____"+this_group_name+'_kegg_enrich_bar', kegg_enrich_down_data, "proteinsetkegg_richBar.json", this_group_name)

				########## DenseBubble pdf

				kegg_enrich_down_file_df['logUnpvalue'] = -np.log10(kegg_enrich_down_file_df.iloc[:,6])
				kegg_enrich_down_data_ = kegg_enrich_down_file_df.sort_values(by = [kegg_enrich_down_file_df.columns[6],kegg_enrich_down_file_df.columns[11]], ascending=True).iloc[:20,[13,12,1,3,0,11,6]]
				#
				suqence_list = ['Environmental Information Processing','Genetic Information Processing','Metabolism','Cellular Processes','Organismal Systems','Human Diseases']
				kegg_enrich_down_data_['tmp_category'] = pd.Categorical(
					kegg_enrich_down_data_.iloc[:,5], 
					categories=suqence_list, 
					ordered=True
				)
				kegg_enrich_down_data____sort_by_category = kegg_enrich_down_data_.sort_values(by=['tmp_category',kegg_enrich_down_data_.columns[6]]).iloc[:,[0,1,2,3,4,5,6]]
				kegg_enrich_down_data = []
				for fc in suqence_list:
					this_fc_dicts_list = []
					for i in kegg_enrich_down_data____sort_by_category.values.tolist():
						if i[5] == fc:
							this_fc_dicts_list.append({'y':i[0],'x':i[1],'desc':i[2],'name':i[3],'size':i[4]})
					kegg_enrich_down_data.append(this_fc_dicts_list)
				#
				# kegg_enrich_down_data__ = kegg_enrich_down_data_.sort_values(by = [kegg_enrich_down_data_.columns[5],kegg_enrich_down_data_.columns[6]], ascending=True).groupby(kegg_enrich_down_data_.columns[5],as_index=False).apply(get_sort)
				# kegg_enrich_down_data = kegg_enrich_down_data__.values.tolist()
				self.chart_proteinsetcluster_keggrichDenseBubble_tojs(this_group_name+"_down____"+this_group_name+'_kegg_enrich_bubble', kegg_enrich_down_data, False,'densebubble', "proteinsetkegg_DenseBubble.json", this_group_name)
				self.chart_proteinsetcluster_keggrichDenseBubble_tojs(this_group_name+"_down____"+this_group_name+'_kegg_enrich_bubble2', kegg_enrich_down_data, True, 'scatterbubble',"proteinsetkegg_DenseBubble.json", this_group_name)

			#########################################################################################################
			if os.path.exists(kegg_enrich_up_file):
				kegg_enrich_up_file_df = pd.read_table(kegg_enrich_up_file)
				kegg_enrich_up_file_df['new_percent'] = kegg_enrich_up_file_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/kegg_enrich_up_file_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
				kegg_enrich_up_data_ = kegg_enrich_up_file_df.sort_values(by = [kegg_enrich_up_file_df.columns[6],kegg_enrich_up_file_df.columns[11]], ascending=[True,False]).iloc[:20,[11,1,12,6]]
				kegg_enrich_up_data_['typeI'] = kegg_enrich_up_data_['typeI'].map(abbrev_dic)
				kegg_enrich_up_data_.iloc[:,1] = kegg_enrich_up_data_.iloc[:,1]+':'+kegg_enrich_up_data_.iloc[:,0]
				kegg_enrich_up_data = kegg_enrich_up_data_.values.tolist()
				self.chart_proteinsetcluster_keggrichBar_tojs(this_group_name+"_up____"+this_group_name+'_kegg_enrich_bar', kegg_enrich_up_data, "proteinsetkegg_richBar.json", this_group_name)

				########## DenseBubble pdf

				kegg_enrich_up_file_df['logUnpvalue'] = -np.log10(kegg_enrich_up_file_df.iloc[:,6])
				kegg_enrich_up_data_ = kegg_enrich_up_file_df.sort_values(by = [kegg_enrich_up_file_df.columns[6],kegg_enrich_up_file_df.columns[11]], ascending=True).iloc[:20,[13,12,1,3,0,11,6]]
				#
				suqence_list = ['Environmental Information Processing','Genetic Information Processing','Metabolism','Cellular Processes','Organismal Systems','Human Diseases']
				kegg_enrich_up_data_['tmp_category'] = pd.Categorical(
					kegg_enrich_up_data_.iloc[:,5], 
					categories=suqence_list, 
					ordered=True
				)
				kegg_enrich_up_data____sort_by_category = kegg_enrich_up_data_.sort_values(by=['tmp_category',kegg_enrich_up_data_.columns[6]]).iloc[:,[0,1,2,3,4,5,6]]
				kegg_enrich_up_data = []
				for fc in suqence_list:
					this_fc_dicts_list = []
					for i in kegg_enrich_up_data____sort_by_category.values.tolist():
						if i[5] == fc:
							this_fc_dicts_list.append({'y':i[0],'x':i[1],'desc':i[2],'name':i[3],'size':i[4]})
					kegg_enrich_up_data.append(this_fc_dicts_list)
				#
				# kegg_enrich_up_data__ = kegg_enrich_up_data_.sort_values(by = [kegg_enrich_up_data_.columns[5],kegg_enrich_up_data_.columns[6]], ascending=True).groupby(kegg_enrich_up_data_.columns[5],as_index=False).apply(get_sort)
				# kegg_enrich_up_data = kegg_enrich_up_data__.values.tolist()
				self.chart_proteinsetcluster_keggrichDenseBubble_tojs(this_group_name+"_up____"+this_group_name+'_kegg_enrich_bubble', kegg_enrich_up_data, False,'densebubble', "proteinsetkegg_DenseBubble.json", this_group_name)
				self.chart_proteinsetcluster_keggrichDenseBubble_tojs(this_group_name+"_up____"+this_group_name+'_kegg_enrich_bubble2', kegg_enrich_up_data, True,'scatterbubble', "proteinsetkegg_DenseBubble.json", this_group_name)

			#########################################################################################################

	def chart_proteinsetcluster_keggrich_web(self, kegg_enrich_all_file, group_name):

		abbrev_dic = {'Cellular Processes':'CP','Environmental Information Processing':'EIP','Genetic Information Processing':'GIP','Human Diseases':'HD','Metabolism':'M','Organismal Systems':'OS'}

		def get_sort(x):
		  df = x.sort_values(by = x.columns[6], ascending=True)
		  return [{'y':i[0],'x':i[1],'desc':i[2],'name':i[3],'size':i[4]} for i in df.values.tolist()]

		if os.path.exists(kegg_enrich_all_file):
			# Bar pdf
			kegg_enrich_all_file_df = pd.read_table(kegg_enrich_all_file)
			kegg_enrich_all_file_df['new_percent'] = kegg_enrich_all_file_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/kegg_enrich_all_file_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
			kegg_enrich_all_data_ = kegg_enrich_all_file_df.sort_values(by = [kegg_enrich_all_file_df.columns[6],kegg_enrich_all_file_df.columns[11]], ascending=[True,False]).iloc[:20,[11,1,12,6]]
			kegg_enrich_all_data_['typeI'] = kegg_enrich_all_data_['typeI'].map(abbrev_dic)
			kegg_enrich_all_data_.iloc[:,1] = kegg_enrich_all_data_.iloc[:,1]+':'+kegg_enrich_all_data_.iloc[:,0]
			kegg_enrich_all_data = kegg_enrich_all_data_.values.tolist()
			self.chart_proteinsetcluster_keggrichBar_tojs('enrichkegg', kegg_enrich_all_data, "proteinsetkegg_richBar.json", group_name)

			########## DenseBubble pdf

			kegg_enrich_all_file_df['logUnpvalue'] = -np.log10(kegg_enrich_all_file_df.iloc[:,6])
			kegg_enrich_all_data_ = kegg_enrich_all_file_df.sort_values(by = [kegg_enrich_all_file_df.columns[6],kegg_enrich_all_file_df.columns[11]], ascending=True).iloc[:20,[13,12,1,3,0,11,6]]
			#
			suqence_list = ['Environmental Information Processing','Genetic Information Processing','Metabolism','Cellular Processes','Organismal Systems','Human Diseases']
			kegg_enrich_all_data_['tmp_category'] = pd.Categorical(
				kegg_enrich_all_data_.iloc[:,5],
				categories=suqence_list,
				ordered=True
			)
			kegg_enrich_all_data____sort_by_category = kegg_enrich_all_data_.sort_values(by=['tmp_category',kegg_enrich_all_data_.columns[6]]).iloc[:,[0,1,2,3,4,5,6]]
			kegg_enrich_all_data = []
			for fc in suqence_list:
				this_fc_dicts_list = []
				for i in kegg_enrich_all_data____sort_by_category.values.tolist():
					if i[5] == fc:
						this_fc_dicts_list.append({'y':i[0],'x':i[1],'desc':i[2],'name':i[3],'size':i[4]})
				kegg_enrich_all_data.append(this_fc_dicts_list)
			#
			# kegg_enrich_all_data__ = kegg_enrich_all_data_.sort_values(by = [kegg_enrich_all_data_.columns[5],kegg_enrich_all_data_.columns[6]], ascending=True).groupby(kegg_enrich_all_data_.columns[5],as_index=False).apply(get_sort)
			# kegg_enrich_all_data = kegg_enrich_all_data__.values.tolist()
			self.chart_proteinsetcluster_keggrichDenseBubble_tojs('enrichkegg', kegg_enrich_all_data, False,'densebubble', "proteinsetkegg_DenseBubble.json", group_name)
			self.chart_proteinsetcluster_keggrichDenseBubble_tojs('enrichkegg', kegg_enrich_all_data, True,'scatterbubble', "proteinsetkegg_DenseBubble.json", group_name)

	def chart_proteinsetcluster_keggrichBar_tojs(self, name, data, json_mode, group_name):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".shadowbar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a['params']['linearGradient'] = OrderedDict((k, a['params']['linearGradient'].get(k)) for k in ["0%","40%","65%","100%"])
			a['params']["text"] = a['params']["text"].format(group_name)
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".shadowbar.js", {"model":"highchart","highchart_type":"shadow_bar","width":850,"height":550}])

	def chart_proteinsetcluster_keggrichDenseBubble_tojs(self, name, data, scatter_or_not,bubble_type, json_mode, group_name):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + "."+bubble_type+".js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a['params']['disperse'] = scatter_or_not
			a['params']["title"] = a['params']["title"].format(group_name)
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4, sort_keys=True))
			self.js_list.append([self.output_dir + name + "."+bubble_type+".js", {"model":"highchart","highchart_type":"graph.bubble","width":900,"height":550}])

	def chart_proteinset_circ(self,diff_group_list, geneset_circ_choose_files, geneset_circ_zscore_files):
		for index,this_group_name in enumerate(diff_group_list):
			for ite in [['all','go'],['all','kegg'],['down','go'],['down','kegg'],['up','go'],['up','kegg']]:
				geneset_circ_choose = geneset_circ_choose_files[index].format(all_up_down=ite[0],anno_type=ite[1])
				geneset_circ_zscore = geneset_circ_zscore_files[index].format(all_up_down=ite[0],anno_type=ite[1])
				if os.path.exists(geneset_circ_choose) and os.path.exists(geneset_circ_zscore):
					body = pd.read_table(geneset_circ_choose).values.tolist()
					pdf_height_index = pd.read_table(geneset_circ_choose).shape[0]
					try:
						header = pd.read_table(geneset_circ_zscore, header=None, index_col=0).values.tolist()
						self.chart_proteinset_circ_tojs(
							this_group_name + '_' + ite[0] + "____" + this_group_name + '_' + ite[1] + '_chord', body,
							header, pdf_height_index, "proteinsetcirc.json", ite[1])
					except:
						pass

	def chart_proteinset_circ_web(self, geneset_circ_choose, geneset_circ_zscore, annot_type):
		if os.path.exists(geneset_circ_choose) and os.path.exists(geneset_circ_zscore):
			body = pd.read_table(geneset_circ_choose).values.tolist()
			pdf_height_index = pd.read_table(geneset_circ_choose).shape[0]
			header = pd.read_table(geneset_circ_zscore, header=None, index_col=0).values.tolist()
			self.chart_proteinset_circ_tojs('chord', body, header, pdf_height_index, "proteinsetcirc.json", annot_type)

	def chart_proteinset_circ_tojs(self, name, body, header,pdf_height_index, json_mode, annot_type):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".circ.js", 'w') as fo:
			a = json.loads(f.read())
			a["table"]["body"] = body
			a["table"]["header"] = header
			if annot_type == 'go':
				a['title'] = 'GO Term'
				a['legendTitle'] = 'GO Term'
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			# 41*2*k=838*3.14, k=32, so height=(rows*2*32)/3.14=rows*20.38~~rows*21
			self.js_list.append([self.output_dir + name + ".circ.js", {"model":"highchart","highchart_type":"circ","width":55*pdf_height_index,"height":21*pdf_height_index}])

	def chart_ppi_centrality(self,diff_group_list, ppi_centrality_all_files, ppi_centrality_up_files, ppi_centrality_down_files):
		for index,this_group_name in enumerate(diff_group_list):
			ppi_centrality_all_df = pd.read_table(ppi_centrality_all_files[index], dtype={'Node_ID':object})
			categories = ppi_centrality_all_df.iloc[:,0].values.tolist()
			data = [ppi_centrality_all_df.iloc[:,2].values.tolist(), ppi_centrality_all_df.iloc[:,3].values.tolist(), ppi_centrality_all_df.iloc[:,4].values.tolist()]
			self.chart_ppi_centrality_tojs(this_group_name+"_all____"+"ppinetwork_topology____ppi.centrality.line", categories, data, "ppi_centrality.json")
			#
			ppi_centrality_up_df = pd.read_table(ppi_centrality_up_files[index], dtype={'Node_ID':object})
			categories = ppi_centrality_up_df.iloc[:,0].values.tolist()
			data = [ppi_centrality_up_df.iloc[:,2].values.tolist(), ppi_centrality_up_df.iloc[:,3].values.tolist(), ppi_centrality_up_df.iloc[:,4].values.tolist()]
			self.chart_ppi_centrality_tojs(this_group_name+"_up____"+"ppinetwork_topology____ppi.centrality.line", categories, data, "ppi_centrality.json")
			#
			ppi_centrality_down_df = pd.read_table(ppi_centrality_down_files[index], dtype={'Node_ID':object})
			categories = ppi_centrality_down_df.iloc[:,0].values.tolist()
			data = [ppi_centrality_down_df.iloc[:,2].values.tolist(), ppi_centrality_down_df.iloc[:,3].values.tolist(), ppi_centrality_down_df.iloc[:,4].values.tolist()]
			self.chart_ppi_centrality_tojs(this_group_name+"_down____"+"ppinetwork_topology____ppi.centrality.line", categories, data, "ppi_centrality.json")

	def chart_ppi_centrality_tojs(self, name, categories, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".showCurve.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".showCurve.js", {"model":"highchart","highchart_type":"showCurve","width":700,"height":430}])

	def chart_ppi_degree(self,diff_group_list, ppi_degree_all_files, ppi_degree_up_files, ppi_degree_down_files):
		for index,this_group_name in enumerate(diff_group_list):
			ppi_degree_all_df = pd.read_table(ppi_degree_all_files[index])
			categories = ppi_degree_all_df.iloc[:,0].values.tolist()
			data = ppi_degree_all_df.iloc[:,1].values.tolist()
			self.chart_ppi_degree_tojs(this_group_name+"_all____"+"ppinetwork_topology____ppi.degree.line", categories, data, "ppi_degree.json")
			#
			ppi_degree_up_df = pd.read_table(ppi_degree_up_files[index])
			categories = ppi_degree_up_df.iloc[:,0].values.tolist()
			data = ppi_degree_up_df.iloc[:,1].values.tolist()
			self.chart_ppi_degree_tojs(this_group_name+"_up____"+"ppinetwork_topology____ppi.degree.line", categories, data, "ppi_degree.json")
			#
			ppi_degree_down_df = pd.read_table(ppi_degree_down_files[index])
			categories = ppi_degree_down_df.iloc[:,0].values.tolist()
			data = ppi_degree_down_df.iloc[:,1].values.tolist()
			self.chart_ppi_degree_tojs(this_group_name+"_down____"+"ppinetwork_topology____ppi.degree.line", categories, data, "ppi_degree.json")

	def chart_ppi_degree_tojs(self, name, categories, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".showCurve.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [data]
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".showCurve.js", {"model":"highchart","highchart_type":"showCurve","width":650,"height":430}])

	def chart_ppi_centrality_degree_web(self, ppi_centrality_all_file, ppi_degree_all_file):
		ppi_centrality_all_df = pd.read_table(ppi_centrality_all_file, dtype={'Node_ID':object})
		categories = ppi_centrality_all_df.iloc[:,0].values.tolist()
		data = [ppi_centrality_all_df.iloc[:,2].values.tolist(), ppi_centrality_all_df.iloc[:,3].values.tolist(), ppi_centrality_all_df.iloc[:,4].values.tolist()]
		self.chart_ppi_centrality_tojs("ppi.centrality.line", categories, data, "ppi_centrality.json")

		ppi_degree_all_df = pd.read_table(ppi_degree_all_file)
		categories = ppi_degree_all_df.iloc[:,0].values.tolist()
		data = ppi_degree_all_df.iloc[:,1].values.tolist()
		self.chart_ppi_degree_tojs("ppi.degree.line", categories, data, "ppi_degree.json")

	def chart_cog(self,diff_group_list, cog_all_files, cog_up_down_files):
		for index,this_group_name in enumerate(diff_group_list):
			cog_all_df = pd.read_table(cog_all_files[index])
			categories = [i[1] for i in cog_all_df.iloc[:,1].values.tolist()]
			data = [cog_all_df.iloc[:,2].values.tolist()]
			legend = [this_group_name+'_all']
			self.chart_cog_tojs(this_group_name+"____"+this_group_name+'_all_cog', categories, data, legend, "proteinsetcog_class.json")
			#
			cog_up_down_df = pd.read_table(cog_up_down_files[index])
			categories = [i[1] for i in cog_up_down_df.iloc[:,1].values.tolist()]
			data = [cog_up_down_df.iloc[:,2].values.tolist(), cog_up_down_df.iloc[:,4].values.tolist()]
			legend = [this_group_name+'_up', this_group_name+'_down']
			self.chart_cog_tojs(this_group_name+"____"+this_group_name+'_up_down_cog', categories, data, legend,"proteinsetcog_class.json")

	def chart_cog_web(self, cog_all_files):
		cog_all_df = pd.read_table(cog_all_files, index_col=False)
		categories = [i[1] for i in cog_all_df.iloc[:, 1].values.tolist()]
		if cog_all_df.shape[1] < 5:
			data = [cog_all_df.iloc[:, 2].values.tolist()]
			legend = [list(cog_all_df.columns)[2][:-4]]
		else:
			data = [cog_all_df.iloc[:, 2].values.tolist(), cog_all_df.iloc[:, 4].values.tolist()]
			legend = [list(cog_all_df.columns)[2][:-4], list(cog_all_df.columns)[4][:-4]]
		self.chart_cog_tojs('cog', categories, data, legend, "proteinsetcog_class.json")

	def chart_cog_tojs(self, name, categories, data, legend, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["categories"] = categories
			a["legend"] = legend
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":700,"height":500}])

	def chart_pfam(self,diff_group_list, pfam_all_files, pfam_up_down_files):
		for index,this_group_name in enumerate(diff_group_list):
			if os.path.exists(pfam_all_files[index]):
				pfam_all_df = pd.read_table(pfam_all_files[index])
				pfam_all_df = pfam_all_df.sort_values(by = [pfam_all_df.columns[3],pfam_all_df.columns[1]], ascending=[False, True]).iloc[:20,:]
				categories = pfam_all_df.iloc[:,1].values.tolist()
				data = [pfam_all_df.iloc[:,3].values.tolist()]
				legend = [this_group_name+'_all']
				self.chart_pfam_tojs(this_group_name+"____"+this_group_name+'_all_pfam', categories, data, legend, "proteinsetpfam_class.json")
			#
			if os.path.exists(pfam_up_down_files[index]):
				pfam_up_down_df = pd.read_table(pfam_up_down_files[index])
				pfam_up_down_df = pfam_up_down_df.sort_values(by = [pfam_up_down_df.columns[3],pfam_up_down_df.columns[5],pfam_up_down_df.columns[1]], ascending=[False,False,True]).iloc[:20,:]
				categories = pfam_up_down_df.iloc[:,1].values.tolist()
				data = [pfam_up_down_df.iloc[:,3].values.tolist(), pfam_up_down_df.iloc[:,5].values.tolist()]
				legend = [this_group_name+'_up', this_group_name+'_down']
				self.chart_pfam_tojs(this_group_name+"____"+this_group_name+'_up_down_pfam', categories, data, legend,"proteinsetpfam_class.json")

	def chart_pfam_web(self, pfam_up_down_file):
		pfam_up_down_df = pd.read_table(pfam_up_down_file)
		if pfam_up_down_df.shape[1] < 6:
			pfam_up_down_df = pfam_up_down_df.sort_values(by = [pfam_up_down_df.columns[3],pfam_up_down_df.columns[1]], ascending=[False,True]).iloc[:20,:]
			categories = pfam_up_down_df.iloc[:,1].values.tolist()
			data = [pfam_up_down_df.iloc[:,3].values.tolist()]
			legend = [list(pfam_up_down_df.columns)[3].split('_num')[0]]
			self.chart_pfam_tojs('pfam', categories, data, legend, "proteinsetpfam_class.json")
			#
		if pfam_up_down_df.shape[1] > 6:
			pfam_up_down_df = pfam_up_down_df.sort_values(by = [pfam_up_down_df.columns[3],pfam_up_down_df.columns[5],pfam_up_down_df.columns[1]], ascending=[False,False,True]).iloc[:20,:]
			categories = pfam_up_down_df.iloc[:,1].values.tolist()
			data = [pfam_up_down_df.iloc[:,3].values.tolist(), pfam_up_down_df.iloc[:,5].values.tolist()]
			legend = [list(pfam_up_down_df.columns)[3].split('_num')[0], list(pfam_up_down_df.columns)[5].split('_num')[0]]
			self.chart_pfam_tojs('pfam', categories, data, legend,"proteinsetpfam_class.json")

	def chart_pfam_tojs(self, name, categories, data, legend, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["categories"] = categories
			a["legend"] = legend
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":700,"height":500}])

	def chart_subloc(self,diff_group_list, subloc_all_files, subloc_up_down_files):
		for index,this_group_name in enumerate(diff_group_list):
			subloc_all_df = pd.read_table(subloc_all_files[index])
			categories = subloc_all_df.iloc[:,0].values.tolist()
			data = [subloc_all_df.iloc[:,1].values.tolist()]
			legend = [this_group_name+'_all']
			self.chart_subloc_tojs(this_group_name+"____"+this_group_name+'_all_subloc', categories, data, legend, "proteinsetsubloc_class.json")
			#
			subloc_up_down_df = pd.read_table(subloc_up_down_files[index])
			categories = subloc_up_down_df.iloc[:,0].values.tolist()
			data = [subloc_up_down_df.iloc[:,1].values.tolist(), subloc_up_down_df.iloc[:,3].values.tolist()]
			legend = [this_group_name+'_up', this_group_name+'_down']
			self.chart_subloc_tojs(this_group_name+"____"+this_group_name+'_up_down_subloc', categories, data, legend,"proteinsetsubloc_class.json")

	def chart_subloc_web(self, subloc_up_down_file):
		subloc_up_down_df = pd.read_table(subloc_up_down_file)
		if subloc_up_down_df.shape[1] < 4:
			categories = subloc_up_down_df.iloc[:,0].values.tolist()
			data = [subloc_up_down_df.iloc[:,1].values.tolist()]
			legend = [list(subloc_up_down_df.columns)[1].split('_num')[0]]
			self.chart_subloc_tojs('subloc', categories, data, legend, "proteinsetsubloc_class.json")
			#
		if subloc_up_down_df.shape[1] > 4:
			categories = subloc_up_down_df.iloc[:,0].values.tolist()
			data = [subloc_up_down_df.iloc[:,1].values.tolist(), subloc_up_down_df.iloc[:,3].values.tolist()]
			legend = [list(subloc_up_down_df.columns)[1].split('_num')[0], list(subloc_up_down_df.columns)[3].split('_num')[0]]
			self.chart_subloc_tojs('subloc', categories, data, legend,"proteinsetsubloc_class.json")

	def chart_subloc_tojs(self, name, categories, data, legend, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["categories"] = categories
			a["legend"] = legend
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bar.js", {"model":"highchart","highchart_type":"showBar","width":700,"height":500}])

	# 下面画的两个图所在的页面 位于 样本比较分析（样本相关性热图，pca分析），这个分析模块是可选的
	# def chart_sample_corr(self, sample_corr_matrix_file, sample_corr_tree_file):
	# 	this_sample_tree_file_str = open(sample_corr_tree_file, 'r').read()
	# 	this_tree_2 = this_sample_tree_file_str.split('\n')[0]
	# 	this_tree_1 = this_tree_2
	# 	this_tree_columns = this_sample_tree_file_str.split('\n')[1].strip().split(';')
	# 	this_tree_rows = this_tree_columns
	#
	# 	this_expression_matrix_df = pd.read_table(sample_corr_matrix_file, index_col=0)
	# 	heatmap_data = this_expression_matrix_df[this_expression_matrix_df.index.isin(this_tree_rows)].reindex(this_tree_rows)[this_tree_rows].values.tolist()
	# 	self.chart_sample_corr_tojs("sam_corr", this_tree_1, this_tree_2, this_tree_rows, this_tree_columns, heatmap_data, "express_corr.json")
	#
	# def chart_sample_corr_tojs(self, name, this_tree_1, this_tree_2, this_tree_rows, this_tree_columns, heatmap_data, json_mode):
	# 	json_mode = self.mode_dir + "/" + json_mode
	# 	with open(json_mode, 'r') as f, open(self.output_dir + name + ".heatmap.js", 'w') as fo:
	# 		a = json.loads(f.read())
	# 		a["tree_2"] = this_tree_2
	# 		a["tree_1"] = this_tree_1
	# 		a["rows"] = this_tree_rows
	# 		a["columns"] = this_tree_columns
	# 		a["heatmap_data"] = heatmap_data
	# 		fo.write("var options = ")
	# 		fo.write(json.dumps(a, indent=4))
	# 		self.js_list.append([self.output_dir + name + ".heatmap.js", {"model":"highchart","highchart_type":"tree_heatmap","width":800,"height":800}])

	def chart_exp_corr(self, sample_corr_file, sample_tree, group_dict=None):
		with open(sample_tree, 'r') as f:
			sample_tree = f.readline().strip()
			samples = f.readline().strip().split(";")

		print samples
		corr_result = pd.read_table(sample_corr_file, index_col=0, header=0)
		corr_result = corr_result.loc[samples, :]
		corr_source = [[""] + samples]
		for line_dict in corr_result.iterrows():
			corr_source.append(
				[line_dict[0]] + [line_dict[1][s] for s in samples]
			)
		tree_source = sample_tree

		sample2group_source = [["name", "group"]]
		if group_dict:
			group_dict = json.loads(group_dict)
			for g, ss in group_dict.items():
				for s in ss:
					sample2group_source.append([s, g])
		else:
			for s in samples:
				sample2group_source.append([s, s])

		self.chart_heat_tree("exp", ".heatmap", corr_source, sample_tree, sample_tree, sample2group_source,
							 sample2group_source, "exp.relation.heat_tree.json")
		self.chart_corr_tree('exp', '.tree', sample_tree, sample2group_source, 'exp.relation.tree.json')

	def chart_heat_tree(self, name, out, corr_source, sample_tree, gene_tree, sample2group_source, gene2group_source,
						json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + out + ".heat_corr.js", 'w') as fo:
			a = json.loads(f.read())
			a["dataset"][0]["source"] = corr_source
			a["dataset"][1]["source"] = gene2group_source
			if gene_tree:
				a["dataset"][1]["categories"] = gene_tree
			else:
				del a["dataset"][1]["categories"]
			a["dataset"][2]["source"] = sample2group_source
			if sample_tree:
				a["dataset"][2]["categories"] = sample_tree
			else:
				del a["dataset"][2]["categories"]
				a["series"] = a["series"][:2]
				a["dataset"] = a["dataset"][:2]
				a["legend"] = a["legend"][:2]

			if len(gene2group_source) <= 1:
				##无基因聚类
				if sample_tree:
					a["dataset"] = [a["dataset"][0], a["dataset"][2]]
					a["legend"] = [a["legend"][0], a["legend"][2]]
					a["series"] = [a["series"][0], a['series'][2]]
					a["series"][1]["datasetIndex"] = 1
					a["series"][1]["visualMap"][0]["legendIndex"] = 1
					a["xAxis"][0]["length"] = "95%"
				else:
					a["dataset"] = [a["dataset"][0]]
					a["legend"] = [a["legend"][0]]
					a["series"] = [a["series"][0]]
					a["xAxis"][0]["length"] = "95%"

			# 类型待确定
			a["chart"]["type"] = "heatmap_tree"
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + out + ".heat_corr.js", {}])

	def chart_corr_tree(self, name, out, sample_tree, sample2group_source, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + out + ".corr_tree.js", 'w') as fo:
			a = json.loads(f.read())
			a['dataset'][0]['source'] = sample2group_source
			a['dataset'][0]["categories"] = sample_tree
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + out + ".corr_tree.js", {}])

	def chart_sample_pca(self, sample_pca_sites_file, sample_pca_propor_file, group_table):
		this_pca_sites_df = pd.read_table(sample_pca_sites_file)
		shape = ["circle", "triangle-up", "diamond", "square", "cross", "triangle-down"]
		color = ["#388E3C", "#F44336", "#0288D1", "#FF9800", "#87CEEB", '#E91E63', '#673AB7', '#006400', '#FFA500']
		# shape = ["circle", "triangle-up", "diamond", "square", "del", "plus", "triangle-down", "circle", "triangle", "diamond", "square", "del", "plus", "triangle-down", "circle", "triangle"]
		# color = ["#388E3C", "#F44336", "#0288D1", "#FF9800", "#C0C0C0","#808080","#000000","#800000","#FFFF00","#808000","#00FF00","#00FFFF","#008080","#000080","#FF00FF","#800080"]
		group_dict = dict()
		with open(group_table, 'r') as g:
			g.readline()
			for line in g:
				items = line.strip().split('\t')
				group_dict[items[0]] = items[1]
		sites = this_pca_sites_df.iloc[:, 0].values.tolist()
		site_list = list(set(sites))
		site_list.sort(key=sites.index)
		# shape_color_dic = {category:[shape[index],color[index]] for index,category in  enumerate(site_list)}
		shape_color_dic = dict()
		for index, category in enumerate(list(set(group_dict.values()))):
			if index < len(shape):
				shape_color_dic[category] = [shape[index], color[index]]
			elif index < len(color):
				shape_color_dic[category] = [shape[index-len(shape)], color[index]]
			else:
				shape_color_dic[category] = [shape[index - len(shape)], color[index - len(color)]]
		data = [{"name":i[0], "symbol_category_name":group_dict[i[0]], "color_category_name":group_dict[i[0]], "color":shape_color_dic[group_dict[i[0]]][1], "symbol":shape_color_dic[group_dict[i[0]]][0], "value":[i[1],i[2]]} for index,i in enumerate(this_pca_sites_df.iloc[:,[0,1,2]].values.tolist())]
		this_pca_propor_df = pd.read_table(sample_pca_propor_file)
		x_y_label = [i[0]+'('+str('%.2f%%' % (i[1] * 100))+')' for i in this_pca_propor_df.values.tolist()]
		self.chart_sample_pca_tojs("sam_pca", data, x_y_label, "express_pca.json")

	def chart_sample_pca_tojs(self, name, data, x_y_label, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".scatter.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["params"]["x_label"] = x_y_label[0]
			a["params"]["y_label"] = x_y_label[1]
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".scatter.js", {"model":"highchart","highchart_type":"scatter","width":800,"height":600}])

	# 4 dia & labelfree
	def chart_exp_venn(self, exp_venn_path):
		colors = ["#FFA500", "#388E3C", "#33CCFF", "#673AB7", "#FF2020", "#FF3366"]
		venn_df = pd.read_table(exp_venn_path)
		num_row = venn_df.shape[0]
		if 1 < num_row <= 6:
			df_to_venn = venn_df
		elif num_row > 6:
			df_to_venn = venn_df.iloc[:6, :]
		else:
			df_to_venn = pd.DataFrame()
		if not df_to_venn.empty:
			venn_color = colors[:num_row]
			all_data = [{"name": i[0], "data": i[1].split(',')} for i in df_to_venn.values.tolist()]
			self.chart_proteinset_venn_tojs("exp_venn", all_data, venn_color, "proteinset_venn.json")

	# 这个venn图 页面以前都没有。 这里根据diff文件解析proteinset，大于1并小于6时则画venn，蛋白集数目大于6则挑前6个画venn
	def chart_proteinset_venn(self, diff_path_for_proteinset):
		colors = ["#FFA500", "#388E3C", "#33CCFF", "#673AB7", "#FF2020", "#FF3366"]
		diff_list = sorted([f for f in os.listdir(diff_path_for_proteinset) if f.startswith('diffcmp')])
		if 1 < len(diff_list) <= 6:
			diff_list_to_venn = diff_list
		elif len(diff_list) > 6:
			diff_list_to_venn = diff_list[:6]
		else:
			diff_list_to_venn = []
		if diff_list_to_venn:
			venn_color = colors[:len(diff_list_to_venn)]
			all_data = []
			for diff_file in [os.path.join(diff_path_for_proteinset,i) for i in diff_list_to_venn]:
				this_diff_df = pd.read_table(diff_file)
				this_diff_df = this_diff_df[(this_diff_df.significant.isin(['yes']))&(this_diff_df.regulate.isin(['up','down']))]
				this_proteiset_name = this_diff_df.compare.values.tolist()[0].replace("|", "_vs_")+"_all"
				this_protein_names = this_diff_df.iloc[:,0].values.tolist()
				all_data.append({"name":this_proteiset_name , "data":this_protein_names})
			self.chart_proteinset_venn_tojs("venn", all_data, venn_color, "proteinset_venn.json")
		else:
			print "当前项目的蛋白集不足两个，因此没有足够数据来绘制venn图。"

	def chart_proteinset_venn_web(self, pset_venn_list):
		colors = ["#FFA500", "#388E3C", "#33CCFF", "#673AB7", "#FF2020", "#FF3366"]
		num_row = len(pset_venn_list)
		if 1 < num_row <= 6:
			list_to_venn = pset_venn_list
		elif num_row > 6:
			list_to_venn = pset_venn_list[:6]
		else:
			list_to_venn = []
		if list_to_venn:
			venn_color = colors[:num_row]
			all_data = [{"name": i[0], "data": i[1]} for i in list_to_venn]
			self.chart_proteinset_venn_tojs("venn", all_data, venn_color, "proteinset_venn.json")

	def chart_proteinset_venn_tojs(self, name, all_data, venn_color, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".venn.js", 'w') as fo:
			a = json.loads(f.read())
			a["params"]["color"] = venn_color
			a["data"] = all_data
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".venn.js", {"model":"highchart","highchart_type":"venn","width":500,"height":695}])

	# 下面是 蛋白转录联合分析 交互页面的 静态图
	def chart_protein_transcript_overview(self, protein_transcript_overview_file):
		protein_transcript_overview_str = open(protein_transcript_overview_file, 'r').read()
		tmp_l = protein_transcript_overview_str.split('\n')
		if len(tmp_l[0][7:].strip())>0:
			both = len(tmp_l[0].split(';'))
		else:
			both = 0
		protein = len(tmp_l[1].split(';')) + both
		gene = len(tmp_l[2].split(';')) + both
		# fix venn gene num showed is different from transcriptome data ,which caused by one gene correspond to multiple proteins. By xuxi on 20210621.
		try:
			gene_list_tmp = list(map(list, zip(*[i.split('|') for i in tmp_l[0].strip().split()[1].split(';')])))[1]
			gene -= (len(gene_list_tmp)-len(set(gene_list_tmp)))
		except:
			pass
		#
		self.chart_protein_transcript_overview_tojs("venn", both, protein, gene, "protein_transcript_overview.json")

	def chart_protein_transcript_overview_tojs(self, name, both, protein, gene, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".venn2.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [{"color": "#00CC66", "name": "gene", "value": gene, "sets": ["gene"]}, {"color": "#FF0033", "name": "protein", "value": protein, "sets": ["protein"]}, {"name": "both", "value": both, "sets": ["gene", "protein"]}]
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".venn2.js", {"model":"highchart","highchart_type":"venn","width":700,"height":600}])

	def chart_protein_transcript_baseinfo(self, protein_transcript_baseinfo_file):
		protein_transcript_baseinfo_str = open(protein_transcript_baseinfo_file, 'r').read()
		tmp_l = protein_transcript_baseinfo_str.split('\n')
		if len(tmp_l[0][7:].strip())>0:
			both = len(tmp_l[0].split(';'))
		else:
			both = 0
		protein = len(tmp_l[1].split(';')) + both
		gene = len(tmp_l[2].split(';')) + both
		# fix venn gene num showed is different from transcriptome data ,which caused by one gene correspond to multiple proteins. By xuxi on 20210621.
		try:
			gene_list_tmp = list(map(list, zip(*[i.split('|') for i in tmp_l[0].strip().split()[1].split(';')])))[1]
			gene -= (len(gene_list_tmp)-len(set(gene_list_tmp)))
		except:
			pass
		#
		self.chart_protein_transcript_baseinfo_tojs("protein_transcript_baseinfo", both, protein, gene, "protein_transcript_baseinfo.json")

	def chart_protein_transcript_baseinfo_tojs(self, name, both, protein, gene, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".venn2.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [{"color": "#00CC66", "name": "gene", "value": gene, "sets": ["gene"]}, {"color": "#FF0033", "name": "protein", "value": protein, "sets": ["protein"]}, {"name": "both", "value": both, "sets": ["gene", "protein"]}]
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".venn2.js", {"model":"highchart","highchart_type":"venn","width":700,"height":600}])

	def chart_protein_transcript_diff_allvenn(self, protein_transcript_diff_allvenn_file):
		protein_transcript_diff_allvenn_str = open(protein_transcript_diff_allvenn_file, 'r').read()
		tmp_l = protein_transcript_diff_allvenn_str.split('\n')
		both = len(tmp_l[0].split(';'))
		protein = len(tmp_l[1].split(';')) + len(tmp_l[0].split(';'))
		gene = len(tmp_l[2].split(';')) + len(tmp_l[0].split(';'))
		self.chart_protein_transcript_diff_allvenn_tojs("protein_transcript_diff_allvenn", both, protein, gene, "protein_transcript_diff_allvenn.json")

	def chart_protein_transcript_diff_allvenn_tojs(self, name, both, protein, gene, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".venn2.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [{"color": "#00CC66", "name": "gene", "value": gene, "sets": ["gene"]}, {"color": "#FF0033", "name": "protein", "value": protein, "sets": ["protein"]}, {"name": "both", "value": both, "sets": ["gene", "protein"]}]
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".venn2.js", {"model":"highchart","highchart_type":"venn","width":700,"height":600}])

	def chart_protein_transcript_diff_relatevenn(self, protein_transcript_diff_relatevenn_file):
		protein_transcript_diff_relatevenn_str = open(protein_transcript_diff_relatevenn_file, 'r').read()
		tmp_l = protein_transcript_diff_relatevenn_str.split('\n')
		def return_list(tmp_l,i):
			if len(tmp_l[i].strip().split('\t')) >1:
				return tmp_l[i].strip().split('\t')[1].split(';')
			else:
				return []
		protein_up = return_list(tmp_l,0)
		protein_down = return_list(tmp_l,1)
		transcript_up = return_list(tmp_l,2)
		transcript_down = return_list(tmp_l,3)
		data = [{"name": "DEG UP", "data":transcript_up},{"name": "DEG DOWN", "data":transcript_down},{"name": "DEP UP", "data":protein_up},{"name": "DEP DOWN", "data":protein_down}]
		self.chart_protein_transcript_diff_relatevenn_tojs("protein_transcript_diff_relatevenn", data, "protein_transcript_diff_relatevenn.json")

	def chart_protein_transcript_diff_relatevenn_tojs(self, name, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".venn.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".venn.js", {"model":"highchart","highchart_type":"venn","width":500,"height":695}])

	def chart_protein_transcript_diff_relatenine(self, protein_transcript_diff_relatenine_file,relatenine_params_from_mongo_maintable_params):
		result_pd = pd.read_table(protein_transcript_diff_relatenine_file, index_col=False)
		result_pd = result_pd.iloc[:,[0,1,2,3,4]].fillna('_')
		corr_value, corr_p = pearsonr(result_pd['protein_fc'], result_pd['transcript_fc'])
		tmp_param = relatenine_params_from_mongo_maintable_params.split(';')
		y_values = [-abs(np.log2(float(tmp_param[1]))),abs(np.log2(float(tmp_param[0])))]
		x_values = [-abs(np.log2(float(tmp_param[3]))),abs(np.log2(float(tmp_param[2])))]
		color_dic = {"no change on gene level-no change on protein level":"#212121","up on gene level-no change on protein level":"#006400","down on gene level-no change on protein level":"#F44336","no change on gene level-up on protein level":"#E91E63","down on gene level-up on protein level":"#0288D1","up on gene level-up on protein level":"#FFA500","no change on gene level-down on protein level":"#FF9800","up on gene level-down on protein level":"#673AB7","down on gene level-down on protein level":"#388E3C"}
		data = [{'name':i[0], 'value':[i[2],i[1]], 'symbol_category_name':i[4]+" on gene level-"+i[3]+" on protein level", 'color_category_name':i[4]+" on gene level-"+i[3]+" on protein level", 'color':color_dic[i[4]+" on gene level-"+i[3]+" on protein level"]} for i in result_pd.values.tolist()]
		self.chart_protein_transcript_diff_relatenine_tojs("protein_transcript_diff_relatenine", data, corr_value, corr_p, x_values, y_values, "protein_transcript_diff_relatenine.json")

	def chart_protein_transcript_diff_relatenine_tojs(self, name, data, corr_value, corr_p, x_values, y_values, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".nine.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["params"]["corr_value"] = corr_value
			a["params"]["corr_p"] = corr_p
			a["params"]["x_values"] = x_values
			a["params"]["y_values"] = y_values
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".nine.js", {"model":"highchart","highchart_type":"venn","width":750,"height":600}])

	def chart_pt_proteinsetcluster(self, expression_matrix_file, seq_tree_file):
		#______子聚类趋势图________
		all_sub = {}
		this_group_dir = os.path.dirname(expression_matrix_file)
		sub_files_path = []
		for file in os.listdir(this_group_dir):
			if file.startswith('seq.subcluster'):
				sub_files_path.append(os.path.join(this_group_dir,file))
		for sub_file in sub_files_path:
			tmp_n = os.path.basename(sub_file)[4:-4].split('_')
			title = tmp_n[0]+'_'+tmp_n[1]+'('+tmp_n[2]+')'
			this_sub_file_df = pd.read_table(sub_file, index_col=0)
			for i in this_sub_file_df.index:
				all_sub[i] = tmp_n[1]
			this_sub_file_df.loc['mean'] = this_sub_file_df.mean()
			data = this_sub_file_df.values.tolist()

			symbols = [""]*int(tmp_n[2])+['circle']
			colors = ["#F0F0F0"]*int(tmp_n[2])+['#000099']
			legend = this_sub_file_df.index.tolist()
			categories = [i.replace('_protein','').replace('_rna','_gene') for i in this_sub_file_df.columns.tolist()]
			self.chart_pt_proteinsetcluster_sub_tojs("pt_proteinsetcluster_"+'_'+tmp_n[0]+'_'+tmp_n[1], data, title, symbols, colors, legend, categories, "pt_proteinsetcluster_SubLine.json")

		seq_tree_file_str = open(seq_tree_file, 'r').read()
		this_tree = seq_tree_file_str.split('\n')[0]
		this_tree_rows = seq_tree_file_str.split('\n')[1].strip().split(';')
		def randomcolor():
			colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
			color = ""
			for i in range(6):
				color += colorArr[random.randint(0,14)]
			return "#"+color
		color_dic = {}
		for v in set(all_sub.values()):
			color_dic[v] = randomcolor()
		left_group_colors = [color_dic[all_sub[i]] for i in this_tree_rows]

		expression_matrix_df = pd.read_table(expression_matrix_file, index_col=0)
		height = int(int(expression_matrix_df.shape[0])*12.82+240)
		this_tree_columns = [i.replace('_protein','').replace('_rna','_gene') for i in list(expression_matrix_df.columns)]
		heatmap_data = expression_matrix_df[expression_matrix_df.index.isin(this_tree_rows)].reindex(this_tree_rows).values.tolist()
		self.chart_pt_proteinsetcluster_tojs("pt_proteinsetcluster", this_tree, this_tree_rows, left_group_colors, this_tree_columns, heatmap_data, height, "pt_proteinsetcluster_heatmap.json")
				

	def chart_pt_proteinsetcluster_tojs(self, name, this_tree, this_tree_rows, left_group_colors,this_tree_columns, heatmap_data, height, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".heatmap_new.js", 'w') as fo:
			a = json.loads(f.read())
			a["tree_1"] = this_tree
			a["rows"] = this_tree_rows
			a["params"]["left_group_colors"] = left_group_colors
			a["columns"] = this_tree_columns
			a["heatmap_data"] = heatmap_data
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".heatmap_new.js", {"model":"highchart","highchart_type":"tree_heatmap_new","width":950,"height":height}])

	def chart_pt_proteinsetcluster_sub_tojs(self, name, data, title, symbols, colors, legend, categories, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".line.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["params"]['title'] = title
			a["params"]["symbols"] = symbols
			a["params"]["colors"] = colors
			a["legend"] = legend
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".line.js", {"model":"highchart","highchart_type":"showCurve","width":600,"height":400}])

	def chart_pt_proteinsetcluster_corr(self, pt_corrScatter_file, pt_corrScatter_pvalue_file):
		pt_corrScatter_df = pd.read_table(pt_corrScatter_file)
		with open(pt_corrScatter_pvalue_file) as cor_r:
			corr_value = float(cor_r.readline().strip().split('\t')[1])
			p_value = float(cor_r.readline().strip().split('\t')[1])
		pt_corrScatter_list = pt_corrScatter_df.values.tolist()
		data = [{"name":i[0], "value":[i[2],i[1]], "density":i[3]} for i in pt_corrScatter_list]
		self.chart_pt_proteinsetcluster_corr_tojs("pt_proteinsetcluster_corr", data, corr_value, p_value, "pt_proteinsetcluster_corrScatter.json")

	def chart_pt_proteinsetcluster_corr_tojs(self, name, data, corr_value, p_value, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".scatter_new.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["params"]["p_value"] = p_value
			a["params"]["corr_value"] = corr_value
			a['params']['linearGradient'] = OrderedDict((k, a['params']['linearGradient'].get(k)) for k in ["0%","40%","65%","100%"])
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".scatter_new.js", {"model":"highchart","highchart_type":"scatter","width":750,"height":600}])

	def chart_pt_goclass(self, pt_goclass_file):
		pt_goclass_df = pd.read_table(pt_goclass_file)
		all_data = pt_goclass_df.iloc[:,[0,1,7,8,4,5]]
		data = all_data.sort_values(by = all_data.columns[5],ascending=False).iloc[:10,:].values.tolist()
		self.chart_pt_goclass_tojs("pt_goclass", data, "pt_goclass.json")

		data_pie = all_data.sort_values(by = all_data.columns[5],ascending=False).iloc[:10,[1,4,2]].values.tolist()
		data_pie_ = [{'term':i[0], 'protein_num':i[1], 'gene_num':i[2]} for i in data_pie]
		self.chart_pt_goclass_pie_tojs("pt_goclass", data_pie_, "pt_goclass_pie.json")		

	def chart_pt_goclass_tojs(self, name, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".go_bar.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["colors"] = OrderedDict([("biological_process","#FF00FF"),("cellular_component","#388E3C"),("molecular_function","#F44336")])
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".go_bar.js", {"model":"highchart","highchart_type":"go_bar","width":850,"height":550}])

	def chart_pt_goclass_pie_tojs(self, name, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".double_nest.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".double_nest.js", {"model":"highchart","highchart_type":"double_nest","width":850,"height":550}])

	def chart_pt_goenrich(self, pt_goenrich_dir):
		def delete_col(a_df,a_colname):
			if a_colname in a_df.columns.values.tolist():
				a_df.drop(a_colname, axis=1,inplace=True)
			for colname in a_df.columns.values.tolist():
				if colname.startswith(a_colname):
					a_df.drop(colname, axis=1,inplace=True)
		gene_xls_file = pt_goenrich_dir
		protein_xls_file = pt_goenrich_dir
		for file in os.listdir(pt_goenrich_dir):
			if "_protein.xls" in file:
				protein_xls_file = os.path.join(pt_goenrich_dir,file)
			elif "_gene.xls" in file:
				gene_xls_file = os.path.join(pt_goenrich_dir,file)
		if gene_xls_file != pt_goenrich_dir and protein_xls_file != pt_goenrich_dir:
			protein_xls_df = pd.read_table(protein_xls_file)
			delete_col(protein_xls_df,"p_sm")
			protein_xls_df = protein_xls_df[(protein_xls_df.iloc[:,2] == 'e')]
			protein_xls_df['enrich_factor'] = protein_xls_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/protein_xls_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
			protein_xls_data_ = protein_xls_df.sort_values(by = [protein_xls_df.columns[6],protein_xls_df.columns[3]], ascending=True).iloc[:20,[0,1,3,11,10]].values.tolist()
			
			gene_xls_df = pd.read_table(gene_xls_file)
			delete_col(gene_xls_df,"p_sm")
			gene_xls_df = gene_xls_df[(gene_xls_df.iloc[:,2] == 'e')]
			gene_xls_df['enrich_factor'] = gene_xls_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/gene_xls_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
			gene_xls_data_ = gene_xls_df.iloc[:,[0,11,10]].values.tolist()
			gene_xls_data_dic = {i[0]:i[-2:] for i in gene_xls_data_}
			data = []
			for i in protein_xls_data_:
				if i[0] in gene_xls_data_dic.keys():
					data.append(i[1:]+gene_xls_data_dic[i[0]])
				else:
					data.append(i[1:]+[0,1])
			self.chart_pt_goenrich_tojs("pt_goenrich", data, "pt_goenrich.json")

			data_protein = []
			data_gene = []
			for i in protein_xls_data_:
				for ii in protein_xls_df.iloc[:,[-2,11,3,8,0]].values.tolist():
					if i[0] == ii[4]:
						data_protein.append({'color':ii[0], 'x':ii[1], 'name':ii[2], 'size':ii[3]})
				for iii in gene_xls_df.iloc[:,[-2,11,3,8,0]].values.tolist():
					if i[0] == iii[4]:
						data_gene.append({'color':iii[0], 'x':iii[1], 'name':iii[2], 'size':iii[3]})
			data = [data_protein, data_gene]
			self.chart_pt_goenrich_bubble_tojs("pt_goenrich_bubble", data, "pt_goenrich_bubble.json")

	def chart_pt_goenrich_tojs(self, name, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".shadow_bar_groups.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a['params']["linearGradient"] = OrderedDict([("0%","#66CCFF"),("50%","#FFFFFF"), ("100%","#FF66CC")])
			a['params']["linearGradient2"] = OrderedDict([("0%","#FF9966"),("50%","#CCCCCC"), ("100%","#66FFFF")])
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".shadow_bar_groups.js", {"model":"highchart","highchart_type":"shadow_bar_groups","width":850,"height":850}])

	def chart_pt_goenrich_bubble_tojs(self, name, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bubble_groups.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bubble_groups.js", {"model":"highchart","highchart_type":"bubble_groups","width":1300,"height":650}])

	def chart_pt_goenrich_cluster(self, pt_goenrich_cluster_matrix_file, pt_goenrich_cluster_tree_file):
		#______子聚类趋势图________
		all_sub = {}
		this_group_dir = os.path.dirname(pt_goenrich_cluster_matrix_file)
		sub_files_path = []
		for file in os.listdir(this_group_dir):
			if file.startswith('seq.subcluster'):
				sub_files_path.append(os.path.join(this_group_dir,file))
		for sub_file in sub_files_path:
			tmp_n = os.path.basename(sub_file)[4:-4].split('_')
			title = tmp_n[0]+'_'+tmp_n[1]+'('+tmp_n[2]+')'
			this_sub_file_df = pd.read_table(sub_file, index_col=0)
			for i in this_sub_file_df.index:
				all_sub[i] = tmp_n[1]
			this_sub_file_df.loc['mean'] = this_sub_file_df.mean()
			data = this_sub_file_df.values.tolist()

			symbols = [""]*int(tmp_n[2])+['circle']
			colors = ["#F0F0F0"]*int(tmp_n[2])+['#000099']
			legend = this_sub_file_df.index.tolist()
			categories = [i.replace('_protein','').replace('_rna','_gene') for i in this_sub_file_df.columns.tolist()]
			self.chart_pt_goenrich_cluster_sub_tojs("pt_goenrich_cluster_"+'_'+tmp_n[0]+'_'+tmp_n[1], data, title, symbols, colors, legend, categories, "pt_proteinsetcluster_SubLine.json")

		pt_goenrich_cluster_tree_file_str = open(pt_goenrich_cluster_tree_file, 'r').read()
		this_tree = pt_goenrich_cluster_tree_file_str.split('\n')[0].replace('GO-','GO:')
		this_tree_rows = pt_goenrich_cluster_tree_file_str.split('\n')[1].strip().replace('GO-','GO:').split(';')
		def randomcolor():
			colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
			color = ""
			for i in range(6):
				color += colorArr[random.randint(0,14)]
			return "#"+color
		color_dic = {}
		for v in set(all_sub.values()):
			color_dic[v] = randomcolor()
		left_group_colors = [color_dic[all_sub[i]] for i in this_tree_rows]

		expression_matrix_df = pd.read_table(pt_goenrich_cluster_matrix_file, index_col=0)
		height = int(int(expression_matrix_df.shape[0])*8.82+240)
		# this_tree_columns = list(expression_matrix_df.columns[:2])
		# heatmap_data_ = expression_matrix_df[expression_matrix_df.index.isin(this_tree_rows)].reindex(this_tree_rows).iloc[:,[0,1,4]]
		# heatmap_data = heatmap_data_.iloc[:,[0,1]].values.tolist()
		# this_tree_rows____ = heatmap_data_.iloc[:,2].values.tolist()
		this_tree_columns = list(expression_matrix_df.columns[:-3])
		heatmap_data_ = expression_matrix_df[expression_matrix_df.index.isin(this_tree_rows)].reindex(this_tree_rows)
		heatmap_data = heatmap_data_.iloc[:,:-3].values.tolist()
		this_tree_rows____ = heatmap_data_.iloc[:,-1].values.tolist()
		self.chart_pt_goenrich_cluster_tojs("pt_goenrich_cluster", this_tree, this_tree_rows____, left_group_colors, this_tree_columns, heatmap_data, height, "pt_goenrich_heatmap.json")
				

	def chart_pt_goenrich_cluster_tojs(self, name, this_tree, this_tree_rows____, left_group_colors,this_tree_columns, heatmap_data, height, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".heatmap_new.js", 'w') as fo:
			a = json.loads(f.read())
			a["tree_1"] = this_tree
			a["rows"] = this_tree_rows____
			a["params"]["left_group_colors"] = left_group_colors
			a["columns"] = this_tree_columns
			a["heatmap_data"] = heatmap_data
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".heatmap_new.js", {"model":"highchart","highchart_type":"tree_heatmap_new","width":950,"height":height}])

	def chart_pt_goenrich_cluster_sub_tojs(self, name, data, title, symbols, colors, legend, categories, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".line.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["params"]['title'] = title
			a["params"]["symbols"] = symbols
			a["params"]["colors"] = colors
			a["legend"] = legend
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".line.js", {"model":"highchart","highchart_type":"showCurve","width":600,"height":400}])

	def chart_pt_keggclass(self, proteinset_file, kegg_stat_file, gene_kegg_level_table_xls):
		# 读入基因集列表
		with open(proteinset_file, 'r') as pset:
			pset_list = [line.split("\t")[0] for line in pset.readlines()]
		stat = pd.read_table(kegg_stat_file, header=0)
		stat.columns=stat.columns.str.replace('Unnamed.*','link')
		level = pd.read_table(gene_kegg_level_table_xls, header=0)
		stat_class = pd.merge(stat, level, on='Pathway_id')

		# 按照kegg官网进行一级分类的排序
		list_custom = ['Metabolism',
					   'Genetic Information Processing',
					   'Environmental Information Processing',
					   'Cellular Processes',
					   'Organismal Systems',
					   'Human Diseases',
					   'Drug Development']
		first_class_index = dict(zip(list_custom,range(len(list_custom))))
		stat_class['first_rank']= stat_class['first_category'].map(first_class_index)
		stat_class.sort_values(['first_rank', 'second_category'], ascending = [True, True], inplace = True)

		stat_class.drop(['graph_id', 'hyperlink', 'graph_png_id', 'first_rank'], axis=1, inplace=True)
		stat_class.rename(columns={'Ko_ids':'ko_ids','Pathway_id':'pathway_id'}, inplace=True)
		stat_class.replace(np.nan, '', regex=True, inplace=True)
		def len_(x):
			while '' in x:
				x.remove('')
			if not x:
				return 0
			return len(x)
		for gene_set in pset_list:
			stat_class.rename(columns={gene_set + '_genes': gene_set + '_geneko'}, inplace=True)
			stat_class[gene_set + '_str'] = stat_class[gene_set + '_geneko'].replace(r'\([^\)]*\)', '', regex=True)
			stat_class[gene_set + '_genes'] = stat_class[gene_set + '_str'].map(lambda x: list(set(x.split(";"))))
			stat_class[gene_set + '_numbers'] = stat_class[gene_set + '_genes'].map(lambda x: len_(x))
			stat_class[gene_set + '_str'] = stat_class[gene_set + '_genes'].map(lambda x: ",".join(x))

		# 导入统计数据
		all_second_cate = list(stat_class['second_category'])
		kegg_class_list = [j for i,j in enumerate(all_second_cate) if all_second_cate.index(j) == i]
		kegg_2to1 = dict(zip(stat_class['second_category'],stat_class['first_category']))

		data_list = list()
		for kegg_class in kegg_class_list:
			data = [
				('first_category', kegg_2to1[kegg_class]),
				('second_category', kegg_class)
			]
			for gene_set in pset_list:
				class_genes = []
				genes_list = list(stat_class[stat_class['second_category'] == kegg_class][gene_set + '_genes'])
				for genes in genes_list:
					class_genes.extend(genes)
				class_genes = list(set(class_genes))
				while '' in class_genes:
					class_genes.remove('')
				data.extend([
					(gene_set + '_genes', class_genes),
					(gene_set + '_genes_num', len(class_genes))
				])
			# data = SON(data)
			data_list.append(data)
		categories = ["".join(map(lambda y:y[0], x.split(' '))) for x in list(set(stat_class['first_category']))]
		#上面代码来自导表
		bar_data = []
		for i in data_list:
		  ic = [ii[1] for ii in i]
		  bar_data.append(ic)
		bar_data_df = pd.DataFrame(bar_data, columns = ['first_category', 'second_category', 'protein_genes', 'protein_genes_num', 'gene_genes', 'gene_genes_num'])
		categories = bar_data_df.iloc[:,1].values.tolist()
		data = [bar_data_df.iloc[:,3].values.tolist(), bar_data_df.iloc[:,5].values.tolist()]
		color_dic = {"Metabolism":"#FF9800", "Genetic Information Processing":"#F44336", "Environmental Information Processing":"#0288D1", "Cellular Processes":"#FF00FF", "Organismal Systems":"#ff0", "Human Diseases":"#388E3C"}
		data_first_category = bar_data_df.iloc[:,0].values.tolist()
		taxon_data = []
		data_colors = []
		for i in sorted(set(data_first_category), key=data_first_category.index):
			tmp_index = [idx for idx, e in enumerate(data_first_category) if e==i]
			taxon_data.append({
				"name":i,
				"color":color_dic[i],
				"start":tmp_index[0],
				"end":tmp_index[-1]
			})
			data_colors.extend([color_dic[i]]*len(tmp_index))
		self.chart_pt_keggclass_tojs("pt_keggclass", data, categories, taxon_data, data_colors, "pt_keggclass.json")	

	def chart_pt_keggclass_tojs(self, name, data, categories, taxon_data, data_colors, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".barline.js", 'w') as fo:
			a = json.loads(f.read())
			a["data_colors"] = data_colors
			a["taxon_data"] = taxon_data
			a["data"] = data
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".barline.js", {"model":"highchart","highchart_type":"showBarLine","width":1000,"height":1060}])

	def chart_pt_keggenrich(self, pt_keggenrich_dir):
		gene_xls_file = pt_keggenrich_dir
		protein_xls_file = pt_keggenrich_dir
		for file in os.listdir(pt_keggenrich_dir):
			if "_protein.list.DE.list.check.kegg_enrichment.xls" in file:
				protein_xls_file = os.path.join(pt_keggenrich_dir,file)
			elif "_gene.list.DE.list.check.kegg_enrichment.xls" in file:
				gene_xls_file = os.path.join(pt_keggenrich_dir,file)
		if gene_xls_file != pt_keggenrich_dir and protein_xls_file != pt_keggenrich_dir:
			protein_xls_df = pd.read_table(protein_xls_file)
			# protein_xls_df = protein_xls_df[(protein_xls_df.iloc[:,2] == 'e')]
			protein_xls_df['enrich_factor'] = protein_xls_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/protein_xls_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
			protein_xls_data_ = protein_xls_df.sort_values(by = [protein_xls_df.columns[6],protein_xls_df.columns[1]], ascending=True).iloc[:20,[3,1,12,6]].values.tolist()
			
			gene_xls_df = pd.read_table(gene_xls_file)
			# gene_xls_df = gene_xls_df[(gene_xls_df.iloc[:,2] == 'e')]
			gene_xls_df['enrich_factor'] = gene_xls_df.iloc[:,4].str.split(pat="/",expand=True).iloc[:,0].astype(int)/gene_xls_df.iloc[:,5].str.split(pat="/",expand=True).iloc[:,0].astype(int)
			gene_xls_data_ = gene_xls_df.iloc[:,[3,12,6]].values.tolist()
			gene_xls_data_dic = {i[0][-5:]:i[-2:] for i in gene_xls_data_}
			data = []
			for i in protein_xls_data_:
				if i[0][-5:] in gene_xls_data_dic.keys():
					data.append([None]+i[1:]+gene_xls_data_dic[i[0][-5:]])
				else:
					data.append([None]+i[1:]+[0,0])
			self.chart_pt_keggenrich_tojs("pt_keggenrich", data, "pt_keggenrich.json")

			data_protein = []
			data_gene = []
			for i in protein_xls_data_:
				for ii in protein_xls_df.iloc[:,[6,12,1,0,3,7]].values.tolist():
					if i[0][-5:] == ii[4][-5:]:
						if float(ii[0]) <= 0.5: 
							data_protein.append({'color':ii[-1], 'x':ii[1], 'name':ii[2], 'size':ii[3]})
				for iii in gene_xls_df.iloc[:,[6,12,1,0,3,7]].values.tolist():
					if i[0][-5:] == iii[4][-5:]:
						if float(iii[0]) <= 0.5: 
							data_gene.append({'color':iii[-1], 'x':iii[1], 'name':iii[2], 'size':iii[3]})
			data = [data_protein, data_gene]
			self.chart_pt_keggenrich_bubble_tojs("pt_keggenrich_bubble", data, "pt_keggenrich_bubble.json")

	def chart_pt_keggenrich_tojs(self, name, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".shadow_bar_groups.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a['params']["linearGradient"] = OrderedDict([("0%","#66CCFF"),("50%","#FFFFFF"), ("100%","#FF66CC")])
			a['params']["linearGradient2"] = OrderedDict([("0%","#FF9966"),("50%","#CCCCCC"), ("100%","#66FFFF")])
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".shadow_bar_groups.js", {"model":"highchart","highchart_type":"shadow_bar_groups","width":850,"height":650}])

	def chart_pt_keggenrich_bubble_tojs(self, name, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".bubble_groups.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".bubble_groups.js", {"model":"highchart","highchart_type":"bubble_groups","width":850,"height":650}])

	def chart_pt_keggenrich_cluster(self, pt_keggenrich_cluster_matrix_file, pt_keggenrich_cluster_tree_file):
		#______子聚类趋势图________
		all_sub = {}
		this_group_dir = os.path.dirname(pt_keggenrich_cluster_matrix_file)
		sub_files_path = []
		for file in os.listdir(this_group_dir):
			if file.startswith('seq.subcluster'):
				sub_files_path.append(os.path.join(this_group_dir,file))
		for sub_file in sub_files_path:
			tmp_n = os.path.basename(sub_file)[4:-4].split('_')
			title = tmp_n[0]+'_'+tmp_n[1]+'('+tmp_n[2]+')'
			this_sub_file_df = pd.read_table(sub_file, index_col=0)
			for i in this_sub_file_df.index:
				all_sub[i] = tmp_n[1]
			this_sub_file_df.loc['mean'] = this_sub_file_df.mean()
			data = this_sub_file_df.values.tolist()

			symbols = [""]*int(tmp_n[2])+['circle']
			colors = ["#F0F0F0"]*int(tmp_n[2])+['#000099']
			legend = this_sub_file_df.index.tolist()
			categories = [i.replace('_protein','').replace('_rna','_gene') for i in this_sub_file_df.columns.tolist()]
			self.chart_pt_keggenrich_cluster_sub_tojs("pt_keggenrich_cluster_"+'_'+tmp_n[0]+'_'+tmp_n[1], data, title, symbols, colors, legend, categories, "pt_proteinsetcluster_SubLine.json")

		pt_keggenrich_cluster_tree_file_str = open(pt_keggenrich_cluster_tree_file, 'r').read()
		this_tree = pt_keggenrich_cluster_tree_file_str.split('\n')[0]
		this_tree_rows = pt_keggenrich_cluster_tree_file_str.split('\n')[1].strip().split(';')
		def randomcolor():
			colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
			color = ""
			for i in range(6):
				color += colorArr[random.randint(0,14)]
			return "#"+color
		color_dic = {}
		for v in set(all_sub.values()):
			color_dic[v] = randomcolor()
		left_group_colors = [color_dic[all_sub[i]] for i in this_tree_rows]

		expression_matrix_df = pd.read_table(pt_keggenrich_cluster_matrix_file, index_col=0)
		height = int(int(expression_matrix_df.shape[0])*8.82+240)
		# this_tree_columns = list(expression_matrix_df.columns[:2])
		# heatmap_data_ = expression_matrix_df[expression_matrix_df.index.isin(this_tree_rows)].reindex(this_tree_rows).iloc[:,[0,1,4]]
		# heatmap_data = heatmap_data_.iloc[:,[0,1]].values.tolist()
		# this_tree_rows____ = heatmap_data_.iloc[:,2].values.tolist()
		this_tree_columns = list(expression_matrix_df.columns[:-3])
		heatmap_data_ = expression_matrix_df[expression_matrix_df.index.isin(this_tree_rows)].reindex(this_tree_rows)
		heatmap_data = heatmap_data_.iloc[:,:-3].values.tolist()
		this_tree_rows____ = heatmap_data_.iloc[:,-1].values.tolist()
		self.chart_pt_keggenrich_cluster_tojs("pt_keggenrich_cluster", this_tree, this_tree_rows____, left_group_colors, this_tree_columns, heatmap_data, height, "pt_goenrich_heatmap.json")
				

	def chart_pt_keggenrich_cluster_tojs(self, name, this_tree, this_tree_rows____, left_group_colors,this_tree_columns, heatmap_data, height, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".heatmap_new.js", 'w') as fo:
			a = json.loads(f.read())
			a["tree_1"] = this_tree
			a["rows"] = this_tree_rows____
			a["params"]["left_group_colors"] = left_group_colors
			a["columns"] = this_tree_columns
			a["heatmap_data"] = heatmap_data
			a["size"]['height'] = 9810
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".heatmap_new.js", {"model":"highchart","highchart_type":"tree_heatmap_new","width":950,"height":height}])

	def chart_pt_keggenrich_cluster_sub_tojs(self, name, data, title, symbols, colors, legend, categories, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".line.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = data
			a["params"]['title'] = title
			a["params"]["symbols"] = symbols
			a["params"]["colors"] = colors
			a["legend"] = legend
			a["categories"] = categories
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".line.js", {"model":"highchart","highchart_type":"showCurve","width":600,"height":400}])

	# preprocess
	def chart_preprocess_assessment(self, assessment_file):
		datalist = [['item'], ['value'], ['category']]
		name = os.path.basename(assessment_file).split('_')[0]
		with open(assessment_file, 'r') as f:
			cat = f.readline().strip().split('\t')
			datalist[0].extend(cat)
			value = f.readline().strip().split('\t')
			datalist[1].extend(value)
			datalist[2].extend(value)
		self.chart_preprocess_assessment_tojs(name, datalist, 'preprocess.assessment.json')

	def chart_preprocess_assessment_tojs(self, name, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + "_assessment.preprocess.js", 'w') as fo:
			a = json.loads(f.read())
			a['dataset'][0]['source'] = data
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
		self.js_list.append([self.output_dir + name + "_assessment.preprocess.js", {}])

	def chart_preprocess_inter_cv(self, inter_cv_file):
		name = 'raw.' if 'raw' in inter_cv_file else ''
		intercv_df = pd.read_table(inter_cv_file, header=0)
		intercv_df = intercv_df.iloc[:-1, :]
		datalist = [{'category': i[0], 'data': sorted(i[4:9]), "name": i[0]} for i in intercv_df.values.tolist()]
		outliers = [{'name': i[0], 'data': i[-2].split(';')} for i in intercv_df.values.tolist()]
		outlier_list = [["name", "x", "y", "value", "category"]]
		for each in outliers:
			tmp_list = [[each['name'], each['name'], float(i), float(i), each['name']] for i in each['data']]
			outlier_list.extend(tmp_list)
		self.chart_preprocess_inter_cv_tojs(name, datalist, outlier_list, 'preprocess.intergroup.cv.json')

	def chart_preprocess_inter_cv_tojs(self, name, data, points, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + "inter-cv.preprocess.js", 'w') as fo:
			a = json.loads(f.read())
			a['dataset'][0]['source'] = data
			a['dataset'][1]['source'] = points
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
		self.js_list.append([self.output_dir + name + "inter-cv.preprocess.js", {}])

	def chart_preprocess_intra_cv(self, intra_cv_file):
		name = '_raw' if 'raw' in intra_cv_file else ''
		intracv_df = pd.read_table(intra_cv_file, header=0)
		intracv_df['groups'] = intracv_df.index
		for each in intracv_df.values.tolist():
			valuelist = ['value'] + each[0:12]
			datalist = each[24:-2]
			group_name = str(each[-1]).split('cv_')[1]
			self.chart_preprocess_intra_cv_tojs(group_name + name, valuelist, datalist, 'preprocess.intragroup.cv.json')

	def chart_preprocess_intra_cv_tojs(self, name, value, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + "_intra-cv.preprocess.js", 'w') as fo:
			a = json.loads(f.read())
			a['dataset'][0]['source'][1] = value
			a['dataset'][1]['source'][0]['data'] = data
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
		self.js_list.append([self.output_dir + name + "_intra-cv.preprocess.js", {}])

	##
	##
	## 下面是生成WGCNA的静态图pdf的代码
	##

	def chart_wgcna_sample_tree(self, sample_tree_table, group_dict=None):
		linkage_pd = pd.read_table(sample_tree_table, header=None)
		linkage_list = [list(x) for x in linkage_pd.as_matrix()]
		ordered_leaf = pd.read_table(sample_tree_table[:-3] + 'order.txt', header=None)[0]
		ordered_seqs = [x.strip() for x in open(sample_tree_table[:-3] + 'labels.txt')]
		newick_tree = convert2(sample_tree_table, ordered_leaf, ordered_seqs)
		source = list()

		sample2group_source = [["name", "group"]]
		if group_dict:
			for g,ss in group_dict.items():
				for s in ss:
					sample2group_source.append([s, g])
		else:
			for s in samples:
				sample2group_source.append([s, s])
		self.chart_tree("wgcna", ".sample_tree", sample_tree=newick_tree, sample2group_source = sample2group_source, json_mode="wgcna.prepare_cluster.tree.json")

	def chart_wgcna_module_tree(self, module_tree_table, group_dict=None):
		linkage_pd = pd.read_table(module_tree_table, header=None)
		linkage_list = [list(x) for x in linkage_pd.as_matrix()]
		ordered_leaf = pd.read_table(module_tree_table[:-3] + 'order.txt', header=None)[0]
		ordered_seqs = [x.strip() for x in open(module_tree_table[:-3] + 'labels.txt')]
		newick_tree = convert2(module_tree_table, ordered_leaf, ordered_seqs)
		source = list()

		module2group_source = [["name", "group"]]
		if group_dict:
			for g,ss in group_dict.items():
				for s in ss:
					module2group_source.append([s, g])
		else:
			modules = ordered_seqs
			for s in modules:
				module2group_source.append([s, s])

		self.chart_tree("wgcna", ".module_tree", newick_tree, sample2group_source = module2group_source, json_mode="wgcna.module_cluster.tree.json")

	def chart_tree(self, name, out, sample_tree, sample2group_source, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + out + ".heat_corr.js", 'w') as fo:
			a = json.loads(f.read())
			a["dataset"][0]["source"] = sample2group_source
			a["dataset"][0]["categories"] = sample_tree
			if str(sample2group_source[1][0]).startswith('ME') and str(sample2group_source[2][0]).startswith('ME'):
				a["series"][0]["visualMap"][0]["visualColorValue"] = [str(self.rcolor_dict.get(c[0][2:]).upper()) for c in sample2group_source[1:]]

			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + out + ".heat_corr.js", {}])


	def chart_wgcna_module_corr(self, module_corr, module_tree):
		with open(module_tree, 'r') as tree_f:
			newick_tree = tree_f.readline().strip()[:-1]
			new_order = re.findall('[(,]([^(]*?):', newick_tree)
		source = list()
		cluster_pd = pd.read_table(module_corr, index_col=0, header=0)
		cluster_pd = cluster_pd.reindex(index = new_order)
		cluster_pd = cluster_pd[new_order]
		ordered_seqs = list(cluster_pd.columns)
		cluster_pd["seq_id"] = cluster_pd.index
		corr_heat = [[""] + ordered_seqs]
		for line_dict in cluster_pd.iterrows():
			corr_heat.append(
				[line_dict[0]] + [line_dict[1][s] for s in ordered_seqs]
			)

		module2group_source = [["name", "group"]]
		for m in ordered_seqs:
			module2group_source.append([m, m])
		self.chart_heat_tree_wgcna("wgcna", ".mofule_corr", corr_heat, newick_tree, newick_tree, module2group_source, module2group_source, new_order, "wgcna.module_relation.heatmap_tree.json")

	def chart_heat_tree_wgcna(self, name, out, corr_source, sample_tree, gene_tree, sample2group_source, gene2group_source, new_order, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + out + ".heat_corr.js", 'w') as fo:
			a = json.loads(f.read())
			a["dataset"][0]["source"] = corr_source
			a["dataset"][1]["source"] = gene2group_source
			if gene_tree:
				a["dataset"][1]["categories"] = gene_tree
			else:
				del a["dataset"][1]["categories"]
			a["dataset"][2]["source"] = sample2group_source
			if sample_tree:
				a["dataset"][2]["categories"] = sample_tree
			else:
				del a["dataset"][2]["categories"]
				a["series"] = a["series"][:2]
				a["dataset"] = a["dataset"][:2]
				a["legend"] = a["legend"][:2]

			if len(gene2group_source) <= 1:
				##无基因聚类
				if sample_tree:
					a["dataset"] = [a["dataset"][0], a["dataset"][2]]
					a["legend"] = [a["legend"][0],  a["legend"][2]]
					a["series"] = [a["series"][0], a['series'][2]]
					a["series"][1]["datasetIndex"] = 1
					a["series"][1]["visualMap"][0]["legendIndex"] = 1
					a["xAxis"][0]["length"] = "95%"
				else:
					a["dataset"] = [a["dataset"][0]]
					a["legend"] = [a["legend"][0]]
					a["series"] = [a["series"][0]]
					a["xAxis"][0]["length"] = "95%"
			a["series"][0]["visualMap"][0]["visualValue"] = ["#009900", "#FFFFFF", "#FF0000"]
			a["series"][1]["visualMap"][0]["visualValue"] = [str(self.rcolor_dict.get(c[2:]).upper()) for c in new_order]
			a["series"][2]["visualMap"][0]["visualValue"] = [str(self.rcolor_dict.get(c[2:]).upper()) for c in new_order]
			# 类型待确定
			a["chart"]["type"] = "heatmap_tree"
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + out + ".heat_corr.js", {}])

	def chart_wgcna_relation_corr(self, relation_corr, relation_corr_pvalue, module_stat):
		cluster_pd = pd.read_table(relation_corr, index_col=0, header=0)
		ordered_seqs = list(cluster_pd.columns)
		visual_color = [str(self.rcolor_dict.get(c[2:]).upper()) for c in cluster_pd.index]

		cluster_pd["module"] = cluster_pd.index
		corr_heat = [[""] + ordered_seqs]
		for line_dict in cluster_pd.iterrows():
			corr_heat.append(
				[line_dict[0]] + list([float("%0.4f" %line_dict[1][x]) for x in ordered_seqs])
			)
		
		cluster_pd_c = pd.read_table(relation_corr_pvalue, index_col=0, header=0)
		ordered_seqs = list(cluster_pd_c.columns)
		cluster_pd_c["module"] = cluster_pd_c.index
		corr_p = [[""] + ordered_seqs]
		for line_dict in cluster_pd_c.iterrows():
			corr_p.append(
				[line_dict[0]] + list([float("%0.4f" %line_dict[1][x]) for x in ordered_seqs])
			)

		group_source = [["name", "group"]]
		size_stat_pd = pd.read_table(module_stat, header=0)
		for mm in cluster_pd.index.tolist():
			for s,m in zip(size_stat_pd["size"], size_stat_pd["module"]):
				if m == mm[2:]:
					group_source.append([s, "ME" + m])

		# print corr_heat
		# print group_source
		# print corr_p
		self.chart_heat("wgcna", ".relation_heat", corr_heat, group_source, corr_p, visual_color, "wgcna.module_trait.heat.json")

	def chart_heat(self, name, out, corr_source, group_source, p_source, color, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + out + ".heat_corr.js", 'w') as fo:
			a = json.loads(f.read())
			a["dataset"][0]["source"] = corr_source
			a["dataset"][1]["source"] = group_source
			a["dataset"][2]["source"] = p_source
			a["series"][1]["visualMap"][0]["visualValue"] = color
			a["series"][1]["visualMap"][0]["visualColorValue"] = color
			a["legend"][1]["color"] = color
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + out + ".heat_corr.js", {}])

	def get_rcolor_dict(self):
		rcolor_dict = dict()
		rcolor_file = os.path.join(self.mode_mode, 'r_color2hex_ME.tsv')
		if os.path.exists(rcolor_file):
			with open(rcolor_file, 'r') as f:
				for color_dict in csv.DictReader(f, delimiter='\t'):
					rcolor_dict[color_dict['color']] = color_dict['hex']
		return rcolor_dict

	def chart_wgcna_module_column(self, module_stat):
		size_stat_pd = pd.read_table(module_stat, header=0)
		
		source = [
			["item"] + list(size_stat_pd["module"]),
			["series"] + list(size_stat_pd["size"]),
			["category"] + list(size_stat_pd["module"])
		]

		visual_color = [str(self.rcolor_dict.get(c).upper()) for c in size_stat_pd["module"]]
		# print source
		# print visual_color
		self.chart_column_color("wgcna", ".module_stat", source, visual_color, "wgcna.module_stat.column.json")

	def chart_column_color(self, name, out, source, color, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + out + ".column.js", 'w') as fo:
			a = json.loads(f.read())
			a["title"]["text"] = a["title"]["text"].format(sample_name = name)
			a["dataset"][0]["source"] = source
			a["series"][0]["visualMap"][0]["visualColorValue"] = color
			a["legend"][0]["color"] = color
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + out + ".column.js", {}])

	def chart_wgcna_relation_ms(self, gene_trait_corr, seq_annot, module_trait_corr, module_trait_corr_pvalue):
		seq_annot_pd = pd.read_table(seq_annot,index_col=0,header=0)
		corr_pd = pd.read_table(gene_trait_corr,index_col=0,header=0)
		corr_pd.index.name = 'seq_id'
		corr_all = corr_pd.join(seq_annot_pd)
		final_corr = corr_all.reset_index()

		corr = pd.read_table(module_trait_corr, index_col=0, header=0)
		module_list = list(corr.index)
		trait_list = list(corr.columns)

		corr_p = pd.read_table(module_trait_corr_pvalue, index_col=0, header=0)
		
		for module, module_df in corr_all.groupby("module"):
			for trait in trait_list:
				source = [["name", "x", "y", "value", "category"]]
				# module_df.to_csv("test.tsv",sep = "\t", header=True)
				for seq_id, gene_dict in module_df.to_dict("index").items():
					source.append([
						seq_id,
						gene_dict["kme"],
						abs(gene_dict[trait]),
						seq_id,
						"all"
					])
				x_title = "Module Membership in {} module".format(module)
				y_title = "Gene significance"
				title = "Module membership vs. gene significance (cor={}, p<{})".format(corr.loc["ME" + module, trait], corr_p.loc["ME" + module, trait])
				self.chart_scatter("wgcna", "relation_ms." + module + "." + trait, source, x_title, y_title, title, "wgcna.module_trait_ms_gs.scatter.json")

		corr_all = corr_pd.abs().join(seq_annot_pd.loc[:, ["module"]])
		mean_gs = corr_all.groupby("module").mean()
		gs_std = corr_all.groupby("module").std()
		gs_all = mean_gs.join(gs_std, lsuffix="_gs", rsuffix="_gs_std")
		gs_all.index = ["ME"+x for x in gs_all.index]

		corr_p = pd.read_table(module_trait_corr_pvalue, index_col=0, header=0)
		corr_all2 = corr.join(corr_p, lsuffix='_corr', rsuffix="_pvalue")
		corr_all2.index.name = "module"

		module_size_pd = seq_annot_pd.groupby("module").count().iloc[:, [0]]
		module_size_pd.columns = ["size"]
		module_size_pd.index = ["ME"+x for x in module_size_pd.index]
		corr_all2 = corr_all2.join(module_size_pd)
		corr_all2 = corr_all2.join(gs_all)
		corr_all2.reset_index(inplace=True)
		corr_all2_dict_list = corr_all2.to_dict("records")
		for trait in trait_list:
			source = [
				["item"] + list(corr_all2["module"]),
				["series"] + list(corr_all2[trait + "_gs"]),
				["category"] + list(corr_all2["module"])
			]

			source2 = [
				["item"] + list(corr_all2["module"]),
				["series"] + [[g, s, s, ""] for g,s in zip(list(corr_all2[trait + "_gs"]), list(corr_all2[trait + "_gs_std"]))],
				["category"] + list(corr_all2["module"])
			]

			# visual_color = [self.rcolor_dict.get(c[2:]).upper() for c in list(corr_all2["module"])]
			visual_color = ['#008000' for i in list(corr_all2["module"])]
			# print source
			print source2
			# print visual_color

			self.chart_column_and_conf("wgcna", ".{}.relation_ms".format(trait), source, source2, visual_color, "wgcna.module_trait_ms.box.json")

	def chart_scatter(self, name, out, source, x_title, y_title, title=None, json_mode=None):
		# 旧版插件
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + out + ".scatter.js", 'w') as fo:
			a = json.loads(f.read())
			a["dataset"][0]["source"] = source
			a["xAxis"][0]["text"] = x_title
			a["xAxis"][0]["title"]["text"] = x_title
			a["yAxis"][0]["text"] = y_title
			a["yAxis"][0]["title"]["text"] = y_title
			a["title"]["text"] = title
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + out + ".scatter.js", {}])


	def chart_column_and_conf(self, name, out, source, source2, color, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + out + ".column_conf.js", 'w') as fo:
			a = json.loads(f.read())
			a["dataset"][0]["source"] = source
			a["dataset"][1]["source"] = source2
			a["series"][0]["visualMap"][0]["visualColorValue"] = color
			a["legend"][0]["color"] = color
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + out + ".column_conf.js", {}])

	def chart_wgcna_prepare_curve(self, wgcna_prepare_curve_file):
		curve_file_df = pd.read_table(wgcna_prepare_curve_file)
		sploe = [-1 if i > 0 else 1 for i in list(curve_file_df['slope'])]
		adapt_data_ = map(lambda x, y: x * y, list(curve_file_df['SFT.R.sq']), sploe)
		adapt_data = [list(i) for i in zip(curve_file_df.iloc[:,0].values.tolist(),adapt_data_)]
		average_data = curve_file_df.iloc[:,[0,4]].values.tolist()
		self.chart_wgcna_prepare_curve_tojs("wgcna_prepare_adapt", adapt_data, "wgcna_prepare_adapt.json")	
		self.chart_wgcna_prepare_curve_tojs("wgcna_prepare_average", average_data, "wgcna_prepare_average.json")		

	def chart_wgcna_prepare_curve_tojs(self, name, data, json_mode):
		json_mode = self.mode_dir + "/" + json_mode
		with open(json_mode, 'r') as f, open(self.output_dir + name + ".wgcnascatter.js", 'w') as fo:
			a = json.loads(f.read())
			a["data"] = [data]
			fo.write("var options = ")
			fo.write(json.dumps(a, indent=4))
			self.js_list.append([self.output_dir + name + ".wgcnascatter.js", {"model":"highchart","highchart_type":"plot_big_scatter","width":450,"height":500}])

	#### 生成cmd脚本
	def generate_html_sh(self, para=True):
		self.command_list = list()
		for js in self.js_list:
			js_file = os.path.abspath(js[0])
			para_dict = js[1]
			html_file = os.path.splitext(js_file)[0] + ".html"
			pdf_file = os.path.splitext(js_file)[0] + ".pdf"
			if "model" in para_dict and para_dict["model"] == "highchart":
				mode_dir = self.mode_mode2  # self.mode_mode2 = chart_dir + "/itraq_tmt/highchart_model"
			else:
				mode_dir = self.mode_mode    #  self.mode_mode = chart_dir + "/model"
			html_mode = os.path.join(mode_dir, "sg_chart_model.html")
			#一些页面需要特殊处理的(如添加js或css)，则需为其指定特制的模板文件
			if js_file.endswith('pep_num.bar.js') or js_file.endswith('pep_len.bar.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_FontDeSize.html")
			if js_file.endswith('pep_error.scatter.js') or js_file.endswith('_scatter.scatter.js') or js_file.endswith('_volcano.volcano.js') or js_file.endswith('wgcnascatter.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_BigScatter.html")
			# if js_file.endswith('level_statistics.go_bar.js'):
			if js_file.endswith('.go_bar.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_GoBar.html")
			if js_file.endswith('_pie.go_pie.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_GoPie.html")
			if js_file.endswith('.barline.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_BarLine.html")
			if js_file.endswith('.heatmap.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_TreeHeatmap.html")
			if js_file.endswith('.shadowbar.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_ShadowBar.html")
			if js_file.endswith('.densebubble.js') or js_file.endswith('.scatterbubble.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_Bubble.html")
			if js_file.endswith('.circ.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_Circ.html")
			if js_file.endswith('.showCurve.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_Curve.html")
			if js_file.endswith('sam_pca.scatter.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_PcaScatter.html")
			if js_file.endswith('.venn.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_Venn.html")
			if js_file.endswith('.venn2.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_Venn2.html")
			if js_file.endswith('.nine.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_scatternine.html")
			if js_file.endswith('.heatmap_new.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_TreeHeatmap_new.html")
			if js_file.endswith('.scatter_new.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_ScatterNew.html")
			if js_file.endswith('.double_nest.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_double_nest.html")
			if js_file.endswith('.shadow_bar_groups.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_double_shadow_bar_groups.html")
			if js_file.endswith('.bubble_groups.js'):
				html_mode = os.path.join(mode_dir, "sg_chart_model_bubble_groups.html")
			html = lxml.html.parse(html_mode)
			root = html.getroot()
			scripts = root.getchildren()
			for script in scripts:
				if script.tag == "body":
					div1 = script.getchildren()[0]
					# 修改图片大小
					size = ""
					if "height" in para_dict:
						size += 'height: {}px;'.format(para_dict["height"])
					if "width" in para_dict:
						size += 'width: {}px;'.format(para_dict["width"])

					div1.attrib['style'] += size

				if script.tag == "script" and "src" in script.attrib:
					# print script.attrib
					if script.attrib["src"] == "./sg_chart_model.js":
						script.set('src', js_file)
						script.text = ""
					elif script.attrib["src"].startswith("./"):
						script.set('src', mode_dir + script.attrib["src"][1:])
						script.text = ""
				else:
					if "model" in para_dict and para_dict["model"] == "highchart":
						if "highchart_type" in para_dict:
							#老版画图插件(即highchart)需要指定画当前类型图的函数名，可在对应页面源代码调取查看(位于jQuery.highcharts_report.后面)
							script.text = script.text.replace("showScatterMarkBig", para_dict["highchart_type"])
						if "type" in para_dict and para_dict["type"] == "network":
							script.text += "\n" + '$("#customize_button").remove();'
			html.write(html_file)

			if "delay" in para_dict:
				delay = para_dict["delay"]
			else:
				delay = 5000

			if "convert" in para_dict:
				cmd = "{}/wkhtmltopdf --no-stop-slow-scripts --enable-local-file-access --page-width {} --page-height {} --zoom {} {} {}".format(
					self.wkhtml_dir,
					para_dict["width"],
					para_dict["height"],
					para_dict["zoom"],
					html_file,
					pdf_file
				)
			elif "use_puppeteer" in para_dict and 'sanger-dev' not in str(Config().SOFTWARE_DIR):
				cmd = "{} {} {} {}".format(self.node_path,self.puppeteer_path, html_file, pdf_file)
			else:
				cmd = "{}/phantomjs {}/sg_chart_phantome.js {} {} {}".format(self.phantomjs_dir, self.mode_mode, html_file, pdf_file, delay)

			self.command_list.append(cmd)

		if para:
			with open(self.work_dir + "para_run.sh", 'w') as fo:
				fo.write("\n".join(self.command_list))
				fo.write("\n")

	def para_to_pdf(self):
		'''
		# run in in tools
		'''
		pass

	def to_pdf(self):
		self.generate_html_sh(para=False)
		for command in self.command_list:
			print "command{}".format(command)
			os.system(command)
			# return_code = subprocess.call(command, shell=True)

	def chart_json_batch(self, chart_json):
		def files_all_exit(file_list):
			if sum([1 for i in file_list if os.path.exists(i)]) == len(file_list):
				return True
			else:
				return False
		if "protein_mw_distribution" in chart_json:
			self.chart_protein_mw(chart_json["protein_mw_distribution"])

		if "protein_seq_cover_distribution" in chart_json:
			self.chart_protein_seq_cover(chart_json["protein_seq_cover_distribution"])
		
		if "protein_infomation" in chart_json:
			self.chart_protein_infomation(chart_json["protein_infomation"])

		if "peptite_number_distribution" in chart_json:
			self.chart_peptite_number_distribution(chart_json["peptite_number_distribution"])

		if "peptite_length_distribution" in chart_json:
			self.chart_peptite_length_distribution(chart_json["peptite_length_distribution"])

		# if "peptide_dmass" in chart_json:
		# 	self.chart_peptide_dmass(chart_json["peptide_dmass"])

		if 'preprocess_assess' in chart_json and 'preprocess_box_origin' in chart_json and 'preprocess_summary_origin' in chart_json:
			self.chart_preprocess_assessment(chart_json["preprocess_assess"])
			self.chart_preprocess_inter_cv(chart_json["preprocess_box_origin"])
			self.chart_preprocess_intra_cv(chart_json["preprocess_summary_origin"])

		if 'preprocess_box_raw' in chart_json and 'preprocess_summary_raw' in chart_json:
			self.chart_preprocess_inter_cv(chart_json["preprocess_box_raw"])
			self.chart_preprocess_intra_cv(chart_json["preprocess_summary_raw"])

		if "exp_venn" in chart_json:
			self.chart_exp_venn(chart_json["exp_venn"])

		if "all_annotation_stat" in chart_json:
			self.chart_all_annotation_stat(chart_json["all_annotation_stat"])

		if "go12level_statistics" in chart_json:
			self.chart_go12level_statistics(chart_json["go12level_statistics"])

		if "go123level_statistics" in chart_json:
			self.chart_go123level_statistics(chart_json["go123level_statistics"])

		if "go1234level_statistics" in chart_json:
			self.chart_go1234level_statistics(chart_json["go1234level_statistics"])

		if "kegg_layer" in chart_json:
			self.chart_kegg_layer(chart_json["kegg_layer"])

		if "kegg_pathway" in chart_json:
			self.chart_kegg_pathway(chart_json["kegg_pathway"])

		if "cog_summary" in chart_json:
			self.chart_cog_summary(chart_json["cog_summary"])

		if "pfam_domain" in chart_json:
			self.chart_pfam_domain(chart_json["pfam_domain"])

		if "subloc_stat" in chart_json:
			self.chart_subloc_stat(chart_json["subloc_stat"])

		if "diff" in chart_json and "num_summary" in chart_json and "allsummary" in chart_json:
			diff_files = [chart_json["diff"].format(group_name=group) for group in chart_json["diff_group"]]
			if files_all_exit(diff_files):
				self.chart_diff(chart_json["diff_group"], chart_json["num_summary"],chart_json["allsummary"], diff_files)

		if "expression_matrix" in chart_json and "seq_tree" in chart_json and "sample_tree" in chart_json:
			proteinsetcluster_group_list = os.listdir(os.path.abspath(os.path.join(os.path.dirname(chart_json["expression_matrix"]), "..")))
			expression_matrix_files = [chart_json["expression_matrix"].format(group_name=group) for group in proteinsetcluster_group_list]
			seq_tree_files = [chart_json["seq_tree"].format(group_name=group) for group in proteinsetcluster_group_list]
			sample_tree_files = [chart_json["sample_tree"].format(group_name=group) for group in proteinsetcluster_group_list]
			if files_all_exit(expression_matrix_files) and files_all_exit(seq_tree_files) and files_all_exit(sample_tree_files):
				self.chart_proteinsetcluster(proteinsetcluster_group_list, expression_matrix_files, seq_tree_files, sample_tree_files)

		if "proteinsetgo_class_all" in chart_json and "proteinsetgo_class_updown" in chart_json:
			proteinsetgo_class_all_files = [chart_json["proteinsetgo_class_all"].format(group_name=group) for group in chart_json["diff_group"]]
			proteinsetgo_class_updown_files = [chart_json["proteinsetgo_class_updown"].format(group_name=group) for group in chart_json["diff_group"]]
			if files_all_exit(proteinsetgo_class_all_files) and files_all_exit(proteinsetgo_class_updown_files):
				self.chart_proteinsetcluster_goclass(chart_json["diff_group"], proteinsetgo_class_all_files, proteinsetgo_class_updown_files)

		if "go_enrich_all" in chart_json and "go_enrich_down" in chart_json and "go_enrich_up" in chart_json:
			go_enrich_all_files = [chart_json["go_enrich_all"].format(group_name=group) for group in chart_json["diff_group"]]
			go_enrich_down_files = [chart_json["go_enrich_down"].format(group_name=group) for group in chart_json["diff_group"]]
			go_enrich_up_files = [chart_json["go_enrich_up"].format(group_name=group) for group in chart_json["diff_group"]]
			if files_all_exit(go_enrich_all_files) and files_all_exit(go_enrich_down_files) and files_all_exit(go_enrich_up_files):
				self.chart_proteinsetcluster_gorich(chart_json["diff_group"], go_enrich_all_files, go_enrich_down_files, go_enrich_up_files)

		if "kegg_class_level" in chart_json and "kegg_class_stat_all" in chart_json and "kegg_class_stat_updown" in chart_json:
			kegg_class_stat_all_files = [chart_json["kegg_class_stat_all"].format(group_name=group) for group in chart_json["diff_group"]]
			kegg_class_stat_updown_files = [chart_json["kegg_class_stat_updown"].format(group_name=group) for group in chart_json["diff_group"]]
			# if files_all_exit(kegg_class_stat_all_files) and files_all_exit(kegg_class_stat_updown_files):
			self.chart_proteinsetcluster_keggclass(chart_json["diff_group"], kegg_class_stat_all_files, kegg_class_stat_updown_files, chart_json["kegg_class_level"])

		if "kegg_enrich_all" in chart_json and "kegg_enrich_down" in chart_json and "kegg_enrich_up" in chart_json:
			kegg_enrich_all_files = [chart_json["kegg_enrich_all"].format(group_name=group) for group in chart_json["diff_group"]]
			kegg_enrich_down_files = [chart_json["kegg_enrich_down"].format(group_name=group) for group in chart_json["diff_group"]]
			kegg_enrich_up_files = [chart_json["kegg_enrich_up"].format(group_name=group) for group in chart_json["diff_group"]]
			# if files_all_exit(kegg_enrich_all_files) and files_all_exit(kegg_enrich_down_files) and files_all_exit(kegg_enrich_up_files):
			self.chart_proteinsetcluster_keggrich(chart_json["diff_group"], kegg_enrich_all_files, kegg_enrich_down_files, kegg_enrich_up_files)

		if "geneset_circ_choose" in chart_json and "geneset_circ_zscore" in chart_json:
			geneset_circ_choose_files = [chart_json["geneset_circ_choose"].format(group_name=group) for group in chart_json["diff_group"]]
			geneset_circ_zscore_files = [chart_json["geneset_circ_zscore"].format(group_name=group) for group in chart_json["diff_group"]]
			self.chart_proteinset_circ(chart_json["diff_group"], geneset_circ_choose_files, geneset_circ_zscore_files)

		if "ppi_centrality_all" in chart_json and "ppi_centrality_up" in chart_json and "ppi_centrality_down" in chart_json:
			ppi_centrality_all_files = [chart_json["ppi_centrality_all"].format(group_name=group) for group in chart_json["diff_group"]]
			ppi_centrality_up_files = [chart_json["ppi_centrality_up"].format(group_name=group) for group in chart_json["diff_group"]]
			ppi_centrality_down_files = [chart_json["ppi_centrality_down"].format(group_name=group) for group in chart_json["diff_group"]]
			if files_all_exit(ppi_centrality_all_files) and files_all_exit(ppi_centrality_up_files) and files_all_exit(ppi_centrality_down_files):
				self.chart_ppi_centrality(chart_json["diff_group"], ppi_centrality_all_files, ppi_centrality_up_files, ppi_centrality_down_files)

		if "ppi_degree_all" in chart_json and "ppi_degree_up" in chart_json and "ppi_degree_down" in chart_json:
			ppi_degree_all_files = [chart_json["ppi_degree_all"].format(group_name=group) for group in chart_json["diff_group"]]
			ppi_degree_up_files = [chart_json["ppi_degree_up"].format(group_name=group) for group in chart_json["diff_group"]]
			ppi_degree_down_files = [chart_json["ppi_degree_down"].format(group_name=group) for group in chart_json["diff_group"]]
			if files_all_exit(ppi_degree_all_files) and files_all_exit(ppi_degree_up_files) and files_all_exit(ppi_degree_down_files):
				self.chart_ppi_degree(chart_json["diff_group"], ppi_degree_all_files, ppi_degree_up_files, ppi_degree_down_files)

		if "cog_all" in chart_json and "cog_up_down" in chart_json:
			cog_all_files = [chart_json["cog_all"].format(group_name=group) for group in chart_json["diff_group"]]
			cog_up_down_files = [chart_json["cog_up_down"].format(group_name=group) for group in chart_json["diff_group"]]
			if files_all_exit(cog_all_files) and files_all_exit(cog_up_down_files):
				self.chart_cog(chart_json["diff_group"], cog_all_files, cog_up_down_files)

		if "pfam_all" in chart_json and "pfam_up_down" in chart_json:
			pfam_all_files = [chart_json["pfam_all"].format(group_name=group) for group in chart_json["diff_group"]]
			pfam_up_down_files = [chart_json["pfam_up_down"].format(group_name=group) for group in chart_json["diff_group"]]
			# if files_all_exit(pfam_all_files) and files_all_exit(pfam_up_down_files):
			self.chart_pfam(chart_json["diff_group"], pfam_all_files, pfam_up_down_files)

		if "subloc_all" in chart_json and "subloc_up_down" in chart_json:
			subloc_all_files = [chart_json["subloc_all"].format(group_name=group) for group in chart_json["diff_group"]]
			subloc_up_down_files = [chart_json["subloc_up_down"].format(group_name=group) for group in chart_json["diff_group"]]
			if files_all_exit(subloc_all_files) and files_all_exit(subloc_up_down_files):
				self.chart_subloc(chart_json["diff_group"], subloc_all_files, subloc_up_down_files)

		if "sample_corr_matrix" in chart_json and "sample_corr_tree" in chart_json:
			sample_corr_matrix_file = chart_json["sample_corr_matrix"]
			sample_corr_tree_file = chart_json["sample_corr_tree"]
			if 'group_table' in chart_json:
				self.chart_exp_corr(sample_corr_matrix_file, sample_corr_tree_file, chart_json['protein_group'])
			else:
				self.chart_exp_corr(sample_corr_matrix_file, sample_corr_tree_file)

		if "sample_pca_sites" in chart_json and "sample_pca_propor" in chart_json and 'group_file' in chart_json:
			sample_pca_sites_file = chart_json["sample_pca_sites"]
			sample_pca_propor_file = chart_json["sample_pca_propor"]
			self.chart_sample_pca(sample_pca_sites_file, sample_pca_propor_file, chart_json['group_file'])

		if "diff_path_for_proteinset_for_venn" in chart_json:
			diff_path_for_proteinset = chart_json["diff_path_for_proteinset_for_venn"]
			self.chart_proteinset_venn(diff_path_for_proteinset)

		if "protein_transcript_overview" in chart_json:
			protein_transcript_overview_file = chart_json["protein_transcript_overview"]
			self.chart_protein_transcript_overview(protein_transcript_overview_file)

		if "protein_transcript_baseinfo" in chart_json:
			protein_transcript_baseinfo_file = chart_json["protein_transcript_baseinfo"]
			self.chart_protein_transcript_baseinfo(protein_transcript_baseinfo_file)

		if "protein_transcript_diff_allvenn" in chart_json:
			protein_transcript_diff_allvenn_file = chart_json["protein_transcript_diff_allvenn"]
			self.chart_protein_transcript_diff_allvenn(protein_transcript_diff_allvenn_file)

		if "protein_transcript_diff_relatevenn" in chart_json:
			protein_transcript_diff_relatevenn_file = chart_json["protein_transcript_diff_relatevenn"]
			self.chart_protein_transcript_diff_relatevenn(protein_transcript_diff_relatevenn_file)

		if "protein_transcript_diff_relatenine" in chart_json and "protein_transcript_diff_relatenine_relateparams" in chart_json:
			protein_transcript_diff_relatenine_file = chart_json["protein_transcript_diff_relatenine"]
			relatenine_params_from_mongo_maintable_params = chart_json["protein_transcript_diff_relatenine_relateparams"]
			self.chart_protein_transcript_diff_relatenine(protein_transcript_diff_relatenine_file,relatenine_params_from_mongo_maintable_params)
		
		if "pt_expression_matrix" in chart_json and "pt_seq_tree" in chart_json:
			expression_matrix_files = chart_json["pt_expression_matrix"]
			seq_tree_files = chart_json["pt_seq_tree"]
			self.chart_pt_proteinsetcluster(expression_matrix_files, seq_tree_files)

		if "pt_corrScatter" in chart_json and "pt_corrScatter_pvalue" in chart_json:
			pt_corrScatter_file = chart_json["pt_corrScatter"]
			pt_corrScatter_pvalue_file = chart_json["pt_corrScatter_pvalue"]
			self.chart_pt_proteinsetcluster_corr(pt_corrScatter_file, pt_corrScatter_pvalue_file)

		if "pt_goclass" in chart_json:
			pt_goclass_file = chart_json["pt_goclass"]
			self.chart_pt_goclass(pt_goclass_file)

		if "pt_goenrich" in chart_json:
			pt_goenrich_dir = chart_json["pt_goenrich"]
			self.chart_pt_goenrich(pt_goenrich_dir)

		if "pt_goenrich_cluster_matrix" in chart_json and "pt_goenrich_cluster_tree" in chart_json:
			pt_goenrich_cluster_matrix_file = chart_json["pt_goenrich_cluster_matrix"]
			pt_goenrich_cluster_tree_file = chart_json["pt_goenrich_cluster_tree"]
			self.chart_pt_goenrich_cluster(pt_goenrich_cluster_matrix_file, pt_goenrich_cluster_tree_file)

		if "pt_keggclass_level" in chart_json and "kegg_class_stat_all" in chart_json and "kegg_class_stat_updown" in chart_json:
			kegg_class_stat_all_files = [chart_json["kegg_class_stat_all"].format(group_name=group) for group in chart_json["diff_group"]]
			kegg_class_stat_updown_files = [chart_json["kegg_class_stat_updown"].format(group_name=group) for group in chart_json["diff_group"]]
			if files_all_exit(kegg_class_stat_all_files) and files_all_exit(kegg_class_stat_updown_files):
				self.chart_proteinsetcluster_keggclass(chart_json["diff_group"], kegg_class_stat_all_files, kegg_class_stat_updown_files, chart_json["pt_keggclass_level"])

		if "pt_proteinset_file" in chart_json and "pt_kegg_stat_file" in chart_json and "pt_gene_kegg_level_table_xls" in chart_json:
			self.chart_pt_keggclass(chart_json["pt_proteinset_file"], chart_json["pt_kegg_stat_file"], chart_json["pt_gene_kegg_level_table_xls"])

		if "pt_keggenrich" in chart_json:
			pt_keggenrich_dir = chart_json["pt_keggenrich"]
			self.chart_pt_keggenrich(pt_keggenrich_dir)

		if "pt_keggenrich_cluster_matrix" in chart_json and "pt_keggenrich_cluster_tree" in chart_json:
			pt_keggenrich_cluster_matrix_file = chart_json["pt_keggenrich_cluster_matrix"]
			pt_keggenrich_cluster_tree_file = chart_json["pt_keggenrich_cluster_tree"]
			self.chart_pt_keggenrich_cluster(pt_keggenrich_cluster_matrix_file, pt_keggenrich_cluster_tree_file)

		if "sample_tree" in chart_json and "group_dict" in chart_json:
			self.chart_wgcna_sample_tree(chart_json["sample_tree"],chart_json["group_dict"])

		if "module_tree" in chart_json:
			self.chart_wgcna_module_tree(chart_json["module_tree"])

		if "module_corr" in chart_json and "module_corr_tree" in chart_json:
			self.chart_wgcna_module_corr(chart_json["module_corr"],chart_json["module_corr_tree"])

		if "module_stat" in chart_json:
			self.chart_wgcna_module_column(chart_json["module_stat"])

		if "relation_corr" in chart_json and "relation_corr_pvalue" in chart_json and "module_stat" in chart_json:
			self.chart_wgcna_relation_corr(chart_json["relation_corr"],chart_json["relation_corr_pvalue"],chart_json["module_stat"])

		if "gene_trait_corr" in chart_json and "seq_annot" in chart_json and "relation_corr" in chart_json and "relation_corr_pvalue" in chart_json:
			self.chart_wgcna_relation_ms(chart_json["gene_trait_corr"],chart_json["seq_annot"],chart_json["relation_corr"],chart_json["relation_corr_pvalue"])

		if "wgcna_prepare_curve" in chart_json:
			self.chart_wgcna_prepare_curve(chart_json["wgcna_prepare_curve"])

		self.generate_html_sh()  #这个是在tool的用法中，仅生成para_run.sh文件，由tool并行启动
		# self.to_pdf()    #这个是在package的用法中，测试本package用的
