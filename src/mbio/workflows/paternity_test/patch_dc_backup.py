# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
# last modified:201703

'''医学检验所-无创产前亲子鉴定流程'''
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import re
import datetime
from bson.objectid import ObjectId
import json
import shutil
import gevent


class PatchDcBackupWorkflow(Workflow):
	def __init__(self, wsheet_object):
		'''
        :param wsheet_object:
        '''
		self._sheet = wsheet_object
		super(PatchDcBackupWorkflow, self).__init__(wsheet_object)
		options = [
			{"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"},  # fastq所在路径(文件夹
			{"name": "cpu_number", "type": "int", "default": 4},  # cpu个数
			{"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 参考序列
			{"name": "targets_bedfile", "type": "infile", "format": "paternity_test.rda"},

			{"name": "err_min", "type": "int", "default": 2},  # 允许错配数
			{"name": "ref_point", "type": "infile", "format": "paternity_test.rda"},  # 参考位点
			{"name": "dedup_num", "type": "int", "default": 2},  # 查重样本数
			{"name": "batch_id", "type": "string"},
			{"name": "update_info", "type": "string"},
			{"name": "member_id", "type": "string"},
			{"name":"direct_get_path", "type":"string"}

		]
		self.add_option(options)
		self.tools = []
		self.tool = []
		self.tools_analysis = []
		self.tools_result = []
		self.tools_dedup_f = []
		self.rdata = []
		self.set_options(self._sheet.options())
		self.step.add_steps("pt_analysis", "result_info", "retab",
		                    "de_dup1", "de_dup2")
		# self.update_status_api = self.api.pt_update_status

	def check_options(self):
		'''
        检查参数设置
        '''
		if not self.option('fastq_path').is_set:
			raise OptionError('必须提供fastq文件所在的路径')
		return True

	def set_step(self, event):
		if 'start' in event['data'].keys():
			event['data']['start'].start()
		if 'end' in event['data'].keys():
			event['data']['end'].finish()
		self.step.update()

	def finish_update(self, event):
		step = getattr(self.step, event['data'])
		step.finish()
		self.step.update()

	def fastq2mongo_run(self):
		api_read_tab = self.api.tab_file
		fastq = os.listdir(self.option('fastq_path').prop['path'])
		file = []
		file_type = []
		for j in fastq:
			m = re.match('(.*)_R1.fastq.gz', j)
			if m:
				type = api_read_tab.type(m.group(1))
				# if type == "pt" or type == "ppt":
				# 	file.append(m.group(1))
				# else:
				# 	pass
				file.append(m.group(1))
				file_type.append(type)
		n = 0
		for i in file:
			x = api_read_tab.tab_exist(i)
			if x:
				if self.option('direct_get_path') == 'True':
					self.logger.info('{}样本已存在于数据库'.format(i))
				elif self.option('direct_get_path') == 'False':
					raise Exception('请确认{}样本是否重名'.format(i))
			else:
				fastq2mongo = self.add_module("paternity_test.fastq2mongo_dc")
				self.step.add_steps('fastq2mongo{}'.format(n))
				fastq2mongo.set_options({
					"sample_id": i,
					"fastq_path": self.option("fastq_path"),
					"cpu_number": self.option("cpu_number"),
					"ref_fasta": self.option("ref_fasta"),
					"targets_bedfile": self.option("targets_bedfile"),
					"batch_id": self.option('batch_id'),
					"type": file_type[file.index(i)]
					# "type": 'pt'
				}
				)
				step = getattr(self.step, 'fastq2mongo{}'.format(n))
				step.start()
				fastq2mongo.on('end', self.finish_update, 'fastq2mongo{}'.format(n))
				self.tools.append(fastq2mongo)
				n += 1

		for j in range(len(self.tools)):
			self.tools[j].on('end', self.set_output, 'fastq2mongo')

		if self.tools:
			if len(self.tools) > 1:
				self.on_rely(self.tools, self.pt_analysis_run)
			elif len(self.tools) == 1:
				self.tools[0].on('end', self.pt_analysis_run)
				# self.result_info.on('end', self.dedup_run)
		else:
			# self.result_info.on('end', self.dedup_run)
			self.pt_analysis_run()

		for tool in self.tools:
			tool.run()

	# gevent.spawn(self.run_tools, self.tools)
	#
	# def run_tools(self, tools):
	#     for i in tools:
	#         gevent.sleep(1)
	#         self.logger.info(str(i)+ ':' + str(datetime.datetime.now()))
	#         i.run()

	def pt_analysis_run(self):
		api_read_tab = self.api.tab_file
		self.family_id = api_read_tab.family_unanalysised()  # tuple
		if not self.family_id:
			raise Exception("没有符合条件的家系")

		api_read_tab = self.api.tab_file
		n = 0
		self.num_list = []
		self.dedup_list = []
		self.father = []
		self.preg =[]
		self.mother = []
		for p in range(len(self.family_id)):
			# temp = re.match('WQ([1-9].*)-F.*', self.family_id[p][0])
			# num = int(temp.group(1))
			# self.num_list = range(num - self.option('dedup_num'), num + self.option('dedup_num') + 1)
			name_list = []
			api_read_tab.export_tab_file(self.family_id[p][1], self.output_dir)
			api_read_tab.export_tab_file(self.family_id[p][2], self.output_dir)


			x = api_read_tab.dedup_sample()
			if len(x):  # 如果库中能取到前后的样本
				for k in range(len(x)):
					api_read_tab.export_tab_file(x[k], self.output_dir)
					if x[k] != self.family_id[p][0] and x[k] != self.family_id[p][0] + '1':
						name_list.append(x[k])
			if name_list == []:
				pass
			else:
				self.father.append(self.family_id[p][0])
				self.mother.append(self.family_id[p][1])
				self.preg.append(self.family_id[p][2])
				self.dedup_list.append(name_list)
				self.tool.append([])

		for i in range(len(self.family_id)):
			dad_id = self.family_id[i][0]
			mom_id = self.family_id[i][1]
			preg_id = self.family_id[i][2]

			pt_analysis = self.add_module("paternity_test.pt_analysis")
			self.step.add_steps('pt_analysis{}'.format(n))
			# dad_tab = api_read_tab.export_tab_file(dad_id, self.output_dir)
			# mom_tab = api_read_tab.export_tab_file(mom_id, self.output_dir)
			# preg_tab =  api_read_tab.export_tab_file(preg_id, self.output_dir)

			pt_analysis.set_options({
				"dad_tab": self.output_dir + '/' + dad_id + '.tab',  # 数据库的tab文件
				"mom_tab": self.output_dir + '/' + mom_id + '.tab',
				"preg_tab": self.output_dir + '/' + preg_id + '.tab',
				"ref_point": self.option("ref_point"),
				"err_min": self.option("err_min")
			})
			# self.rdata = self.work_dir + '/PtAnalysis/FamilyMerge/output/family_joined_tab.Rdata'
			step = getattr(self.step, 'pt_analysis{}'.format(n))
			step.start()
			pt_analysis.on('end', self.finish_update, 'pt_analysis{}'.format(n))
			self.tools_analysis.append(pt_analysis)
			n += 1
		for j in range(len(self.tools_analysis)):
			self.tools_analysis[j].on('end', self.set_output, 'pt_analysis')
			self.tools_analysis[j].on('end', self.result_info_run, j)

		# if len(self.tools_analysis) > 1:
		# 	self.on_rely(self.tools_analysis, self.result_info_run,len(self.tools_analysis))
		# elif len(self.tools_analysis) == 1:
		# 	self.tools_analysis[0].on('end', self.result_info_run,len(self.tools_analysis))

		for t in self.tools_analysis:
			t.run()

	def result_info_run(self, event):
		n = event['data']
		result_info = self.add_tool("paternity_test.result_info")
		self.step.add_steps('result_info{}'.format(n))
		if int(n) == 0:
			results = os.listdir(self.work_dir + "/PtAnalysis/FamilyMerge/output")
			for f in results:
				if re.match(r'.*family_joined_tab\.Rdata$', f):
					rdata = f
				else:
					print "Oops!"
			self.rdata = self.work_dir + '/PtAnalysis/FamilyMerge/output/' + rdata
			self.father_sample = rdata.split('_')[0]
			self.mom_sample = rdata.split('_')[1]
			self.preg_sample = rdata.split('_')[2]
		else:
			results = os.listdir(self.work_dir + "/PtAnalysis{}/FamilyMerge/output".format(n))
			for f in results:
				if re.match(r'.*family_joined_tab\.Rdata$', f):
					rdata = f
				elif re.match(r'.*family_joined_tab\.txt$', f):
					pass
			self.rdata = self.work_dir + '/PtAnalysis{}/FamilyMerge/output/'.format(n) + rdata
			self.father_sample = rdata.split('_')[0]
			self.mom_sample = rdata.split('_')[1]
			self.preg_sample = rdata.split('_')[2]

		# self.logger.info(self.father_sample)
		# self.logger.info(self.father)
		# self.logger.info(self.mom_sample)
		# self.logger.info(self.mother)
		# self.logger.info(self.preg_sample)
		# self.logger.info(self.preg)

		if self.father_sample in self.father:
			q = self.father.index(self.father_sample)
			result_info.on('end', self.dedup_run, q)
		elif self.mom_sample in self.mother:
			q = self.mother.index(self.mom_sample)
			result_info.on('end', self.dedup_run, q)
		elif self.preg_sample in self.preg:
			q = self.preg.index(self.preg_sample)
			result_info.on('end', self.dedup_run, q)
		else:
			result_info.on('end', self.end)

		result_info.set_options({
			"tab_merged": self.rdata
		})
		step = getattr(self.step, 'result_info{}'.format(n))
		step.start()
		result_info.on('end', self.finish_update, 'result_info{}'.format(n))
		result_info.on('end', self.set_output, 'result_info')
		result_info.run()

	def dedup_run(self, event):
		n = 0
		q = event['data']
		name_list = self.dedup_list[q]
		self.tools_dedup = []
		self.logger.info(q)
		self.logger.info(name_list)
		# family = self.father[q].split('-')[0]
		for i in name_list:
			pt_analysis_dedup = self.add_module("paternity_test.pt_analysis")
			self.step.add_steps('dedup_{}'.format(n))
			pt_analysis_dedup.set_options({
				"dad_tab": self.output_dir + '/'+ i + '.tab',  # 数据库的tab文件
				"mom_tab": self.output_dir + '/' + self.mother[q] +'.tab',
				"preg_tab": self.output_dir +'/' + self.preg[q]+'.tab',
				"ref_point": self.option("ref_point"),
				"err_min": self.option("err_min")
			}
			)
			step = getattr(self.step, 'dedup_{}'.format(n))
			step.start()
			pt_analysis_dedup.on('end', self.finish_update, 'dedup_{}'.format(n))
			self.tools_dedup.append(pt_analysis_dedup)
			n += 1
		self.tool[q] = self.tools_dedup

		for j in range(len(self.tool[q])):
			self.tool[q][j].on('end', self.set_output, 'dedup')

		x = len(self.tool) - 1

		if self.list_2D(self.tool):
			if len(self.tool[x]) > 1:
				self.on_rely(self.tool[x], self.end)
			else:
				self.tool[x][0].on("end", self.end)

			for tool in self.tool:
				for tool_i in tool:
					tool_i.run()

	def list_2D(self, name):
		for m in name:
			if m == []:
				temp = False
				break
			else:
				temp = True
		return temp

	def linkdir(self, dirpath, dirname):
		"""
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
		allfiles = os.listdir(dirpath)
		newdir = os.path.join(self.output_dir, dirname)
		if not os.path.exists(newdir):
			os.mkdir(newdir)
		oldfiles = [os.path.join(dirpath, i) for i in allfiles]
		newfiles = [os.path.join(newdir, i) for i in allfiles]
		for newfile in newfiles:
			if os.path.exists(newfile):
				if os.path.isfile(newfile):
					os.remove(newfile)
				else:
					os.system('rm -r %s' % newfile)
		for i in range(len(allfiles)):
			if os.path.isfile(oldfiles[i]):
				os.link(oldfiles[i], newfiles[i])
			elif os.path.isdir(oldfiles[i]):
				os.system('cp -r %s %s' % (oldfiles[i], newdir))

	def set_output(self, event):
		obj = event["bind_object"]
		if event['data'] == "fastq2mongo":
			self.linkdir(obj.output_dir + '/fastq2tab', self.output_dir)

		if event['data'] == "pt_analysis":
			self.linkdir(obj.output_dir + '/family_analysis', self.output_dir)
			self.linkdir(obj.output_dir + '/family_merge', self.output_dir)
			# api_main = self.api.sg_paternity_test
			# self.flow_id = api_main.add_pt_task_main(err_min=self.option("err_min"), task=None, flow_id=None)

		if event['data'] == "result_info":
			self.linkdir(obj.output_dir, self.output_dir)
			# api_main = self.api.sg_paternity_test
			# api_main.add_pt_figure(obj.output_dir)

		if event['data'] == "dedup":
			self.linkdir(obj.output_dir + '/family_analysis', self.output_dir)

	def run(self):
		self.fastq2mongo_run()
		super(PatchDcBackupWorkflow, self).run()

	def end(self):
		api_main = self.api.sg_paternity_test
		api_read_tab = self.api.tab_file
		api_update_status = self.api.pt_customer

		results = os.listdir(self.output_dir)

		for i in range(len(self.family_id)):
			dad_id = self.family_id[i][0]
			mom_id = self.family_id[i][1]
			preg_id = self.family_id[i][2]

			api_read_tab.update_pt_tab(dad_id)
			self.father_id = api_main.add_sg_father(dad_id, mom_id, preg_id, self.option('batch_id'),
			                                        self.option("member_id"))
			api_main.add_sg_ref_file(self.father_id, self.option('ref_fasta').prop['path'],
			                         self.option('targets_bedfile').prop['path'],
			                         self.option('ref_point').prop['path'], self.option('fastq_path').prop['path'])

			self.pt_father_id = api_main.add_pt_father(father_id=self.father_id, err_min=self.option("err_min"),
			                                           dedup=self.option('dedup_num'))

			dedup = '.*' + mom_id + '_' + preg_id + '_family_analysis.txt'
			dedup1 = '.*_NA_' + preg_id + '_family_analysis.txt'
			dedup2 = '.*' + mom_id + '_NA_family_analysis.txt'
			for f in results:
				if re.search(dedup, f):
					api_main.add_analysis_tab(self.output_dir + '/' + f, self.pt_father_id)
				elif re.search(dedup1, f):
					api_main.add_analysis_tab(self.output_dir + '/' + f, self.pt_father_id)
				elif re.search(dedup2, f):
					api_main.add_analysis_tab(self.output_dir + '/' + f, self.pt_father_id)
				elif f == dad_id + '_' + mom_id + '_' + preg_id + '_family_joined_tab.txt':
					api_main.add_sg_pt_father_detail(self.output_dir + '/' + f, self.pt_father_id)
				elif f == mom_id + '_' + preg_id + '_info_show.txt':
					api_main.add_info_detail(self.output_dir + '/' + f, self.pt_father_id)
				elif f == dad_id + '_' + mom_id + '_' + preg_id + '_test_pos.txt':
					api_main.add_test_pos(self.output_dir + '/' + f, self.pt_father_id)
				elif f == dad_id + '_' + mom_id + '_' + preg_id + '_family.png':
					file_dir = self.output_dir + '/' + dad_id + '_' + mom_id + '_' + preg_id
					api_main.add_pt_father_figure(file_dir, self.pt_father_id)

			#如遇深度较低的样本，在结果处报错
			if dad_id + '_' + mom_id + '_' + preg_id + '_family.png' not in results:
				api_main.has_problem(self.pt_father_id, dad_id)
			if mom_id + '_' + preg_id + '_info_show.txt' not in results:
				api_main.update_infoshow(self.pt_father_id,mom_id,preg_id)


			# 把筛选的内容提取到主表中去
			api_main.add_father_result(self.father_id, self.pt_father_id,dad_id)
			api_main.add_father_qc(self.father_id, self.pt_father_id)
			# 更新单次运行的状态
			api_main.update_sg_pt_father(self.pt_father_id)

		super(PatchDcBackupWorkflow, self).end()
