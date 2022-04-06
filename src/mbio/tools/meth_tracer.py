# -*- coding: utf-8 -*-
"""
Created on Fri Sep 08 15:27:18 2017

@author: Lotus
last modified 20171209
"""
#import urllib
from pymongo import MongoClient
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
import random
import time
from bson.objectid import ObjectId

from email import encoders
from email.header import Header
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.utils import parseaddr, formataddr
import smtplib
import StringIO

_jobID=None
class MethTracerAgent(Agent):
	def __init__(self, parent):
		super(MethTracerAgent,self).__init__(parent)
		options = [
			{"name":"path","type":"string","default":None},
			{"name":"model","type":"string","default":"1"},
			#{"name":"email","type":"string","default":None},
			{"name":"jobID","type":"string","default":None},
			]
		self.add_option(options)

		_jobID = self.option('jobID')
		self.logger.info(_jobID)
		
	def getPath(self):
		return self.option("path")
	
	def getModel(self):
		return self.option("model")
	
	def getEmail(self):
		return self.option("email")
	
	def check_options(self):
		if not self.option('path'):
			raise OptionError("parameter path must be provided!")
		else:
			_path = self.option('path')
			self.logger.info("path:%s" % self.option("path"))
		if not self.option('model'):
			raise OptionError("parameter model must be provided!")
		else:
			_model = self.option("model")
			self.logger.info("model:%s" % self.option("model"))
	def set_resource(self):
		if(self.option('model') == '1'):
			self._cpu = 8
			self._memory = '15G'
		elif(self.option('model')== '5'):
			self._cpu = 1
			self._memory = '6G'
		else:
			self._cpu = 1
			self._memory = '6G'
		
	def end(self):
		# modified by Luo Tong, 20171209
		def send_res_email(job_id, to_addr = None, client = MongoClient("10.100.200.131", 27017), data = None):
			if not isinstance(job_id, ObjectId):
				job_id = ObjectId(job_id)
			jobInfo = client.cmeth.job.find_one({'_id':job_id})
			if not to_addr:
				to_addr = jobInfo['email']
			def _format_addr(s):
				name, addr = parseaddr(s)
				return formataddr(( \
					Header(name, 'utf-8').encode(), \
					addr.encode('utf-8') if isinstance(addr, unicode) else addr))
			
			def get_data(job_id, client = client):
				job_score = client.cmeth.job_detail.find_one({'job_id':ObjectId(job_id)})
				job_pos = client.cmeth.job_bar.find_one({'job_id':ObjectId(job_id)})
				result = {}
				result['score'] = job_score['score']
				result['tumor_burden'] = job_score['tumor_burden']
				result['predicted_type'] = job_score['predicted_type']
				result['position'] = job_pos['position']
				result['proportion'] = job_pos['proportion']
				return result
			if not data:
				data = get_data(job_id, client)

			if to_addr and jobInfo['has_sent_email']==1:
				from_addr = 'ctmethtracer@sina.com'
				password = 'MethTracer2017'
				smtp_server = 'smtp.sina.com'
				
				## email content
				content = ''
				content += 'Dear ctMethTracer user,\r\n\r\n'
				content += 'This email is to inform you that your job {} is finished, and the txt file attached is your result. For an interactive view, you can inquire your result at  http://ctmeth.i-sanger.com.\r\n'.format(str(job_id))
				subject = 'Your job {} at ctMethTracer finished'.format(str(job_id))
				
				## email attachment
				attachment = ''
				attachment += 'predicted_type\t{}\r\n'.format(data['score'][data['predicted_type']]['name'])
				attachment += 'score\t{}\r\n'.format(data['score'][data['predicted_type']]['score'])
				attachment += '\r\n========================================================\r\n\r\n'
				attachment += 'CpG position\r\n'
				for pos in data['position']:
					attachment += '{}\t{}\r\n'.format(pos, data['position'][pos])
				attachment += '\r\n========================================================\r\n\r\n'
				attachment += 'score\r\n'
				for cancer in data['score']:
					attachment += '{}\t{}\r\n'.format(data['score'][cancer]['name'], data['score'][cancer]['score'])
				## generate attachment file 'attach_txt'
				sio = StringIO.StringIO()
				sio.write(attachment)
				attach_txt = sio.getvalue()
				msg = MIMEMultipart()
				msg['From'] = _format_addr(from_addr)
				msg['To'] = _format_addr(to_addr)
				msg['Subject'] = Header(subject, 'utf-8')
				msg.attach(MIMEText(content, 'plain', 'utf-8'))
				mime = MIMEText(attach_txt, 'base64', 'utf-8')
				mime['Content-Type'] = 'application/octet-stream'
				mime['Content-Disposition'] = 'attachment; filename="ctMethTracer_{}.txt"'.format(str(job_id))
				msg.attach(mime)
				
				## send email
				server = smtplib.SMTP(smtp_server, 25)
				server.set_debuglevel(1)
				server.login(from_addr, password)
				try:
					server.sendmail(from_addr, [to_addr], msg.as_string())
					server.quit()
					client.cmeth.job.update_one({'_id':ObjectId(job_id)}, {'$set':{'has_sent_email':2}})
				except Exception as note:
					client.cmeth.job.update_one({'_id':ObjectId(job_id)}, {'$set':{'note':note}})
		try:
			send_res_email(job_id=self.option('jobID'))
		except Exception, e:
			self.logger.info(e)
		super(MethTracerAgent,self).end()

class MethTracerTool(Tool):
	def __init__(self, config):
		super(MethTracerTool, self).__init__(config)
		self._cpu = 1
		if(self.option('model') == '1'):
			self._cpu = 8
			#self._memory = '15G'
		else:
			self._cpu = 1
			#self._memory = '5G'
			
		self.Python_path = '/miniconda2/bin/python '
		self.java_path = self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin/java'
		self.cmd_path = self.config.SOFTWARE_DIR + '/ctMethTracer.py'
		# self.script_path = "bioinfo/meta/scripts/beta_diver.sh"
		self.R_path = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/bin/R')
		self.bedtools_path=self.config.SOFTWARE_DIR+'/bioinfo/rna/bedtools2-master/bin/bedtools'
		self.db_addr = '10.100.200.131'
		self.port = 27017
		self.data_dir = self.config.SOFTWARE_DIR+'/database/ctMethTracer/'
		# self.stable = self.data_dir+'stableSitesMean_33_Normal.txt'
		# self.normal = self.data_dir+'normalMode.txt'
		# self.bed450K = self.data_dir+'450K_probeInfo_intersect.bed'
		# self.type2class = self.data_dir+'type2class.txt'
		# self.trainSet = self.data_dir+'cancer14_intersect.xls'
		# self.ref = self.data_dir+'position_ref.bed'
		# self.regionRef = self.data_dir+'region_cor0.5_cpg5_mdist0.002226.txt'
		# self.trainMat = self.data_dir+'normal_cancer14.xls_trans.xls'
		self.cancerLocator = self.data_dir+'dyh.jar'
		self.outdir = os.path.join(self.work_dir,'output')

	def run(self):
		super(MethTracerTool, self).run()
		self.run_MethTracerTool()
		self.end()


	def run_MethTracerTool(self):
		
		#log_file = self.work_dir + '/cmd.log'
		#_cmd = self.Python_path + self.cmd_path + ' -jobID %s -log %s' % (self.option('jobID'), log_file)
		_cmd = self.Python_path + self.cmd_path + ' -jobID %s -DBaddress %s -port %s -nthreads %i -bedtoolsPath %s -cancerLocatorPath %s -javaPath %s -outdir %s' % (self.option('jobID'),self.db_addr, self.port, self._cpu, self.bedtools_path, self.cancerLocator , self.java_path,self.outdir) #modified by luotong 20171129
		#_cmd = self.Python_path + self.cmd_path + ' -jobID %s -DBaddress %s -port %s -stable %s -normal %s -bed450K %s -type2class %s -trainSet %s -ref %s -nthreads %i -bedtoolsPath %s -regionRef %s -trainMat %s -cancerLocatorPath %s -javaPath %s -outdir %s' % (self.option('jobID'),self.db_addr, self.port, self.stable, self.normal, self.bed450K, self.type2class, self.trainSet, self.ref, self._cpu, self.bedtools_path, self.regionRef, self.trainMat, self.cancerLocator , self.java_path,self.outdir) #modified by luotong 20171027
		#_cmd = self.Python_path + self.cmd_path + ' -jobID %s' % self.option('jobID')
		#+ ' -log {} '.format("cmd1.log") # modify by hongdongxuan 20170922
	
		self.logger.info("run ctMethTracer.py")
		self.logger.info(_cmd)
		_cmd1 = self.add_command("meth_tracer.sh", _cmd).run()

		self.wait(_cmd1)
		if(_cmd1.return_code == 0):
			self.logger.info("run ctMethTracer successfully.")
		else:
			raise('ctMethTracer is down.')
		
		# _path = self.option("path")
		# _model = self.option("model")
		# _id = self.option('jobID')
		# _client = MongoClient(self.db_addr,self.port)
		#_db = _client.ctMethTracer
		#_job = _db.job
		#_posts = {
		#		  'step':0,
		#		  'status':'test',
		#		  'description':None,
		#		  'params':{
		#				'input':_path,
		#				'model':_model,
		#		  },
		#		  'email':'123@127.com'
		#	  }
		#_jobID = _job.insert_one(_posts).inserted_id

