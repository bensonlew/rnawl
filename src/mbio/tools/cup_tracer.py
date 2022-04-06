# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 17:17:18 2017

@author: Luo Tong
last modified 20171219
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
class CupTracerAgent(Agent):
	def __init__(self, parent):
		super(CupTracerAgent,self).__init__(parent)
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
		self._cpu = 4
		self._memory = '8G'
		
	def end(self):
		# modified by Luo Tong, 20171219
		def send_res_email(job_id, to_addr = None, client = MongoClient("10.100.200.131", 27017), data = None):
			if not isinstance(job_id, ObjectId):
				job_id = ObjectId(job_id)
			jobInfo = client.cuptracer.job.find_one({'_id':job_id})
			if not to_addr:
				to_addr = jobInfo['email']
			def _format_addr(s):
				name, addr = parseaddr(s)
				return formataddr(( \
					Header(name, 'utf-8').encode(), \
					addr.encode('utf-8') if isinstance(addr, unicode) else addr))
			
			def get_data(job_id, client = client):
				job_score = client.cuptracer.job_detail.find_one({'job_id':ObjectId(job_id)})
				result = {}
				result['score'] = job_score['score']
				result['predicted_type'] = job_score['predicted_type']
				return result
			if not data:
				data = get_data(job_id, client)

			if to_addr and jobInfo['has_sent_email']==1:
				from_addr = 'cuptracer@sina.com'
				password = 'cuptracer2017'
				smtp_server = 'smtp.sina.com'
				
				## email content
				content = ''
				content += 'Dear CUPtracer user,\r\n\r\n'
				content += 'This email is to inform you that your job {} is finished, and the txt file attached is your result. For an interactive view, you can inquire your result at  http://cuptracer.i-sanger.com.\r\n\r\nCUPtracer Server Team'.format(str(job_id))
				subject = 'Your job {} at CUPtracer finished'.format(str(job_id))
				
				## email attachment, modified by luotong 20180124
				attachment = ''
				attachment += 'predicted_type\t{}\r\n'.format(data['score'][data['predicted_type']]['name'])
				attachment += 'score\t{}\r\n'.format(data['score'][data['predicted_type']]['score'])
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
				mime['Content-Disposition'] = 'attachment; filename="CUPtracer_{}.txt"'.format(str(job_id))
				msg.attach(mime)
				
				## send email
				server = smtplib.SMTP(smtp_server, 25)
				server.set_debuglevel(1)
				server.login(from_addr, password)
				try:
					server.sendmail(from_addr, [to_addr], msg.as_string())
					server.quit()
					client.cuptracer.job.update_one({'_id':ObjectId(job_id)}, {'$set':{'has_sent_email':2}})
				except Exception as note:
					client.cuptracer.job.update_one({'_id':ObjectId(job_id)}, {'$set':{'note':note}})
		try:
			send_res_email(job_id=self.option('jobID'))
		except Exception, e:
			self.logger.info(e)
		super(CupTracerAgent,self).end()

class CupTracerTool(Tool):
	def __init__(self, config):
		super(CupTracerTool, self).__init__(config)
		self._cpu = 4
			
		self.Python_path = '/miniconda2/bin/python '
		self.cmd_path = self.config.SOFTWARE_DIR + '/CUPTracer.py'
		# self.script_path = "bioinfo/meta/scripts/beta_diver.sh"
		self.R_path = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/bin/R')
		self.bedtools_path=self.config.SOFTWARE_DIR+'/bioinfo/rna/bedtools2-master/bin/bedtools'
		self.db_addr = '10.100.200.131'
		self.port = 27017
		self.data_dir = self.config.SOFTWARE_DIR+'/database/CUPtracer/'
		self.outdir = os.path.join(self.work_dir,'output')
		self.train = self.data_dir+'cancer31_trainTrans_5000.xls'
		self.probes = self.data_dir+'probe5000.list'

	def run(self):
		super(CupTracerTool, self).run()
		self.run_CupTracerTool()
		self.end()


	def run_CupTracerTool(self):
		
		#log_file = self.work_dir + '/cmd.log'
		#_cmd = self.Python_path + self.cmd_path + ' -jobID %s -log %s' % (self.option('jobID'), log_file)
		_cmd = self.Python_path + self.cmd_path + ' -jobID %s -DBaddress %s -port %s -nthreads %i -bedtoolsPath %s  -outdir %s -train %s -probes %s' % (self.option('jobID'),self.db_addr, self.port, self._cpu, self.bedtools_path, self.outdir, self.train, self.probes) #modified by luotong 20180122
	
		self.logger.info("run CUPTracer.py")
		self.logger.info(_cmd)
		_cmd1 = self.add_command("cup_tracer.sh", _cmd).run()

		self.wait(_cmd1)
		if(_cmd1.return_code == 0):
			self.logger.info("run CUPTracer successfully.")
		else:
			raise('CUPTracer is down.')
		
		# _path = self.option("path")
		# _model = self.option("model")
		# _id = self.option('jobID')
		# _client = MongoClient(self.db_addr,self.port)
		#_db = _client.CUPTracer
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

