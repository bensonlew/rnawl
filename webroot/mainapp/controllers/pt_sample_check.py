# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
from biocluster.core.function import load_class_by_path
from biocluster.api.file.remote import RemoteFileManager
from mainapp.libs.signature import check_sig
import web
import json
from biocluster.core.exceptions import FileError
import os
from mainapp.models.mongo.submit.paternity_test_mongo import PaternityTest as import_mongodb


class SampleCheck(object):

	@check_sig
	def POST(self):
		data = web.input()
		print "******"
		print data
		print "******"
		client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
		self.sample = []
		self.sample.append(str(data.dad_id))
		self.sample.append(str(data.mom_id))
		self.sample.append(str(data.preg_id))
		for i in self.sample:
			result = import_mongodb().check_sample(i)
			if not result:
				path = data.fastq_path +'/'+ i + '_R1.fastq.gz'
				if not os.path.exists(path):
					info = {"success": False, "info":"数据{}不存在于数据库和fastq文件夹中，请确认该样本已下机".format(i)}
				else:
					info = {"success": True, "info": "数据不在库中，但有fastq文件，现可以转换"}
			else:
				info = {"success": True, "info": "数据在库中，直接提取"}
		return json.dumps(info)



# 测试python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py
		# post pt_sample_check -c client03 -b http://192.168.12.102:8093
		# -n "dad_id;mom_id;preg_id;fastq_path" -d "WQ2193-F;WQ2193-M;WQ2193-S;
		# /mnt/ilustre/users/sanger-dev/sg-users/zhoumoli/pt/Sample_WQ2131"