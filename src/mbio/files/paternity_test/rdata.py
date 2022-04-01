# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'

"""rdata格式文件类（空）"""

from biocluster.iofile import File
from biocluster.core.exceptions import OptionError


class RdataFile(File):
	"""
	txt类
	"""

	def __init__(self):
		super(RdataFile, self).__init__()

	def check(self):
		"""
		检测文件是否满足要求，发生错误时应该触发FileError异常
		:return:
		"""
		if super(RdataFile, self).check():
			return True
		else:
			raise FileError("文件格式错误")
