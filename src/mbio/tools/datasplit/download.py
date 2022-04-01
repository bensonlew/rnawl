# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20180820
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
from biocluster.api.file.lib.s3 import S3TransferManager


class DownloadAgent(Agent):
    """
    拆分将文件从对象存储上下载合并，放到各产品线的服务器上
    """
    def __init__(self, parent=None):
        super(DownloadAgent, self).__init__(parent)
        options = [
            {"name": "download_file", "type": "string"},  # 参数文件
            {"name": "target_path", "type": "string"},  # 下载的文件所放的路径
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("download_file"):
            raise OptionError("缺少下载的参数文件")
        if not self.option("target_path"):
            raise OptionError("缺少下载的文件所放的路径")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"


class DownloadTool(Tool):
    def __init__(self, config):
        super(DownloadTool, self).__init__(config)
        self.to_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/ob_storage/to_path"

    def get_info(self):
        """
        获取下载信息,并将文件下载或链接到目标路径下，在此不进行合并
        """
        self.specimen_list = []
        self.project_info = {}
        with open(self.option("download_file"), "r") as f, open("list.txt", "w") as w:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                project_type = item[-3]
                project_sn = item[0]
                specimen_name = item[4]
                from_path = item[-1]
                to_path = "output/" + project_sn + "/" + os.path.basename(item[-1])
                self.logger.info(from_path)
                self.logger.info(to_path)
                if from_path.startswith("s3:"):
                    self.s3_download(from_path, to_path)
                else:
                    self.link(from_path, os.path.join(self.work_dir, to_path))
                w.write(specimen_name + "\t" + os.path.join(self.work_dir, to_path) + "\n")

    def s3_download(self, from_path, to_path):
        """
        对象存储上的文件的下载
        """
        transfer = S3TransferManager(base_path=self.work_dir)
        transfer.add(from_uri=from_path, to_uri=to_path)
        transfer.wait()

    def link(self, old, new):
        """
        文件链接
        """
        dir = os.path.dirname(new)
        if not os.path.exists(dir):
            os.makedirs(dir)
        os.link(old, new)

    def run_download(self):
        for project_sn in self.project_info.keys():
            project_dir = os.path.join(self.output_dir, project_sn)
            for info in self.project_info[project_sn].keys():
                specimen_name = info.split(":")[-2]
                paths = self.project_info[project_sn][info]
                if len(paths) == 1:
                    from_files = paths[0]
                    # from_files = os.path.dirname(paths[0]) + "/"
                    # from_files = "s3://rerewrweset/files/datasplit/2018/20180703PE300-M1/5b7b9b0cf9f24c8d188b4579_20180821_125458/test/"
                    to_path = "output/" + project_sn + "/" + os.path.basename(paths[0])
                    self.logger.info(from_files)
                    self.logger.info(to_path)
                    # self.download_from_s3(from_files, to_path=to_path, cover=True)
                    transfer = S3TransferManager(base_path=self.work_dir)
                    transfer.add(from_uri=from_files, to_uri=to_path)
                    transfer.wait()
                else:
                    self.logger.info("样本:{}需要合并".format(specimen_name))
        # download_from_s3(from_files, to_path="download/", cover=True)
        # from_path = "s3://rerewrweset/files/datasplit/2018/20180703PE300-M1/5b7b9b0cf9f24c8d188b4579_20180821_125458/meta_qc/MJ180629_10:1:5b7b9b0df9f24c8d188b458a.fq"
        # to_path = "output/MJ180629_10:1:5b7b9b0df9f24c8d188b458a.fq"
        # transfer = S3TransferManager(base_path=self.work_dir)
        # transfer.add(from_uri=from_path, to_uri=to_path)
        # transfer.wait()

    def run(self):
        super(DownloadTool, self).run()
        self.get_info()
        # self.run_download()
        self.end()
