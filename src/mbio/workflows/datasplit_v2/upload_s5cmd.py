# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20201010

import os
import re
import subprocess
import logging

class UploadS5cmd(object):
    """
    将文件上传到对象存储
    """
    def __init__(self, aws_profile="default", http="http://10.5.0.10", cpu=100):
        self.s5cmd_path = "/mnt/ilustre/users/sanger-dev/tsanger/app/s3/s5cmd"
        self.aws_profile = aws_profile  # datasplit
        self.http = http
        self.cpu = cpu
        self.upload_path = {}
        self.logger = logging.getLogger("S5cmdFileTransfer")
        self.connection_reset = 0

    def upload(self, dir_path, s3_dir_path):
        """
        上传结果文件到远程路径
        """
        if not os.path.isdir(dir_path):
            raise Exception("%s路径不是有效的文件夹路径" % dir_path)
        else:
            self._dir_path = os.path.abspath(dir_path)
            if not os.listdir(dir_path):
                # raise Exception("文件夹%s为空，请确认是否已经拷贝？" % dir_path)
                self.logger.info("文件夹%s为空，请确认是否已经拷贝？" % dir_path)
                return
        self.logger.info("开始上传%s到%s" % (dir_path, s3_dir_path))
        os.environ["AWS_SHARED_CREDENTIALS_FILE"] = "/mnt/ilustre/users/sanger-dev/.aws/credentials"
        os.environ["AWS_REGION"] = self.aws_profile
        os.environ["AWS_PROFILE"] = self.aws_profile
        self.get_upload_files(dir_path, s3_dir_path)
        for path in self.upload_path:
            self.s5cmd_upload(path, self.upload_path[path])
        self.logger.info("上传%s到%s完成" % (dir_path, s3_dir_path))

    def get_upload_files(self, dir_path, s3_dir_path):
        for f in os.listdir(dir_path):
            if os.path.isdir(os.path.join(dir_path, f)):
                self.get_upload_files(os.path.join(dir_path, f), os.path.join(s3_dir_path, f))
            else:
                self.upload_path[os.path.join(dir_path, f)] = os.path.join(s3_dir_path, f)
                self.logger.info("添加文件%s到上传列表" % (os.path.join(dir_path, f)))

    def s5cmd_upload(self, path, s3_path):
        """
        使用s5cmd将文件上传到远程对象存储
        """
        cmd = self.s5cmd_path + " --endpoint-url " + self.http + " --no-verify-ssl "
        cmd = cmd + "--numworkers " + str(self.cpu)
        cmd = cmd + " ls " + s3_path
        val = os.system(cmd)
        if val == 0:
            self.logger.info("对象存储上文件%s已存在" % (s3_path))
            return
        cmd = self.s5cmd_path + " --endpoint-url " + self.http + " --no-verify-ssl "
        cmd = cmd + "--numworkers " + str(self.cpu)
        cmd = cmd + " cp " + path + " " + s3_path
        self.logger.info("开始传输文件%s 到 %s, 大小:%s" % (path, s3_path, os.path.getsize(path)))
        try:
            result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
            self.logger.info("File:%s trans success" % (path))
        except subprocess.CalledProcessError as e:
            self.logger.info("File:%s trans: %s" % (path, e.output))
            # if re.search("connection reset by peer", e.output) and self.connection_reset < 3:
            #     self.connection_reset += 1
            #     self.s5cmd_upload(path, s3_path)
            # else:
            #     raise Exception("File:%s trans: %s" % (path, e.output))
            if self.connection_reset < 3:
                self.connection_reset += 1
                self.s5cmd_upload(path, s3_path)
            else:
                raise Exception("File:%s trans: %s" % (path, e.output))

    def s5cmd_download(self, s3_path, path):
        """
        使用s5cmd将文件下载到本地
        """
        cmd = self.s5cmd_path + " --endpoint-url " + self.http + " --no-verify-ssl "
        cmd = cmd + "--numworkers " + str(self.cpu)
        cmd = cmd + " cp " + s3_path + " " + path
        self.logger.info("开始下载文件%s 到 %s" % (s3_path, path))
        try:
            result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
            self.logger.info("File:%s download success" % (path))
        except subprocess.CalledProcessError as e:
            self.logger.info("File:%s download: %s" % (path, e.output))
            if self.connection_reset < 3:
                self.connection_reset += 1
                self.s5cmd_download(s3_path, path)
            else:
                raise Exception("File:%s download: %s" % (s3_path, e.output))
