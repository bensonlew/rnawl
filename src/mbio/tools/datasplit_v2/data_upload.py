# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190327

import os
import commands
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import TransferManager


class DataUploadAgent(Agent):
    """
    数据上传
    """
    def __init__(self, parent=None):
        super(DataUploadAgent, self).__init__(parent)
        options = [
            {"name": "upload_path", "type": "string", "required": True},  # 线下上传文件或文件夹路径
            {"name": "target_path", "type": "string", "required": True},  # 上传的对象存储目标路径
            {"name": "server_type", "type": "string", "default": "other"},  # 上传路径所处的服务器，sanger/other,为other的时候要将ilustre替换成clustre
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.memory_ = 0
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("upload_path"):
            raise OptionError("请设置上传路径upload_path")
        if not self.option("target_path"):
            raise OptionError("请设置上传的对象存储目标路径target_path")
        if not self.option("target_path").startswith("s3://"):
            raise OptionError("请将目标路径target_path设置为对象存储的路径")
        upload_path = self.option("upload_path").replace("ilustre", "clustre") if self.option("server_type") == "other" else self.option("upload_path")
        # if not os.path.exists(upload_path):
        #     raise OptionError("%s上传路径不存在，请检查" % upload_path)
        self.check_gzip(upload_path)
        # if self.memory_ > 200:
        #     if os.path.isdir(upload_path):
        #         raise OptionError("文件夹过于庞大,无法进行压缩,请压缩成一个文件之后再进行上传")
        #     self.memory_ = 200
        if self.memory_ > 100:
            self.memory_ = 100

    def check_gzip(self, path):
        """
        检查文件的压缩情况，不能有未压缩的fasta和fastq文件
        """
        if os.path.isdir(path):
            for f in os.listdir(path):
                try:
                    f_ = os.path.join(path, f)
                except UnicodeDecodeError:
                    raise OptionError("有中文命名的文件，无法解析，请用非中文命名")
                if os.path.isdir(f_):
                    f_ = commands.getoutput("readlink -f {}".format(f_))
                    f_ = f_.replace("ilustre", "clustre") if self.option("server_type") == "other" else f_
                    self.check_gzip(f_)
                else:
                    self.check_size(f_)
        else:
            f_ = path.replace("ilustre", "clustre") if self.option("server_type") == "other" else path
            self.check_size(f_)

    def check_size(self, f_):
        """
        检查文件的大小，进行设置资源
        """
        # if f_.endswith("fasta") or f_.endswith("fastq"):
        #     raise OptionError("上传的文件不能有未压缩的fasta/fastq文件：%s" % f_)
        self.source_file = commands.getoutput("readlink -f {}".format(f_))
        if not self.source_file:
            self.get_source_file(f_)
        else:
            self.source_file = self.source_file.replace("ilustre", "clustre") if self.option("server_type") == "other" else self.source_file
        if self.source_file.startswith("../"):
            num = len(self.source_file.split("../"))
            self.source_file = "/" + "/".join(f_.split("/")[1:-num])
            self.check_size(self.source_file)
        if not os.path.exists(self.source_file):
            raise OptionError("%s文件的源文件没找到，请检查" % f_)
        size = os.path.getsize(self.source_file)
        if size == 0:
            raise OptionError("上传的文件不能有空文件:%s" % f_)
        f_size = float(size) / 1024 / 1024 / 1024
        self.memory_ += f_size

    def get_source_file(self, f):
        """
        循环找源文件
        """
        f_ = commands.getoutput("readlink {}".format(f))
        if not f_:
            self.source_file = f
            return
        else:
            f1 = f_.replace("ilustre", "clustre") if self.option("server_type") == "other" else f_
            self.get_source_file(f1)

    def set_resource(self):
        self._cpu = 2
        self.memory_ = 5 if int(self.memory_) < 5 else int(self.memory_) + 1
        self._memory = str(self.memory_) + "G"

    def end(self):
        super(DataUploadAgent, self).end()


class DataUploadTool(Tool):
    def __init__(self, config):
        super(DataUploadTool, self).__init__(config)
        self.tar_gz_path = "bioinfo/seq/scripts/tar_gz.sh"
        self.islink_status = True
        self.pigz_path = "bioinfo/seq/pigz-2.4/pigz"

    def check_islink(self, old):
        """
        检查上传的文件是否有软链接
        """
        if os.path.islink(old):
            self.islink_status = True
            return True
        if os.path.isdir(old):
            for f in os.listdir(old):
                f_ = os.path.join(old, f)
                f_ = f_.replace("ilustre", "clustre") if self.option("server_type") == "other" else f_
                if os.path.islink(f_):
                    self.islink_status = True
                    return True
                if os.path.isdir(f_):
                    f_ = f_.replace("ilustre", "clustre") if self.option("server_type") == "other" else f_
                    self.check_islink(f_)

    def get_islink(self, old, new):
        """
        将软链接的文件链接到工作路径下
        """
        old_ = commands.getoutput("readlink -f {}".format(old))
        if not old_:
            old_ = os.readlink(old)
        old_ = old_.replace("ilustre", "clustre") if self.option("server_type") == "other" else old_
        if not os.path.exists(new):
            os.mkdir(new)
        if os.path.isdir(old_):
            for f in os.listdir(old_):
                old1 = os.path.join(old_, f)
                self.old1_ = commands.getoutput("readlink -f {}".format(old1))
                if not self.old1_:
                    self.get_source_file(old1)
                else:
                    self.old1_ = self.old1_.replace("ilustre", "clustre") if self.option("server_type") == "other" else self.old1_
                new1 = os.path.join(new, f)
                if os.path.isdir(self.old1_):
                    self.get_islink(self.old1_, new1)
                else:
                    # if os.path.exists(new1):
                    #     self.set_error("文件:%s已经存在，不能有相同名称，请检查" % self.old1_)
                    if self.old1_.startswith("../"):
                        num = len(self.old1_.split("../"))-1
                        self.old1_ = os.path.join("/" + "/".join(old_.split("/")[1:-num]), self.old1_.split("../")[-1])
                        if os.path.isdir(self.old1_):
                            self.get_islink(self.old1_, new1)
                        else:
                            old1 = self.old1_
                            self.old1_ = commands.getoutput("readlink -f {}".format(old1))
                            if not self.old1_:
                                self.get_source_file(old1)
                            else:
                                self.old1_ = self.old1_.replace("ilustre", "clustre") if self.option("server_type") == "other" else self.old1_
                    if os.path.exists(new1):
                        pass
                    else:
                        os.system("cp {} {}".format(self.old1_, new1))
        else:
            # if os.path.exists(new):
            #     self.set_error("文件:%s已经存在，不能有相同名称，请检查" % old_)
            if self.old1_.startswith("../"):
                num = len(self.old1_.split("../"))-1
                self.old1_ = os.path.join("/" + "/".join(old_.split("/")[1:-num]), self.old1_.split("../")[-1])
                old1 = self.old1_
                self.old1_ = commands.getoutput("readlink -f {}".format(old1))
                if not self.old1_:
                    self.get_source_file(old1)
                else:
                    self.old1_ = self.old1_.replace("ilustre", "clustre") if self.option("server_type") == "other" else self.old1_
            if os.path.exists(new1):
                pass
            else:
                os.system("cp {} {}".format(old_, new))

    def get_source_file(self, f):
        """
        循环找源文件
        """
        f_ = commands.getoutput("readlink {}".format(f))
        if not f_:
            self.old1_ = f
            return
        else:
            f1 = f_.replace("ilustre", "clustre") if self.option("server_type") == "other" else f_
            self.get_source_file(f1)

    def tar_gz(self):
        """
        对文件夹进行压缩
        """
        self.gz_path = os.path.join(self.output_dir, os.path.basename(self.upload_path) + ".tar.gz")
        cmd = "{} tar -czvf {} -C {} ".format(self.tar_gz_path, self.gz_path, os.path.dirname(self.upload_path))
        cmd += "{}".format(os.path.basename(self.upload_path))
        self.logger.info(cmd)
        command = self.add_command("tar_gz", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("tar_gz压缩文件运行完成")
        else:
            self.set_error("tar_gz压缩文件运行失败")

    def pigz_gz(self, gz_dir):
        gz_dir1 = os.listdir(gz_dir)
        for i in range(len(gz_dir1)):
            f = os.path.join(gz_dir, gz_dir1[i])
            if os.path.isfile(f):
                if f.endswith("fasta") or f.endswith("fastq") or f.endswith("vcf") or f.endswith("bam"):
                    cmd = "{} {}".format(self.pigz_path, f)
                    command = self.add_command("pigz_file"+str(i), cmd).run()
                    self.wait()
                    if command.return_code == 0:
                        self.logger.info("pigz_file压缩文件运行完成")
                    else:
                        self.set_error("pigz_file压缩文件运行失败")
            else:
                self.pigz_gz(f)

    def get_file_list(self, upload_path, target_path):
        if os.path.isdir(upload_path):
            for f in os.listdir(upload_path):
                f1 = os.path.join(upload_path, f)
                f2 = os.path.join(target_path, f)
                if os.path.isdir(f1):
                    self.get_file_list(f1, f2)
                else:
                    self.upload_files.append(f1)
                    self.target_files.append(f2)
        else:
            self.upload_files.append(upload_path)
            self.target_files.append(target_path)

    def run_upload(self):
        """
        上传文件到对象存储
        """
        self.upload_files, self.target_files = [], []
        # if os.path.isdir(self.upload_path):
        #     self.tar_gz()
        #     self.upload_files.append(self.gz_path)
        #     target_path = os.path.join(self.option("target_path"), os.path.basename(self.gz_path))
        # else:
        #     self.upload_files.append(self.upload_path)
        #     target_path = os.path.join(self.option("target_path"), os.path.basename(self.upload_path))
        # self.target_files.append(target_path)
        if os.path.isdir(self.upload_path):
            self.pigz_gz(self.upload_path)
        self.get_file_list(self.upload_path, os.path.join(self.option("target_path"), os.path.basename(self.upload_path)))
        self.logger.info(self.upload_files)
        self.logger.info(self.target_files)
        transfer = TransferManager()
        for i in range(len(self.upload_files)):
            # transfer.add(from_uri=self.upload_files[i], to_uri=self.target_files[i])
            transfer.add(from_uri=self.upload_files[i], to_uri=self.target_files[i])
            if i / 5 > 0:
                transfer.wait()
                transfer = TransferManager()
        transfer.wait()

    def run(self):
        super(DataUploadTool, self).run()
        _upload_path = self.option("upload_path").replace("ilustre", "clustre") if self.option("server_type") == "other" else self.option("upload_path")
        upload_path = commands.getoutput("readlink -f {}".format(_upload_path))
        if not upload_path:
            self.get_source_file(_upload_path)
            upload_path = self.old1_
        self.check_islink(upload_path)
        # islink_status = self.check_islink(upload_path)
        if self.islink_status:
            self.upload_path = os.path.join(self.output_dir, os.path.basename(upload_path))
            self.get_islink(upload_path, self.upload_path)
        else:
            self.upload_path = upload_path
        self.run_upload()
        self.end()
