#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20210618

"""提取三代需要上传的文件"""
from genericpath import exists
import os
import json
import shutil
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import TransferManager


class UploadPacbioAgent(Agent):
    def __init__(self, parent):
        super(UploadPacbioAgent, self).__init__(parent)
        options = [
            {"name": "project_table", "type": "infile", "format": "datasplit.path"},  # 项目table表
            {"name": "lane_library_ids", "type": "string"},    # mongo数据库需要更新的lane_library_id
            {"name": "is_upload", "type": "bool", "default": False},  # 是否上传文件到对象存储
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("project_table").is_set:
            raise OptionError("必须设置table表")
        if self.option("is_upload") and not self.option("lane_library_ids"):
            raise OptionError("is_upload为true的时候必须设置lane_library_ids")

    def set_resource(self):
        self._cpu = "1"
        self._memory = "10G"

    def end(self):
        # self.add_upload_dir(self.output_dir)
        super(UploadPacbioAgent, self).end()


class UploadPacbioTool(Tool):
    def __init__(self, config):
        super(UploadPacbioTool, self).__init__(config)
        self._version = 1.0
        self.python3 = "program/Python35/bin"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.python3, LD_LIBRARY_PATH=self.lib)
        self.get_rawdata = self.config.PACKAGE_DIR + "/datasplit/getRawdataPath.py"
        # self.pbmerge = "bioinfo/ref_rna_v3/HTSeq/miniconda3/bin/pbmerge"
        pbmerge = self.config.SOFTWARE_DIR + "/program/miniconda3/bin/pbmerge"
        if os.path.exists(pbmerge):
            self.pbmerge = "program/miniconda3/bin/pbmerge"
        else:
            self.pbmerge = "bioinfo/ref_rna_v3/HTSeq/miniconda3/bin/pbmerge"
        self.perl = "program/perl-5.24.0/bin/perl"
        self.trim_fq = self.config.PACKAGE_DIR + "/datasplit/trim_fqSeq_pacbio.pl"
        self.fq_rename = self.config.PACKAGE_DIR + "/datasplit/fastqRename.py"

    def get_rawdatas_path(self):
        """
        获取table表对应的path
        """
        self.out_table = os.path.join(self.work_dir, "project_table.path.xls")
        cmd = "{}/python {} -i {} -o {}".format(self.python3, self.get_rawdata,
              self.option("project_table").prop["path"], self.out_table)
        command = self.add_command("get_rawdata_path", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("getRawdataPath运行成功")
        else:
            self.set_db("falied", "getRawdataPath运行失败")
            self.set_error("getRawdataPath运行失败")

    def get_merge_trim_info(self):
        """
        获取需要进行merge和trim、rename的信息
        """
        raw_paths = []
        with open(self.out_table, "rb") as f:
            lines = f.readlines()
            for i in range(1, len(lines)):
                item = lines[i].strip().split("\t")
                raw_path = self.replace_path(item[-5])
                clean_path = self.replace_path(item[-4])
                merge = item[-3]
                trim = item[-2]
                rename = item[-1]
                new_name = ""
                if item[0] != "全长微生物多样性分析":  # modifed by zengjing 20210929 除了“全长微生物多样性分析”项目外，其余项目的名字规则改为：美吉编号_样品名称_文件名
                    new_name = item[4] + "_" + item[5]
                self.path_info[i] = {
                    "raw_path": raw_path,
                    "clean_path": clean_path,
                    "s3_raw_path": "-",
                    "s3_clean_path": "-",
                    "new_name": new_name,
                }
                if merge != "false":
                    self.merge_info[i] = merge
                if trim != "false":
                    self.trim_info[i] = trim
                if rename != "false":
                    self.rename_info[i] = rename
                if raw_path.endswith("fastq"):
                    if i not in self.gz_info.keys():
                        self.gz_info[i] = {}
                    self.gz_info[i]["raw_path"] = raw_path
                if clean_path.endswith("fastq"):
                    if i not in self.gz_info.keys():
                        self.gz_info[i] = {}
                    self.gz_info[i]["clean_path"] = clean_path
                if raw_path != "" and raw_path != "-":
                    if raw_path in raw_paths:
                        self.set_db("falied", "样本%s的raw文件:%s重复，请检查" % (item[4], raw_path))
                        self.set_error("falied", "样本%s的raw文件:%s重复，请检查" % (item[4], raw_path))
                    raw_paths.append(raw_path)

    def new_split_path(self, split_path):
        if os.path.exists(split_path):
            return split_path
        if "ilustre" in split_path:
            split_path1 = split_path.replace("ilustre", "clustre")
            if os.path.exists(split_path1):
                return split_path1
        if "sglustre" in split_path:
            split_path1 = split_path.replace("sglustre", "ilustre")
            if os.path.exists(split_path1):
                return split_path1
        return split_path

    def replace_path(self, old_path):
        if old_path == "-" or old_path == "":
            return "-"
        new_path = []
        old_paths = old_path.split(";")
        for path in old_paths:
            # path1 = path.replace("ilustre", "clustre")
            path1 = self.new_split_path(path)
            if not os.path.exists(path1):
                self.set_db("falied", "文件:%s没有找到，请检查" % path1)
                self.set_error("文件:%s 没有找到，请检查" % path1)
            new_path.append(path1)
        return ";".join(new_path)

    def run_pbbam_merge(self):
        """
        对table表里的细菌基因组完成图需要merge的bam用pbmerge进行merge
        """
        merge_dir = os.path.join(self.work_dir, "merge")
        if not os.path.exists(merge_dir):
            os.mkdir(merge_dir)
        for i in self.merge_info:
            merge_bam = os.path.join(merge_dir, self.merge_info[i])
            input_bam = " ".join(self.path_info[i]["raw_path"].split(";"))
            cmd = "{} -o {} {}".format(self.pbmerge, merge_bam, input_bam)
            command = self.add_command("pbmerge"+str(i), cmd).run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("pbmerge"+str(i)+"样本合并运行成功")
            else:
                self.set_db("falied", "pbmerge"+str(i)+"样本合并运行失败")
                self.set_error("pbmerge"+str(i)+"样本合并运行失败")
            self.path_info[i]["raw_path"] = merge_bam

    def run_trim_fastq(self):
        """
        对table表里的全长微生物多样性分析需要trim和rename进行trim和rename
        """
        trim_dir = os.path.join(self.work_dir, "trim")
        if not os.path.exists(trim_dir):
            os.mkdir(trim_dir)
        for i in self.trim_info:
            trim_fq = os.path.join(trim_dir, self.trim_info[i].split(":")[0])
            mx = self.trim_info[i].split(":")[1].split("-")
            cmd = "{} {} -i {} -o {} -m {} -x {}".format(self.perl, self.trim_fq,
                   self.path_info[i]["raw_path"], trim_fq, mx[0], mx[1])
            command = self.add_command("trim_fq_seq"+str(i), cmd).run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("trim_fqSeq"+str(i)+"长度质控运行成功")
            else:
                self.set_db("falied", "trim_fqSeq"+str(i)+"长度质控运行失败")
                self.set_error("trim_fqSeq"+str(i)+"长度质控运行失败")
            rename = self.rename_info[i].split(":")
            rename_fq = os.path.join(trim_dir, rename[0])
            cmd = "{}/python {} -i {} -o {} -n {}".format(self.python3, self.fq_rename,
                  trim_fq, rename_fq, rename[1])
            command = self.add_command("fastq_rename"+str(i), cmd).run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("fastqRename"+str(i)+"fq重命名运行成功")
            else:
                self.set_db("falied", "fastqRename"+str(i)+"fq重命名运行失败")
                self.set_error("fastqRename"+str(i)+"fq重命名运行失败")
            self.path_info[i]["clean_path"] = rename_fq
            if i not in self.gz_info.keys():
                self.gz_info[i] = {}
            self.gz_info[i]["clean_path"] = rename_fq

    def run_fastq_gz(self):
        """
        对文件进行压缩
        """
        self.logger.info("开始进行文件压缩")
        gz_dir = os.path.join(self.work_dir, "gz_dir")
        if not os.path.exists(gz_dir):
            os.mkdir(gz_dir)
        new_paths = {}
        for i in self.gz_info:
            for k in self.gz_info[i]:
                path = self.gz_info[i][k]
                # gz_path = path + ".gz"
                # if os.path.exists(gz_path):
                #     os.remove(gz_path)
                gz_path = os.path.join(gz_dir, os.path.basename(path) + ".gz")
                if "clustre" in path:
                    new_path = os.path.join(gz_dir, os.path.basename(path))
                    gz_path = os.path.join(gz_dir, os.path.basename(path) + ".gz")
                    if new_path in new_paths.keys():
                        if new_paths[new_path] == path:
                            self.set_db("falied", "{}文件重复，请检查".format(path))
                            self.set_error("{}文件重复，请检查".format(path))
                            raise Exception("{}文件重复，请检查".format(path))
                        else:
                            gz_dir1 = os.path.join(self.work_dir, "gz_dir1")
                            if not os.path.exists(gz_dir1):
                                os.mkdir(gz_dir1)
                            new_path = os.path.join(gz_dir1, os.path.basename(path))
                            gz_path = os.path.join(gz_dir1, os.path.basename(path) + ".gz")
                    new_paths[new_path] = path
                    os.system("cp {} {}".format(path, new_path))
                    path = new_path
                os.system("gzip -c {} > {}".format(path, gz_path))
                if k == "raw_path":
                    self.path_info[i]["raw_path"] = gz_path
                if k == "clean_path":
                    self.path_info[i]["clean_path"] = gz_path
        self.logger.info("文件压缩完成")

    def upload_output_dir(self):
        """
        上传文件
        """
        data_json = os.path.dirname(self.work_dir) + "/data.json"
        self.s3_output_dir = json.loads(open(data_json).read())["output"]
        self.transfer = TransferManager()
        for i in self.path_info:
            raw_path = self.upload_file(self.path_info[i]["raw_path"], self.path_info[i]["new_name"])
            clean_path = self.upload_file(self.path_info[i]["clean_path"], self.path_info[i]["new_name"])
            self.path_info[i]["s3_raw_path"] = raw_path
            self.path_info[i]["s3_clean_path"] = clean_path
            self.run_md5sum(self.path_info[i]["raw_path"],raw_path)
            self.run_md5sum(self.path_info[i]["clean_path"],clean_path)
        self.transfer.wait()

    def run_md5sum(self,path,s3_path):
        """
        生成md5校验码
        """
        if s3_path != "-":
            md5 = os.popen("md5sum {}".format(path)).readlines()[0].split(" ")[0]
            if self.md5_dict.has_key(os.path.dirname(s3_path)):
                self.md5_dict[s3_path][os.path.basename(s3_path)] = md5
            else:
                self.md5_dict[s3_path] = {}
                self.md5_dict[s3_path][os.path.basename(s3_path)] = md5

    def upload_file(self, old_path, new_name):
        if old_path != "-":
            if "Pacbio/" in old_path:
                file_name = old_path.split("Pacbio/")[1]
                if new_name:
                    file_name = os.path.join(os.path.dirname(file_name), new_name + "_" + os.path.basename(file_name))
                new_path = os.path.join(self.s3_output_dir, file_name)
            else:
                file_name = "/".join(old_path.split("/")[-3:])
                if new_name:
                    file_name = os.path.join(os.path.dirname(file_name), new_name + "_" + os.path.basename(file_name))
                new_path = os.path.join(self.s3_output_dir, file_name)
            if new_path in self.new_s3_path:
                self.set_db("falied", "{}文件重复，请检查".format(new_path))
                self.set_error("{}文件重复，请检查".format(new_path))
                raise Exception("{}文件重复，请检查".format(new_path))
            self.new_s3_path.append(new_path)
            # self.logger.info(old_path)
            # self.logger.info(new_path)

            self.transfer.add(from_uri=old_path, to_uri=new_path)
        else:
            new_path = "-"
        return new_path

    def make_md5sum_file(self):
        md5dir = os.path.join(self.work_dir,"md5sum")
        if not os.path.exists(md5dir):
            os.mkdir(md5dir)
        md5_path = os.path.join(md5dir, "md5sum.txt")
        with open(md5_path, "wb") as w:
            for i in self.md5_dict.keys():
            # md5dir = os.path.join(self.work_dir,os.path.basename(i).split('.')[0])
                for f in self.md5_dict[i].keys():
                    w.write(self.md5_dict[i][f] + "  " + f + "\n")
                self.transfer.add(from_uri=md5_path,to_uri=os.path.join(os.path.dirname(i),"md5sum.txt"))
            self.transfer.wait()

    def create_sample_path(self):
        """
        创建项目的样本对应的path
        """
        self.new_out_table = os.path.join(self.output_dir, "project_table.path.xls")
        with open(self.out_table, "rb") as f, open(self.new_out_table, "wb") as w:
            lines = f.readlines()
            # w.write(lines[0])
            w.write(lines[0].split("\n")[0] + "\traw_path\tclean_path\ts3_raw_path\ts3_clean_path\n")
            for i in range(1, len(lines)):
                item = lines[i].strip().split("\t")
                # w.write("\t".join(item[:-5]) + "\t" + self.path_info[i]["raw_path"] + "\t")
                # w.write(self.path_info[i]["clean_path"] + "\t" + "\t".join(item[-5:]) + "\n")
                w.write("\t".join(item) + "\t" + self.path_info[i]["raw_path"] + "\t" + self.path_info[i]["clean_path"] + "\t")
                w.write(self.path_info[i]["s3_raw_path"] + "\t" + self.path_info[i]["s3_clean_path"] + "\n")

    def set_db(self, status, desc=""):
        if self.option("is_upload"):
            api_db = datasplit_api = self.api.api("datasplit.datasplit_upload")
            if status != "end":
                api_db.update_coll_status_pacbio(status=status, lane_library_ids=self.option("lane_library_ids"), desc=desc)
            else:
                api_db.update_coll_path_pacbio(project_table=self.new_out_table,md5sum=self.md5_dict)

    def run(self):
        super(UploadPacbioTool, self).run()
        self.md5_dict = {}
        self.set_db("start")
        self.get_rawdatas_path()
        self.path_info, self.merge_info, self.trim_info, self.rename_info = {}, {}, {}, {}
        self.gz_info = {}
        self.get_merge_trim_info()
        if self.merge_info.keys() != 0:
            self.run_pbbam_merge()
        if self.trim_info.keys() != 0:
            self.run_trim_fastq()
        if self.gz_info.keys() != 0:
            self.run_fastq_gz()
        self.new_s3_path = []
        # self.link_output_dir()
        if self.option("is_upload"):
            self.upload_output_dir()
            self.make_md5sum_file()
        self.create_sample_path()
        self.set_db("end")
        self.end()
