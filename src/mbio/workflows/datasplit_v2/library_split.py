# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190128

import re
import os
import json
import time
# from src.mbio.workflows.datasplit_v2.submit import Submit
# from submit import Submit
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from upload_s5cmd import UploadS5cmd


class LibrarySplitWorkflow(Workflow):
    """
    文库拆分，bcl2fastq
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LibrarySplitWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "library_params", "type": "infile", "format": "datasplit.library_params", "required": True},  # 进行文库参数的参数文件
            {"name": "run_type", "type": "string", "default": "auto"},  # 是否进行自动拆分,自动进行下面的样本拆分
            {"name": "split_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.end_times = 0

    def check_options(self):
        if not self.option("library_params").is_set:
            raise OptionError("缺少文库拆分参数，请检查")
        f = open(self.option("library_params").prop["path"], "rb")
        try:
            json_dict = json.loads(f.read())
            self.params_json = json_dict["library_split"]
        except:
            raise OptionError("json格式不正确")

    def run_bcl2fastq(self):
        for lane in self.params_json.keys():
            options = {
                "data_path": self.params_json[lane]["data_path"],
                "sample_sheet": self.params_json[lane]["sample_sheet"],
                "bases_mask": self.params_json[lane]["bases_mask"],
                "barcode_mismatch": self.params_json[lane]["barcode_mismatch"]
            }
            self.bcl2fastq = self.add_tool("datasplit_v2.bcl2fastq")
            self.bcl2fastq.set_options(options)
            self.bcl2fastq.on("end", self.set_output, lane)
            self.bcl2fastq.run()

    def links(self, old, new, meta_genomic_list):
        if os.path.isfile(old):
            if os.path.exists(new):
                os.remove(new)
            os.link(old, new)
        else:
            sample_id = os.path.basename(old)
            # if sample_id in meta_genomic_list:
            #     return
            if not os.path.exists(new):
                os.mkdir(new)
            for f in os.listdir(old):
                self.links(os.path.join(old, f), os.path.join(new, f), meta_genomic_list)

    def set_output(self, event):
        obj = event["bind_object"]
        fastq_dir = os.path.join(obj.output_dir, "Fastq")
        output_dir = os.path.join(self.output_dir, event["data"])
        self.file_rename(fastq_dir)
        time.sleep(60)
        meta_genomic_list = []  # 若是宏基因组的样本，则不进行上传
        # meta_genomic_file = os.path.join(self.work_dir, "meta_genomic.txt")
        # if os.path.exists(meta_genomic_file):
        #     f = open(meta_genomic_file, "rb")
        #     meta_genomic_info = json.loads(f.read())
        #     meta_genomic_list = meta_genomic_info[event["data"]]
        self.links(os.path.join(obj.output_dir, "Fastq"), output_dir, meta_genomic_list)
        time.sleep(60)
        s3_upload_dir = os.path.join(self._sheet.output, event["data"])
        self.set_db(obj.output_dir, output_dir, self.params_json[event["data"]]["sample_sheet"], s3_upload_dir)
        # self.set_db(obj.output_dir, fastq_dir, self.params_json[event["data"]]["sample_sheet"], s3_upload_dir)
        self.end_times += 1
        if self.end_times == len(self.params_json.keys()):
            # self.set_db(self.output_dir, "path")
            # if self.option("run_type") == "auto":  # 不在此时post，因为此时可能文件还未上传完成
            #     self.post_sample_split()
            # rename_file = os.path.join(self.work_dir, "rename.txt")
            # if os.path.exists(rename_file):
            #     f = open(rename_file, "rb")
            #     rename_info = json.loads(f.read())
            #     for lane in rename_info.keys():
            #         lane_dir = os.path.join(self.output_dir, lane)
                    # for s_info in rename_info[lane]:
                        # sample_dir = os.path.join(self.output_dir, lane+"/"+s_info["sample_id"])
                        # old_new = {}
                        # for fq in os.listdir(sample_dir):
                        #     if fq.endswith("R1_001.fastq.gz"):
                        #         old_new[fq] = s_info["library_number"]+":"+s_info["specimen_name"]+".R1.raw.fastq.gz"
                        #         os.rename(os.path.join(sample_dir, fq), os.path.join(sample_dir, s_info["library_number"]+":"+s_info["specimen_name"]+".R1.raw.fastq.gz"))
                        #     elif fq.endswith("R2_001.fastq.gz"):
                        #         old_new[fq] = s_info["library_number"]+":"+s_info["specimen_name"]+".R2.raw.fastq.gz"
                        #         os.rename(os.path.join(sample_dir, fq), os.path.join(sample_dir, s_info["library_number"]+":"+s_info["specimen_name"]+".R2.raw.fastq.gz"))
                        # md5_path = os.path.join(sample_dir, "md5sum.txt")
                        # if os.path.exists(md5_path):
                        #     md5_info = {}
                        #     with open(md5_path, "rb") as r:
                        #         for line in r:
                        #             item = line.strip().split("  ")
                        #             if item[1] in old_new.keys():
                        #                 md5_info[old_new[item[1]]] = item[0]
                        #             else:
                        #                 md5_info[item[1]] = item[0]
                        #     with open(md5_path, "wb") as w:
                        #         for f in md5_info.keys():
                        #             w.write(md5_info[f] + "  " + f + "\n")
                    # sample_sheet = os.path.join(self.work_dir, lane+".sample_sheet.csv")
                    # s3_upload_dir = os.path.join(self._sheet.output, lane)
                    # self.update_path(lane_dir, sample_sheet, s3_upload_dir)
            self.end()

    def file_rename(self, bcl_dir):
        rename_file = os.path.join(self.work_dir, "rename.txt")
        if os.path.exists(rename_file):
            f = open(rename_file, "rb")
            rename_info = json.loads(f.read())
            for lane in rename_info.keys():
                for s_info in rename_info[lane]:
                    sample_dir = os.path.join(bcl_dir, s_info["sample_id"])
                    if not os.path.exists(sample_dir):
                        continue
                    old_new = {}
                    for fq in os.listdir(sample_dir):
                        if fq.endswith("R1_001.fastq.gz"):
                            old_new[fq] = s_info["library_number"]+"--"+s_info["specimen_name"]+".R1.raw.fastq.gz"
                            os.rename(os.path.join(sample_dir, fq), os.path.join(sample_dir, s_info["library_number"]+"--"+s_info["specimen_name"]+".R1.raw.fastq.gz"))
                        elif fq.endswith("R2_001.fastq.gz"):
                            old_new[fq] = s_info["library_number"]+"--"+s_info["specimen_name"]+".R2.raw.fastq.gz"
                            os.rename(os.path.join(sample_dir, fq), os.path.join(sample_dir, s_info["library_number"]+"--"+s_info["specimen_name"]+".R2.raw.fastq.gz"))
                    md5_path = os.path.join(sample_dir, "md5sum.txt")
                    if os.path.exists(md5_path):
                        md5_info = {}
                        with open(md5_path, "rb") as r:
                            for line in r:
                                item = line.strip().split("  ")
                                if item[1] in old_new.keys():
                                    md5_info[old_new[item[1]]] = item[0]
                                else:
                                    md5_info[item[1]] = item[0]
                        with open(md5_path, "wb") as w:
                            for f in md5_info.keys():
                                w.write(md5_info[f] + "  " + f + "\n")

    def set_db(self, output_dir, fastq_dir, sample_sheet, s3_upload_dir):
        """
        导表，更新文库路径
        """
        if self.option("split_id"):
            datasplit_api = self.api.api("datasplit.datasplit_new")
            for l in os.listdir(os.path.join(output_dir, "Reports/html/")):
                if os.path.isdir(os.path.join(output_dir, "Reports/html/", l)):
                    lane_barcode_path = output_dir + "/Reports/html/" + l + "/all/all/all/laneBarcode.html"
                    lane_path = output_dir + "/Reports/html/" + l + "/all/all/all/lane.html"
                    break
            datasplit_api.add_flowcell_summary(self.option("split_id"), lane_path, lane_barcode_path)
            datasplit_api.update_lib_path(self.option("split_id"), fastq_dir, sample_sheet, s3_upload_dir)
        self.logger.info("文库信息导表完成")

    def update_path(self, fastq_dir, sample_sheet, s3_upload_dir):
        """
        改名之后更新文库和原始样本路径
        """
        if self.option("split_id"):
            datasplit_api = self.api.api("datasplit.datasplit_new")
            datasplit_api.update_lib_path(self.option("split_id"), fastq_dir, sample_sheet, s3_upload_dir)
        self.logger.info("文库和原始样本路径更新完成")

    def post_sample_split(self):
        """
        启动接口，进行样本拆分
        """
        split_params = {"split_id": self.option("split_id"), "data_source": "specimen"}
        if self._sheet.client01:
            type = "sanger"
        else:
            type = "tsg"
        result = Submit(split_params, "/s/datasplit/datasplit_v2", type).webapitest()
        self.logger.info(result)
        if result["success"]:
            self.logger.info("拆分请求成功！")
        else:
            self.logger.info("拆分请求失败！")

    def run(self):
        self.run_bcl2fastq()
        super(LibrarySplitWorkflow, self).run()

    def end(self):
        self.add_upload_dir(self.output_dir)
        # u = UploadS5cmd()
        # u.upload(self.output_dir, self._sheet.output)
        super(LibrarySplitWorkflow, self).end()
