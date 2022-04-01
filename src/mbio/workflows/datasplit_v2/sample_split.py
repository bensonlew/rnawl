# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190128

import os
import json
import time
import shutil
# from submit import Submit
import gevent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class SampleSplitWorkflow(Workflow):
    """
    文库拆分成样本，用于多样性和动植物基因组
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SampleSplitWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "project_params", "type": "infile", "format": "datasplit.library_params"},  # 参数文件
            {"name": "run_type", "type": "string", "default": "auto"},  # 是否进行自动拆分,自动进行下面的样本拆分
            {"name": "split_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.start_times, self.end_times = 0, 0

    def check_options(self):
        if not self.option("project_params").is_set:
            raise OptionError("请设置参数文件")

    def get_params(self):
        """
        解析输入json文件，得到其中的参数信息
        """
        f = open(self.option("project_params").prop["path"], "rb")
        try:
            self.json_dict = json.loads(f.read())
        except:
            raise OptionError("json格式不正确")
        return self.json_dict

    def run_dna_split(self):
        for dna_info in self.json_dict["dna"]:
            options = {
                "lib_path": dna_info["lib_path"],
                "library_info": dna_info["library_info"]
            }
            # self.start_times += 1
            self.dna_split = self.add_module("datasplit_v2.dna_split")
            self.dna_split.set_options(options)
            self.dna_split.on("end", self.set_output, "dna")
            self.dna_split.run()

    def run_meta_split(self):
        for meta_info in self.json_dict["meta"]:
            options = {
                "lib_path": meta_info["lib_path"],
                "barcode_info": meta_info["barcode_info"],
                "lib_specimen_id": meta_info["lib_specimen_id"],
            }
            if "mismatch" in meta_info.keys():
                options["mismatch"] = meta_info["mismatch"]
            if "split_type" in meta_info.keys():
                options["split_type"] = meta_info["split_type"]
            if "min_lenth" in meta_info.keys():
                options["min_lenth"] = meta_info["min_lenth"]
            if "max_lenth" in meta_info.keys():
                options["max_lenth"] = meta_info["max_lenth"]
            if "mismatch_rate" in meta_info.keys():
                options["mismatch_rate"] = meta_info["mismatch_rate"]
            if "leading" in meta_info.keys():
                options["leading"] = meta_info["leading"]
            if "sliding_window" in meta_info.keys():
                options["sliding_window"] = meta_info["sliding_window"]
            if "tailing" in meta_info.keys():
                options["tailing"] = meta_info["tailing"]
            if "minlen" in meta_info.keys():
                options["minlen"] = meta_info["minlen"]
            if "valid_len" in meta_info.keys():
                options["valid_len"] = meta_info["valid_len"]
            if "min_len" in meta_info.keys():
                options["min_len"] = meta_info["min_len"]
            # self.start_times += 1
            self.meta_split_qc = self.add_module("datasplit_v2.meta_qc")
            self.meta_split_qc.set_options(options)
            self.meta_split_qc.on("end", self.set_output, "meta")
            self.meta_split_qc.run()

    def run_md5sum(self):
        self.logger.info("开始进行md5校验")
        self.md_tool = []
        for f in os.listdir(self.output_dir):
            # if f in ["meta_qc", "dna_raw"]:
            if f in ["meta_raw", "meta_clean", "dna_raw"]:
                options = {"fastq_dir": os.path.join(self.output_dir, f)}
                self.md5sum = self.add_tool("datasplit_v2.md5sum")
                self.md5sum.set_options(options)
                self.md_tool.append(self.md5sum)
        self.logger.info(len(self.md_tool))
        if len(self.md_tool) > 1:
            self.on_rely(self.md_tool, self.set_db)
            for md_tool in self.md_tool:
                md_tool.run()
        elif len(self.md_tool) == 1:
            self.md_tool[0].on("end", self.set_db)
            self.md_tool[0].run()
        else:
            self.end()

    def link_dir(self, old_dir, new_dir):
        if os.path.exists(new_dir):
            shutil.rmtree(new_dir)
        os.mkdir(new_dir)
        for f in os.listdir(old_dir):
            old = os.path.join(old_dir, f)
            new = os.path.join(new_dir, f)
            if os.path.isdir(old):
                self.link_dir(old, new)
            else:
                if os.path.exists(new):
                    raise Exception("{}文件重复，请检查".format(new))
                os.link(old, new)

    def set_output(self, event):
        obj = event["bind_object"]
        if event["data"] == "dna":
            dna_raw = os.path.join(self.output_dir, "dna_raw")
            self.link_dir(obj.output_dir, dna_raw)
        if event["data"] == "meta":
            # meta_qc = os.path.join(self.output_dir, "meta_qc")
            # self.link_dir(obj.output_dir, meta_qc)
            meta_raw = os.path.join(self.output_dir, "meta_raw")
            meta_clean = os.path.join(self.output_dir, "meta_clean")
            self.link_dir(os.path.join(obj.output_dir, "meta_raw"), meta_raw)
            self.link_dir(os.path.join(obj.output_dir, "meta_clean"), meta_clean)
        time.sleep(10)
        self.end_times += 1
        if self.start_times == self.end_times:
            # if self.option("run_type") == "auto":
            #     self.post_sample_qc()
            # self.set_db()
            self.run_md5sum()
            # self.end()

    def set_db(self):
        """
        更新原始样本路径，多样性质控结果导表
        """
        if self.option("split_id"):
            datasplit_api = self.api.api("datasplit.datasplit_new")
            if "meta" in self.json_dict.keys():
                # meta_qc = os.path.join(self.output_dir, "meta_qc")
                meta_raw = os.path.join(self.output_dir, "meta_raw")
                meta_clean = os.path.join(self.output_dir, "meta_clean")
                fastq_stat = os.path.join(meta_clean, "fastq_stat.xls")
                raw_fastq_stat = os.path.join(meta_raw, "fastq_stat.xls")
                lib_qc_path = os.path.join(meta_clean, "lib_qc_stat.xls")
                trim_hist_dir = os.path.join(meta_clean, "trim_hist")
                s3_upload_dir = os.path.join(self._sheet.output, "meta_clean")
                datasplit_api.update_sample_path(self.option("split_id"), meta_clean, s3_upload_dir, "meta_clean")
                s3_upload_dir = os.path.join(self._sheet.output, "meta_raw")
                datasplit_api.update_sample_path(self.option("split_id"), meta_raw, s3_upload_dir, "meta_raw")
                datasplit_api.add_sg_split_clean_qc(self.option("split_id"), fastq_stat, raw_fastq_stat)
                datasplit_api.add_meta_sg_bar(self.option("split_id"), trim_hist_dir)
                datasplit_api.add_sg_split_library_qc(self.option("split_id"), lib_qc_path)
                datasplit_api.delete_sg_split_clean_qc(self.option("split_id"))
                datasplit_api.delete_sg_split_library_qc(self.option("split_id"))
            if "dna" in self.json_dict.keys():
                dna_raw = os.path.join(self.output_dir, "dna_raw")
                s3_upload_dir = os.path.join(self._sheet.output, "dna_raw")
                datasplit_api.update_sample_path(self.option("split_id"), dna_raw, s3_upload_dir, "dna_raw")
            datasplit_api.update_raw_sample_path(self.option("split_id"))
        self.logger.info("原始样本信息导表完成")
        self.end()

    def post_sample_qc(self):
        """
        启动接口，进行样本质控
        """
        split_params = {"split_id": self.option("split_id"), "data_source": "qc"}
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
        self.get_params()
        if not self.json_dict:
            self.start_listener()
            self.fire("start")
            if self.option("split_id"):
                datasplit_api = self.api.api("datasplit.datasplit_new")
                datasplit_api.update_raw_sample_path(self.option("split_id"))
            gevent.spawn_later(5, self.end)
        else:
            if "dna" in self.json_dict.keys():
                self.start_times += len(self.json_dict["dna"])
            if "meta" in self.json_dict.keys():
                self.start_times += len(self.json_dict["meta"])
            if "dna" in self.json_dict.keys():
                self.run_dna_split()
            if "meta" in self.json_dict.keys():
                self.run_meta_split()
        super(SampleSplitWorkflow, self).run()

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SampleSplitWorkflow, self).end()
