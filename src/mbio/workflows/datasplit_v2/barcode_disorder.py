# -*- coding: utf-8 -*-
# __author__: zengjing
# last_modify: 20210223

import os
import re
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class BarcodeDisorderWorkflow(Workflow):
    """
    barcode错乱：一个文库里每个barcode有多少条reads
    """
    def __init__(self, wsheet_object):
        super(BarcodeDisorderWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 文库文件夹路径
            {"name": "barcode_primer_info", "type": "infile", "format": "datasplit.path"},  # 文库中样本barcode引物信息
            {"name": "split_type", "type": "string", "default": "Pair"},  # 拆分样本序列类型 Pair or Single
            {"name": "library_number", "type": "string"},  # 文库名称
            {"name": "lib_insert_size", "type": "int"},  # 文库插入片段长度
            {"name": "lib_type", "type": "string", "default": "no_official"},  # 文库类型，官方建库和非官方建库
            {"name": "update_info", "type": "string"},
            {"name": "verify_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.fq_list = []  # 文库中的序列R1端/R2端
        self.barcode_seq = {}
        self.barcode_path = ""

    def check_options(self):
        if not self.option("fq_dir"):
            raise OptionError("必须输入文库文件夹对应的路径信息")
        if not self.option("barcode_primer_info"):
            raise OptionError("必须输入barcode_primer_info信息表")
        if not self.option("library_number"):
            raise OptionError("必须设置library_number")

    def get_primer_info(self):
        self.barcode_info = os.path.join(self.work_dir, "sample_barcode_info.txt")
        with open(self.option("barcode_primer_info").prop["path"], "rb") as f, open(self.barcode_info, "wb") as w:
            w.write("#Sample\tBarcode-tag\tFbarcode\tRbarcode\n")
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                w.write(item[0] + "\t" + item[1] + "\t" + item[2] + "\t" + item[4] + "\n")
                self.barcode_seq[item[0]] = {"f_barcode": item[2], "r_barcode": item[4]}
        for fq in os.listdir(self.option("fq_dir").prop["path"]):
            if re.search(r"_R1_", fq):
                self.fq_list.insert(0, os.path.join(self.option("fq_dir").prop["path"], fq))
            elif re.search(r".R1.raw.", fq):
                self.fq_list.insert(0, os.path.join(self.option("fq_dir").prop["path"], fq))
            elif re.search(r"_R2_", fq):
                self.fq_list.append(os.path.join(self.option("fq_dir").prop["path"], fq))
            elif re.search(r".R2.raw.", fq):
                self.fq_list.append(os.path.join(self.option("fq_dir").prop["path"], fq))

    def run_getvaild_bybarcode(self):
        options = ({
            "fq1": self.fq_list[0],
            "fq2": self.fq_list[1],
            "barcode_info": self.barcode_info,
            "lib_name": self.option("library_number"),
        })
        self.getvaild_bybarcode = self.add_tool("datasplit_v2.getvaild_bybarcode")
        self.getvaild_bybarcode.set_options(options)
        self.getvaild_bybarcode.on("end", self.get_barcode_count, "valid")
        self.getvaild_bybarcode.run()

    def get_barcode_count(self, event):
        """
        获取seq2sam.stat文件里每个barcode的reads条数
        """
        self.logger.info("计算seq2sam.stat文件里每个barcode的reads条数")
        obj = event["bind_object"]
        seq2sam = ""
        for f in os.listdir(obj.work_dir):
            if f.endswith(".seq2sam.stat"):
                seq2sam = os.path.join(obj.work_dir, f)
                break
        if seq2sam == "":
            raise Exception("没有找到seq2sam.stat文件")
        barcode_count = {}
        with open(seq2sam, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                if item[1] not in barcode_count:
                    barcode_count[item[1]] = 0
                barcode_count[item[1]] = barcode_count[item[1]] + 1
        self.barcode_path = os.path.join(self.work_dir, "barcode_reads.txt")
        with open(self.barcode_path, "wb") as w:
            w.write("barcode\treads\tf_barcode\tr_barcode\n")
            for barcode in barcode_count:
                w.write(barcode + "\t" + str(barcode_count[barcode]) + "\t" + self.barcode_seq[barcode]["f_barcode"] + "\t")
                w.write(self.barcode_seq[barcode]["r_barcode"] + "\n")
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        if self.option("verify_id"):
            datasplit_api = self.api.api("datasplit.datasplit_new")
            verify_id = self.option("verify_id")
            datasplit_api.add_sg_meta_verify_barcode_detail(verify_id, self.barcode_path)
        self.end()

    def run(self):
        self.get_primer_info()
        self.run_getvaild_bybarcode()
        super(BarcodeDisorderWorkflow, self).run()

    def end(self):
        super(BarcodeDisorderWorkflow, self).end()
