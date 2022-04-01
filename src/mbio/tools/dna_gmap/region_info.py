
# -*- coding: utf-8 -*-
# __author__ = 'zegnjing'
# modified 2018.07.09

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import json
import math


class RegionInfoAgent(Agent):
    """
    关联区域详情
    """
    def __init__(self, parent):
        super(RegionInfoAgent, self).__init__(parent)
        options = [
            {"name": "vcf_total", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "trait_file", "type": "infile", "format": "dna_gmap.trait"},  # 性状文件
            {"name": "pop_type", "type": "string", "default": "F2"},  # 群体类型
            {"name": "pid", "type": "string"},  # 父本
            {"name": "mid", "type": "string"},  # 母本
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("vcf_total").is_set:
            raise OptionError("请设置性状的vcf.total文件", code="34801701")
        if not self.option("trait_file").is_set:
            raise OptionError("请设置性状文件", code="34801702")
        if not self.option("pid"):
            raise OptionError("请设置父本", code="34801703")
        if not self.option("mid"):
            raise OptionError("请设置母本", code="34801704")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(RegionInfoAgent, self).end()


class RegionInfoTool(Tool):
    def __init__(self, config):
        super(RegionInfoTool, self).__init__(config)
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.region_info = self.config.PACKAGE_DIR + "/dna_gmap/trit-markerinfo.pl"

    def run_region_info(self):
        self.region_info_dir = os.path.join(self.work_dir, "region_info")
        if os.path.exists(self.region_info_dir):
            os.system("rm -r {}".format(self.region_info_dir))
        os.mkdir(self.region_info_dir)
        cmd = "{} {} -vcf {} ".format(self.perl_path, self.region_info, self.option("vcf_total").prop["path"])
        cmd += "-trit {} -fid {} -mid {}".format(self.option("trait_file").prop["path"], self.option("pid"), self.option("mid"))
        cmd += " -pop {} -out {}".format(self.option("pop_type"), self.region_info_dir)
        command = self.add_command("region_info", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("脚本trit-markerinfo.pl运行完成")
        else:
            self.set_error("脚本trit-markerinfo.pl运行失败", code="34801701")

    def run_box_stat(self):
        """
        统计关联区域箱线图的数据
        """
        # trait = os.path.basename(self.option("trait_file").prop["path"]).split(".txt")[0]
        with open(self.option("trait_file").prop["path"], "r") as f:
            lines = f.readlines()
            item = lines[1].strip().split(",")
            trait = item[0]
        with open(self.output_dir + "/" + trait + ".box.xls", "w") as w:
            w.write("MarkerID\tGenotype\tMax\tQ3\tQ2\tQ1\tMin\tFilter\n")
            for f in os.listdir(self.region_info_dir):
                old = os.path.join(self.region_info_dir, f)
                if f.endswith(".assocation.xls"):
                    new = os.path.join(self.output_dir, f)
                    if os.path.exists(new):
                        os.remove(new)
                    os.link(old, new)
                else:
                    m = re.match(r"{}-(.+)-info.xls".format(trait), f)
                    mark_id = m.group(1)
                    box_dict = {}
                    with open(old, "r") as f:
                        lines = f.readlines()
                        for line in lines[1:]:
                            item = line.strip().split("\t")
                            if item[1] == "--":
                                continue
                            if math.isnan(float(item[3])):
                                continue
                            if item[1] not in box_dict.keys():
                                box_dict[item[1]] = {"specimen_ids": [], "trait_value": []}
                            box_dict[item[1]]["specimen_ids"].append(item[0])
                            box_dict[item[1]]["trait_value"].append(float(item[3]))
                    for genotype in box_dict.keys():
                        specimen_id = box_dict[genotype]["specimen_ids"]
                        max_box, q3, q2, q1, min_box, filter_list = self.get_box_new(box_dict[genotype]["trait_value"])
                        w.write(mark_id + "\t" + genotype + "\t" + str(min_box) + "\t" + str(q1) + "\t" + str(q2) + "\t" + str(q3) + "\t")
                        w.write(str(max_box) + "\t")
                        if filter_list:
                            filter_dict = {}
                            for i in filter_list.keys():
                                filter_dict[specimen_id[i]] = filter_list[i]
                            w.write(json.dumps(filter_dict) + "\n")
                        else:
                            w.write("\t\n")

    def get_box_new(self, mylist):
        """
        计算四分位数和溢出值
        """
        mylist_ = sorted(mylist)
        length = len(mylist_)
        filter_list = {}
        if length > 3:
            q1_index = 0.25 * (length + 1)
            q2_index = 0.5 * (length + 1)
            q3_index = 0.75 * (length + 1)
            q1_index_int = int(q1_index)
            q2_index_int = int(q2_index)
            q3_index_int = int(q3_index)
            q1 = mylist_[q1_index_int - 1] + (mylist_[q1_index_int] - mylist_[q1_index_int - 1]) * (q1_index - q1_index_int)
            q3 = mylist_[q3_index_int - 1] + (mylist_[q3_index_int] - mylist_[q3_index_int - 1]) * (q3_index - q3_index_int)
            q2 = mylist_[q2_index_int - 1] + (mylist_[q2_index_int] - mylist_[q2_index_int - 1]) * (q2_index - q2_index_int)
            qd = q3 - q1
            max_limit = 1.5 * qd + q3
            min_limit = q1 - 1.5 * qd
            max_box, min_box = mylist_[0], mylist_[0]
            for i in range(length):
                data = mylist[i]
                if data >= min_limit and data <= max_limit:
                    if data > max_box:
                        max_box = data
                    if data < min_box:
                        min_box = data
                else:
                    filter_list[i] = data
            return max_box, q3, q2, q1, min_box, filter_list
        elif(length == 3):
            max_box = mylist_[2]
            min_box = mylist_[0]
            q3 = (mylist_[1] + mylist_[2]) / 2
            q2 = mylist_[1]
            q1 = mylist_[0]
            return max_box, q3, q2, q1, min_box, filter_list
        elif(length == 2):
            max_box = mylist_[1]
            min_box = mylist_[0]
            q2 = (mylist_[1] + mylist_[0]) / 2
            q3 = (mylist_[1] + q2) / 2
            q1 = (q2 + mylist_[0]) / 2
            return max_box, q3, q2, q1, min_box, filter_list
        elif(length == 1):
            return mylist_[0], mylist_[0], mylist_[0], mylist_[0], mylist_[0], filter_list

    def run(self):
        super(RegionInfoTool, self).run()
        self.run_region_info()
        self.run_box_stat()
        self.end()
