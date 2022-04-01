# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.03.20

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import json
import subprocess


class CutVcfAgent(Agent):
    """
    生成引物设计需要的输入文件。
    """
    def __init__(self, parent):
        super(CutVcfAgent, self).__init__(parent)
        options = [
            {"name": "marker_type", "type": "string"},  # 标记类型snp/indel,ssr
            {"name": "data_type", "type": "string"},  # 选择标记位置信息 vcf,location(位置文件),custom(自定义)
            {"name": "marker_detail", "type": "string"},  # 选择标记位置信息 location:[{location:chr:start-end},{}] custom:[chr:start-end,....]
            # {"name": "condition", "type": "string"},  # 设置引物设计条件 {snp:{tm1:,tm2:....},indel"{}}
            {"name": "vcf", "type": "string"},  # vcf文件:原始vcf或者运行结果vcf。
            {"name": "ssr", "type": "string"}  # 参考基因组中的ssr.ref.result.xls。
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("marker_type"):
            raise OptionError("请设置marker_type", code="34504401")
        if not self.option("data_type"):
            raise OptionError("请设置data_type", code="34504402")
        if not self.option("data_type") == "vcf" and not self.option("marker_detail"):
            raise OptionError("请设置marker_detail", code="34504404")
        # if not self.option("condition"):
        #     raise OptionError("请设置condition", code="34504406")
        if not self.option("vcf"):
            raise OptionError("请设置vcf", code="34504406")
        if self.option("marker_type") == "ssr" and self.option("data_type") == "custom" and not self.option("ssr"):
            raise OptionError("请设置ssr", code="34504406")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(CutVcfAgent, self).end()


class CutVcfTool(Tool):
    def __init__(self, config):
        super(CutVcfTool, self).__init__(config)
        if self.option("ssr"):
            self.json_path = self.config.SOFTWARE_DIR + "/database/dna_geneome"  # 参考组配置文件
            self.ssr = os.path.join(self.json_path, self.option("ssr"))

    def location(self):
        """
        根据位置信息筛选vcf
        marker_detail:location:[{location:chr:start-end},{}] custom:[chr:start-end,....]
        """
        write_lines = "chr\tpos\tref\talt\n"
        with open(self.option('vcf'), "r")as fr:
            lines = fr.readlines()
            for line in lines:
                if line.startswith("#"):
                    pass
                else:
                    detaiil = line.strip().split("\t")
                    if self.option('marker_type') == "snp/indel":
                        vcf_chr = detaiil[0]
                        vcf_pos = detaiil[1]
                        vcf_ref = detaiil[3]
                        vcf_alt = detaiil[4]
                    else:
                        vcf_chr = detaiil[0]
                        vcf_pos = detaiil[1]
                        vcf_ref = detaiil[4]
                        vcf_alt = detaiil[5]
                    for i in json.loads(self.option('marker_detail')):
                        chr_start_end = i
                        if self.option('data_type') == "location":
                            chr_start_end = i["location"]
                        temp = chr_start_end.strip().split(":")
                        chr = temp[0]
                        tmp = temp[1].strip().split("-")
                        start = tmp[0]
                        end = tmp[1]
                        if vcf_chr == chr and int(start) <= int(vcf_pos) <= int(end):
                            write_lines += vcf_chr + "\t" + vcf_pos + "\t" + vcf_ref + "\t" + vcf_alt + "\n"
                            break
        with open(os.path.join(self.work_dir, "vcf.txt"), "w")as fw:
            fw.write(write_lines)

    def cut_vcf(self):
        """
        如果marker_type是snp/indel，需要把vcf分成两个文件。
        :return:
        """
        snp_write_lines = "chr\tpos\tref\talt\n"
        indel_write_lines = "chr\tpos\tref\talt\n"
        if self.option('data_type') == 'vcf':
            with open(self.option('vcf'), "r")as fr:
                lines = fr.readlines()
                for line in lines:
                    if line.startswith("#"):
                        pass
                    else:
                        temp = line.strip().split("\t")
                        conbine = ",".join([temp[3], temp[4]])
                        tmp = conbine.strip().split(",")
                        line_type = "snp"
                        for i in tmp:
                            if len(i) > 1:
                                line_type = "indel"
                        if line_type == "indel":
                            indel_write_lines += temp[0] + "\t" + temp[1] + "\t" + temp[3] + "\t" + temp[4] + "\n"
                        if line_type == "snp":
                            snp_write_lines += temp[0] + "\t" + temp[1] + "\t" + temp[3] + "\t" + temp[4] + "\n"
        else:
            with open(os.path.join(self.work_dir, "vcf.txt"), "r")as fr:
                lines = fr.readlines()
                for line in lines[1:]:
                    temp = line.strip().split("\t")
                    conbine = ",".join([temp[2], temp[3]])
                    tmp = conbine.strip().split(",")
                    line_type = "snp"
                    for i in tmp:
                        if len(i) > 1:
                            line_type = "indel"
                    if line_type == "indel":
                        indel_write_lines += temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + temp[3] + "\n"
                    if line_type == "snp":
                        snp_write_lines += temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + temp[3] + "\n"
        with open(os.path.join(self.output_dir, "snp_infile.txt"), "w")as fw1:
            fw1.write(snp_write_lines)
        with open(os.path.join(self.output_dir, "indel_infile.txt"), "w")as fw2:
            fw2.write(indel_write_lines)

    def make_infile(self):
        """
        将ssr的vcf转换成输入文件。
        :return:
        """
        write_lines = "chr\tpos\tref\talt\n"
        with open(self.option('vcf'), "r")as fr:
            lines = fr.readlines()
            for line in lines:
                if line.startswith("#"):
                    pass
                else:
                    temp = line.strip().split("\t")
                    write_lines += temp[0] + "\t" + temp[1] + "\t" + temp[4] + "\t" + temp[5] + "\n"
        with open(os.path.join(self.output_dir, "infile.txt"), "w")as fw:
            fw.write(write_lines)

    def find_ssr(self):
        """
        根据参考基因组的ssr文件找到vcf中的ssr，再根据自定义位置信息筛选。
        :return:
        """
        write_lines = "chr\tpos\tref\talt\n"
        ssr_lines_list = []
        with open(self.ssr, "r")as fs:
            ssr_lines = fs.readlines()
            for ssr_line in ssr_lines[1:]:
                ssr_line_tmp = ssr_line.strip().split("\t")
                for i in json.loads(self.option('marker_detail')):
                    chr_start_end = i
                    temp = chr_start_end.strip().split(":")
                    chr = temp[0]
                    tmp = temp[1].strip().split("-")
                    start = tmp[0]
                    end = tmp[1]
                    if chr == ssr_line_tmp[0] and int(start) <= int(ssr_line_tmp[5]) and int(ssr_line_tmp[6]) <= int(end):
                        ssr_lines_list.append(ssr_line)
                        break
        self.logger.info(ssr_lines_list)
        with open(self.option('vcf'), "r")as fr:
            lines = fr.readlines()
            for line in lines:
                if line.startswith("#"):
                    pass
                else:
                    detaiil = line.strip().split("\t")
                    vcf_chr = detaiil[0]
                    vcf_pos = detaiil[1]
                    vcf_ref = detaiil[3]
                    vcf_alt = detaiil[4]
                    for i in json.loads(self.option('marker_detail')):
                        chr_start_end = i
                        temp = chr_start_end.strip().split(":")
                        chr = temp[0]
                        tmp = temp[1].strip().split("-")
                        start = tmp[0]
                        end = tmp[1]
                        if vcf_chr == chr and int(start) <= int(vcf_pos) <= int(end):
                            for ssr_line in ssr_lines_list:
                                ssr_temp = ssr_line.strip().split("\t")
                                if vcf_chr == ssr_temp[0] and int(ssr_temp[5]) <= int(vcf_pos) <= int(ssr_temp[6]):
                                    write_lines += vcf_chr + "\t" + vcf_pos + "\t" + vcf_ref + "\t" + vcf_alt + "\n"
        with open(os.path.join(self.output_dir, "infile.txt"), "w")as fw:
            fw.write(write_lines)

    def run(self):
        super(CutVcfTool, self).run()
        if self.option('data_type') == 'vcf':
            if self.option('marker_type') == "snp/indel":
                self.cut_vcf()
            else:
                self.make_infile()
        elif self.option('data_type') == 'location':
            self.location()
            if self.option('marker_type') == "snp/indel":
                self.cut_vcf()
            else:
                os.link(os.path.join(self.work_dir, "vcf.txt"), os.path.join(self.output_dir, "infile.txt"))
        elif self.option('data_type') == 'custom' and self.option('marker_type') == "snp/indel":
            self.location()
            self.cut_vcf()
        elif self.option('data_type') == 'custom' and self.option('marker_type') == "ssr":
            self.find_ssr()
            self.logger.info(self.option('vcf'))
            self.logger.info(self.option('ssr'))
        self.end()
