# -*- coding: utf-8 -*-
# __author__ = 'wentian'
# modified 2019.02.25

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class SvStatV2Agent(Agent):
    """
    sv 变异类型统计,变异长度统计
    """
    def __init__(self, parent):
        super(SvStatV2Agent, self).__init__(parent)
        options = [
            {"name": "sv_vcf", "type": "infile", "format": "sequence.vcf"}  # vcf
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sv_vcf"):
            raise OptionError("请设置sv注释文件", code="34507001")
        if not os.path.exists(self.option("sv_vcf").prop["path"]):
            raise OptionError("sv注释文件:%s不存在，请检查", variables=(self.option("sv_anno")), code="34507002")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(SvStatV2Agent, self).end()


class SvStatV2Tool(Tool):
    def __init__(self, config):
        super(SvStatV2Tool, self).__init__(config)

    def run_sv_stat(self):
        """
        sv变异类型和长度统计
        """
        stat_type = ["DEL", "INS", "DUP", "INV", "BND"]
        stat_result = {}
        len_result = {}
        sample_len = 0
        samples = []
        with open(self.option("sv_vcf").prop["path"], "r")as fr:
            lines = fr.readlines()
            for line in lines:
                if line.startswith("##"):
                    pass
                elif line.startswith("#CHROM"):
                    samples = line.strip().split("\t")[9: len(line.strip().split("\t"))]
                    sample_len = len(samples)
                    for i in stat_type:
                        stat_result[i] = {}
                        len_result[i] = {}
                        for x in samples:
                            stat_result[i][x] = 0
                            len_result[i][x] = []
                    self.logger.info(len_result)
                else:
                    tmp = line.strip().split("\t")
                    for i in range(sample_len):
                        sample_data = tmp[9 + i]
                        temp = sample_data.strip().split(":")
                        if temp[3] == "PASS" and temp[0] != "0/0":
                            svtype_list = tmp[7].strip().split(";")
                            svtype = svtype_list[1][-3:]
                            stat_result[svtype][samples[i]] += 1
                        m = re.match(r".+;SVTYPE=(.+?);.+END=(\d+);.+", tmp[7])
                        sv_type = m.group(1)
                        end = m.group(2)
                        if temp[3] == "PASS":
                            if sv_type in ["DEL", "DUP", "INV"]:
                                length = int(end) - int(tmp[1])
                                len_result[sv_type][samples[i]].append(int(length))
                            elif sv_type == "INS":
                                n = re.match(r".+;INSLEN=(\d+);.+", tmp[7])
                                length = n.group(1)
                                len_result[sv_type][samples[i]].append(int(length))
                            elif sv_type == "BND":
                                length = len(tmp[3])
                                len_result[sv_type][samples[i]].append(int(length))
        with open(os.path.join(self.output_dir, "stat.txt"), "w")as fw:
            write_lines = "sample\tDEL\tINS\tDUP\tINV\tBND\n"
            for i in samples:
                write_lines += i + "\t"
                for x in stat_type:
                    write_lines += str(stat_result[x][i]) + "\t"
                write_lines += "\n"
            fw.write(write_lines)
        sv_len = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
        sv_len_dict = {}
        self.logger.info(len_result["INS"])
        for i in samples:
            # write_lines = "sv_length\tDEL\tINS\tDUP\tINV\tBND\n"
            write_lines = "style\t0-1k\t1-2k\t2-3k\t3-4k\t4-5k\t5-6k\t6-7k\t7-8k\t8-9k\t9-10k\t>10k\n"
            for l in sv_len:
                sv_len_dict[l] = 0
            with open(os.path.join(self.output_dir, i + ".sv.length.txt"), "w")as fw2:
                for x in stat_type:
                    for y in len_result[x][i]:
                        if int(y) <= 1000:
                            sv_len_dict[0] += 1
                        elif int(y) <= 2000:
                            sv_len_dict[1000] += 1
                        elif int(y) <= 3000:
                            sv_len_dict[2000] += 1
                        elif int(y) <= 4000:
                            sv_len_dict[3000] += 1
                        elif int(y) <= 5000:
                            sv_len_dict[4000] += 1
                        elif int(y) <= 6000:
                            sv_len_dict[5000] += 1
                        elif int(y) <= 7000:
                            sv_len_dict[6000] += 1
                        elif int(y) <= 8000:
                            sv_len_dict[7000] += 1
                        elif int(y) <= 9000:
                            sv_len_dict[8000] += 1
                        elif int(y) <= 10000:
                            sv_len_dict[9000] += 1
                        else:
                            sv_len_dict[10000] += 1
                    write_lines += str(x) + "\t" + str(sv_len_dict[0]) + "\t" + str(sv_len_dict[1000]) + "\t"\
                                   + str(sv_len_dict[2000]) + "\t" + str(sv_len_dict[3000]) + "\t"\
                                   + str(sv_len_dict[4000]) + "\t" + str(sv_len_dict[5000]) + "\t"\
                                   + str(sv_len_dict[6000]) + "\t" + str(sv_len_dict[7000]) + "\t"\
                                   + str(sv_len_dict[8000]) + "\t" + str(sv_len_dict[9000]) + "\t"\
                                   + str(sv_len_dict[10000]) + "\n"
                    for l in sv_len:
                        sv_len_dict[l] = 0
                fw2.write(write_lines)

    def run(self):
        super(SvStatV2Tool, self).run()
        self.run_sv_stat()
        self.end()
