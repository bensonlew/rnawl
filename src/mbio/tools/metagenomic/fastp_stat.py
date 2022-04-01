#-*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import json
from biocluster.core.exceptions import OptionError
from collections import defaultdict

class FastpStatAgent(Agent):
    """
    对fastp软件计算结果进行统计
    """
    def __init__(self, parent):
        super(FastpStatAgent, self).__init__(parent)
        options = [
            {"name": "json_dir", "type": "infile", "format": "meta.qc.json_dir"},
            {"name": "sample_info", "type": "infile", "format": "meta_genomic.specimen_info"},
            {"name": "clean", "type": "bool", "default": False},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("json_dir").is_set:
            raise OptionError("必须提供fastp软件计算json结果文件夹")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "2G"

    def end(self):
        super(FastpStatAgent, self).end()


class FastpStatTool(Tool):
    def __init__(self, config):
        super(FastpStatTool, self).__init__(config)
        self.sample_info = defaultdict(list)

    def stat_info(self):
        """
        统计结果
        :return:
        """
        json_dir = self.option('json_dir').prop['path']
        raw_stat = self.work_dir + "/reads.rawData.stat"
        raw_stat2 = self.work_dir + "/reads.rawData2.stat"
        clean_stat = self.work_dir + "/reads.cleanData.stat"
        with open(raw_stat2, 'w') as w2, open(raw_stat, 'w') as w, open(clean_stat, 'w') as m:
            w.write('#Sample\tReadsNum\tBasesNum\tAverageLength\n')
            w2.write('#Sample\tSource type\tInsert size(bp)\tRead length(bp)\tRaw reads\tRaw bases(bp)\n')
            m.write('#Sample\tReadsNum\tBasesNum\tAverageLength\n')
            for sample_file in os.listdir(json_dir):
                json_path = os.path.join(json_dir, sample_file)
                json_dict = json.load(open(json_path, "r"))
                summary = json_dict["summary"]
                specimen_name = sample_file.rsplit('.')[0]
                if specimen_name not in self.sample_info:
                    sp_info = ['-', '-', summary["before_filtering"]["read1_mean_length"]]
                else:
                    sp_info = self.sample_info[specimen_name]
                raw_average = (int(summary["before_filtering"]["read1_mean_length"]) + int(summary["before_filtering"]["read2_mean_length"])) / 2
                raw_read_num = int(summary["before_filtering"]["total_reads"])
                raw_base = int(summary["before_filtering"]["total_bases"])
                clean_average = (int(summary["after_filtering"]["read1_mean_length"]) + int(summary["after_filtering"]["read2_mean_length"])) / 2
                clean_read_num = int(summary["after_filtering"]["total_reads"])
                clean_base = int(summary["after_filtering"]["total_bases"])
                w.write('{}\t{}\t{}\t{}\n'.format(specimen_name, raw_read_num, raw_base, raw_average))
                w2.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(specimen_name, sp_info[0], sp_info[1], sp_info[2],raw_read_num, raw_base))
                m.write('{}\t{}\t{}\t{}\n'.format(specimen_name, clean_read_num, clean_base, clean_average))

    def get_base_info(self):
        json_dir = self.option('json_dir').prop['path']
        base_info = os.path.join(self.output_dir, "base_info.txt")
        title = ["A_quality", "T_quality", "G_quality", "C_quality", "mean_quality",
                 "A_content", "T_content", "G_content", "C_content", "N_content",
                 "type", "sample"]
        with open(base_info, 'w') as w:
            w.write("pos\t{}\n".format("\t".join(title)))
            for sample_file in os.listdir(json_dir):
                sample = sample_file.split(".json")[0]
                json_path = os.path.join(json_dir, sample_file)
                json_dict = json.load(open(json_path, 'r'))
                for read in ["read1", "read2"]:
                    out_dict = defaultdict(dict)
                    before_f = read + "_before_filtering"
                    for info in ["content", "quality"]:
                        info_t = info + "_curves"
                        for k in json_dict[before_f][info_t]:
                            if k not in ['A', 'T', 'G', 'C', 'N', 'mean']:
                                continue
                            one = json_dict[before_f][info_t][k]
                            for index, n in enumerate(one):
                                out_dict[index]["sample"] = sample
                                out_dict[index]["type"] = read
                                key = k + '_' + info
                                out_dict[index][key] = n
                    for k, v in out_dict.items():
                        w.write("{}".format(k + 1))
                        for t in title:
                            w.write("\t{}".format(v[t]))
                        w.write('\n')

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.logger.info('正在生成文件结果目录')
        if os.path.exists(self.output_dir + "/reads.rawData.stat"):
            os.remove(self.output_dir + "/reads.rawData.stat")
        os.link(self.work_dir + "/reads.rawData.stat", self.output_dir + "/reads.rawData.stat")
        if os.path.exists(self.output_dir + "/reads.cleanData.stat"):
            os.remove(self.output_dir + "/reads.cleanData.stat")
        os.link(self.work_dir + "/reads.cleanData.stat", self.output_dir + "/reads.cleanData.stat")
        self.logger.info('生成文件结果目录成功！')

    def get_sp_info(self):
        with open(self.option("sample_info").path, 'r') as r:
            r.readline()
            for line in r:
                line = line.strip().split('\t')
                self.sample_info[line[0]] = line[1:4]

    def run(self):
        super(FastpStatTool, self).run()
        if self.option("sample_info").is_set:
            self.get_sp_info()
        self.stat_info()
        if not self.option("clean"):
            self.get_base_info()
        self.set_output()
        self.end()

