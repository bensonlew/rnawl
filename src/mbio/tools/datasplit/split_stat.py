# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""在备份完成之后，统计一次拆分和二次拆分的结果，生成json"""
from __future__ import division
import os
import re
import json
import errno
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SplitStatAgent(Agent):
    def __init__(self, parent=None):
        super(SplitStatAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'},  # 样本拆分信息表
            {'name': 'bcl2fastq_path', 'type': "string"},  # bcl2fastq的运行目录
            {'name': 'fastx_path', 'type': "string"},  # fastx模块产生的运行目录
            {'name': 'rawValid_path', 'type': "string"},  # rawValid模块的运行目录
            {'name': 'qualcontrol_path', 'type': "string"},  # 质量控制的运行目录
            {'name': 'merge_path', 'type': "string"},  # merge模块的运行目录
            {'name': 'markSeq_path', 'type': "string"},  # mark_seq模块的运行目录
            {'name': 'lenControl_path', 'type': "string"},  # len_control模块的运行目录
            {'name': 'seqExtract_path', 'type': "string"},  # seq_extract模块的运行目录
            {'name': "time", 'type': "infile", 'format': 'datasplit.backup_time'}  # backup时所使用的month和year
        ]
        self.add_option(options)

    def check_option(self):
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        if not self.option("time").is_set:
            raise OptionError("参数time不能为空")
        if not self.option('fastx_path'):
            raise OptionError("参数fastx_path不能为空")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = ''


class SplitStatTool(Tool):
    def __init__(self, config):
        super(SplitStatTool, self).__init__(config)
        self._version = 1.0
        self.backup_dir = "/mnt/ilustre/users/sanger-dev/workspace/datasplit_tmp/"
        self.option('sample_info').get_info()
        self.option('time').get_info()
        self.json_str = ""
        year = self.option('time').prop['year']
        month = self.option('time').prop['month']
        name = "id_" + str(self.option('sample_info').prop["split_id"]) +\
               "_" + str(self.option('sample_info').prop["sequcing_sn"])
        program = self.option('sample_info').prop["program"]
        self.seq_id = os.path.join(self.backup_dir, program, str(year), str(month), name)
        self.statKeyList = ["q20", "q30", "Raw_pair", "chimeric", "chimeric_rate", "NoFBarcode",
                            "NoRBarcode", "valid_pair", "valid_rate", "Pair_trim", "Trim_rate",
                            "Pair_merge", "merge_rate", "primer_mismatch", "final_reads", "split_rate",
                            "highQuality_rate"]

    def make_ess_dir(self):
        dir_list = list()
        statDir = os.path.join(self.work_dir, "stat")
        dir_list = [statDir]
        for name in dir_list:
            try:
                os.makedirs(name)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(name):
                    pass
                else:
                    raise OSError("创建目录失败")

    def createShLog(self):
        """
        将整个过程中调用了其他软件时所使用的sh命令加以整理， 记录成一个文件。 调用了其他文档的模块有：
        bcl2fastq, fastx, qual_control, merge
        """
        bcl2fastqShPath = os.path.join(self.option("bcl2fastq_path"), "sh.log")
        fastxShPath = os.path.join(self.option("fastx_path"), "sh.log")
        qualControlShPath = os.path.join(self.option("qualcontrol_path"), "sh.log")
        mergeShPath = os.path.join(self.option("merge_path"), "sh.log")
        shList = [bcl2fastqShPath, fastxShPath, qualControlShPath, mergeShPath]
        shLogStat = os.path.join(self.work_dir, "stat", "sh.log")
        for shlog in shList:
            with open(shlog, 'rb') as r, open(shLogStat, 'ab') as a:
                for line in r:
                    a.write(line)

    def create_json(self):
        stat_json = dict()
        stat_json["split_id"] = self.option('sample_info').prop["split_id"]
        stat_json["sequcing_id"] = self.option('sample_info').prop["sequcing_id"]
        stat_json["sequcing_sn"] = self.option('sample_info').prop["sequcing_sn"]
        report_path = os.path.join(self.seq_id, "Report")
        stat_json["html_path"] = report_path
        p_sample = list()
        c_sample = list()
        statPath = os.path.join(self.work_dir, "output", 'stat.xls')
        with open(statPath, "wb") as w:
            w.write("文库id\t文库名\tindex\t")
            w.write("\t".join(self.statKeyList))
            w.write("\n")
        for c_id in self.option('sample_info').prop['child_ids']:
            my_c = self.get_child_sample_info(c_id)
            c_sample.append(my_c)
        for p_id in self.option('sample_info').prop['parent_ids']:
            my_p = self.get_parent_sample_info(p_id)
            p_sample.append(my_p)
        stat_json["parent_sample"] = p_sample
        stat_json["child_sample"] = c_sample
        json_path = os.path.join(self.work_dir, 'output', "stat.json")
        self.json_str = json.dumps(stat_json)
        with open(json_path, 'w') as w:
            json.dump(stat_json, w, indent=4)

    def get_parent_sample_info(self, p_id):
        """
        获取一个父样本的拆分信息
        """
        p_info = dict()
        p_stat = dict()
        statPath = os.path.join(self.work_dir, "output", 'stat.xls')
        if self.option('sample_info').parent_sample(p_id, "has_child"):
            p_stat = self.get_p_stat(p_id)
        else:
            p_stat = self.get_no_child_p_stat(p_id)

        with open(statPath, "ab") as a:
            str_ = "{}\t{}\t{}\t".format(p_id, self.option("sample_info").parent_sample(p_id, "library_name"), self.option("sample_info").parent_sample(p_id, "index"))
            a.write(str_)
            tmp = list()
            for k in self.statKeyList:
                tmp.append(str(p_stat[k]))
            a.write("\t".join(tmp))
            a.write("\n")
        fastq1 = dict()
        fastq2 = dict()
        library_type = self.option('sample_info').parent_sample(p_id, "library_type")
        if library_type is None:
            library_type = "undefine"
        parent_path = os.path.join(self.seq_id, library_type, "Library")
        p_info["sample_id"] = p_id
        fastq1 = self.get_parent_fastq_info(parent_path, p_info["sample_id"], "r1")
        fastq2 = self.get_parent_fastq_info(parent_path, p_info["sample_id"], "r2")
        p_info['stat'] = p_stat
        p_info["fastq_info"] = {"fastq1": fastq1, "fastq2": fastq2}
        return p_info

    def get_no_child_p_stat(self, p_id):
        p_stat = dict()
        rawValidStat = os.path.join(self.option("rawValid_path"), "output", "stat.xls")
        with open(rawValidStat, 'rb') as r:
            for line in r:
                line = line.strip("\n")
                line = re.split("\t", line)
                if line[0] == p_id:
                    p_stat["Raw_pair"] = int(float(line[2]))
                    p_stat["valid_pair"] = int(float(line[3]))
        (p_stat["q20"], p_stat["q30"]) = self.get_parent_q20q30(p_id)
        for k in self.statKeyList:
            if k not in ["Raw_pair", "valid_pair", "q20", "q30"]:
                p_stat[k] = "NA"
        return p_stat

    def get_parent_q20q30(self, p_id):
        q20q30Path_r1 = os.path.join(self.option('fastx_path'), "fastx", p_id + "_r1.fastq.q20q30")
        q20q30Path_r2 = os.path.join(self.option('fastx_path'), "fastx", p_id + "_r2.fastq.q20q30")
        with open(q20q30Path_r1, 'rb') as r1, open(q20q30Path_r2, 'rb') as r2:
            q20q30 = list()
            for line in r1:
                line = re.split("\t", line)
                r1_total = int(line[1])
                q20q30.append(int(line[2]))
            r1_q20 = max(q20q30)
            r1_q30 = min(q20q30)

            q20q30 = list()
            for line in r2:
                line = re.split("\t", line)
                r2_total = int(line[1])
                q20q30.append(int(line[2]))
            r2_q20 = max(q20q30)
            r2_q30 = min(q20q30)
            q20 = "{:.2f}%".format(((r1_q20 + r2_q20) / (r1_total + r2_total)) * 100)
            q30 = "{:.2f}%".format(((r1_q30 + r2_q30) / (r1_total + r2_total)) * 100)
        return (q20, q30)

    def get_parent_fastq_info(self, parent_path, sample_id, suffix):
        fastq = dict()
        fastq["file_path"] = os.path.join(parent_path, sample_id + "_" + suffix + ".fastq.gz")
        fastq["quilty_table"] = os.path.join(parent_path, "fastx",
                                             sample_id + "_" + suffix + ".fastq.fastxstat")
        fastq["boxplot_png"] = os.path.join(parent_path, "fastx",
                                            sample_id + "_" + suffix + ".fastq.fastxstat.box.png")
        fastq["nucl_png"] = os.path.join(parent_path, "fastx",
                                         sample_id + "_" + suffix + ".fastq.fastxstat.nucl.png")
        q20q30 = os.path.join(parent_path, "fastx",
                              sample_id + "_" + suffix + ".fastq.q20q30")
        with open(q20q30, 'r') as r:
            line = r.readline().rstrip('\n')
            line = re.split('\t', line)
            q20 = float(line[3])
            line = r.readline().rstrip('\n')
            line = re.split('\t', line)
            q30 = float(line[3])
        my_q20 = max(q20, q30)
        my_q30 = min(q20, q30)
        fastq['q20'] = my_q20
        fastq['q30'] = my_q30
        fastq['size'] = self._get_size(fastq["file_path"])
        return fastq

    def get_p_stat(self, p_id):
        p_stat = dict()
        rawValidStat = os.path.join(self.option("rawValid_path"), "output", "stat.xls")
        with open(rawValidStat, 'rb') as r:
            for line in r:
                line = line.strip("\n")
                line = re.split("\t", line)
                if line[0] == p_id:
                    p_stat["Raw_pair"] = int(line[2])
                    p_stat["valid_pair"] = int(line[3])
                    p_stat["valid_rate"] = line[4]
                    p_stat["chimeric"] = int(line[5])
                    p_stat["chimeric_rate"] = "{:.2f}%".format((int(line[5]) / int(line[2])) * 100)
                    p_stat["NoRBarcode"] = int(line[6])
                    p_stat["NoFBarcode"] = int(line[7])
                    break
        trimStat = os.path.join(self.option("qualcontrol_path"), "output", "stat.xls")
        with open(trimStat, 'rb') as r:
            for line in r:
                line = line.strip("\n")
                line = re.split("\t", line)
                if line[0] == p_id:
                    p_stat["Pair_trim"] = int(line[2])
                    p_stat["Trim_rate"] = line[3]
        mergeStat = os.path.join(self.option("merge_path"), "output", "stat.xls")
        with open(mergeStat, 'rb') as r:
            for line in r:
                line = line.strip("\n")
                line = re.split("\t", line)
                if line[0] == p_id:
                    p_stat["Pair_merge"] = int(line[2])
                    p_stat["merge_rate"] = line[3]
        markSeqStat = os.path.join(self.option('markSeq_path'), "output", "stat.xls")
        with open(markSeqStat, 'rb') as r:
            for line in r:
                line = line.strip("\n")
                line = re.split("\t", line)
                if line[0] == p_id:
                    p_stat["primer_mismatch"] = int(line[4])
                    p_stat["RbarcodeMissMatch2"] = int(line[5])
                    p_stat["noBarcode2"] = int(line[6])
        lenControlStat = os.path.join(self.option('lenControl_path'), 'output', 'stat.xls')
        with open(lenControlStat, 'rb') as r:
            for line in r:
                line = line.strip("\n")
                line = re.split("\t", line)
                if line[0] == p_id:
                    p_stat["final_reads"] = int(line[3])
        p_stat["split_rate"] = "{:.2f}%".format((p_stat["final_reads"] / p_stat["Pair_merge"]) * 100)
        p_stat["highQuality_rate"] = "{:.2f}%".format((p_stat["final_reads"] / p_stat["valid_pair"]) * 100)
        (p_stat["q20"], p_stat["q30"]) = self.get_parent_q20q30(p_id)
        return p_stat

    def get_child_sample_info(self, c_id):
        """
        获取一个子样本的信息
        """
        c_info = dict()
        c_info["sample_id"] = c_id
        c_stat = dict()
        p_id = self.option("sample_info").find_parent_id(c_id)
        child_file = os.path.join(self.option("seqExtract_path"), "child_sample", "library_" + p_id,
                                  "final", c_id + ".final.fastq")
        num_lines = sum(1 for line in open(child_file))
        c_stat["total"] = int(num_lines / 4)
        fastq1 = dict()
        library_type = self.option('sample_info').parent_sample(p_id, "library_type")
        if library_type is None:
            library_type = "undefine"
        final_path = os.path.join(self.seq_id, library_type, "child_sample", "library_" + p_id,
                                  "final", c_id + ".final.fastq.gz")
        rawValid_r1_path = os.path.join(self.seq_id, library_type, "child_sample", "library_" + p_id,
                                        "rawValid", c_id + ".rawValid_r1.fastq.gz")
        rawValid_r2_path = os.path.join(self.seq_id, library_type, "child_sample", "library_" + p_id,
                                        "rawValid", c_id + ".rawValid_r2.fastq.gz")
        fastq1["final"] = final_path
        fastq1["rawValid_r1"] = rawValid_r1_path
        fastq1["rawValid_r2"] = rawValid_r2_path
        fastq1['size'] = self._get_size(fastq1["final"])
        c_info["stat"] = c_stat
        c_info["fastq_info"] = {"fastq1": fastq1}
        return c_info

    def _get_size(self, file_path):
        """
        统计一个文件的大小
        """
        size = int(os.stat(file_path).st_size)
        """
        level = {
            1: "B", 2: "KB", 3: "MB", 4: "GB", 5: "TB", 6: "PB"
        }
        i = 1
        while(i < 6 and size > 1000):
            size = size / 1000
            i += 1
        read_size = "{:.2f} {}".format(size, level[i])
        """
        return str(size)

    def run(self):
        super(SplitStatTool, self).run()
        self.create_json()
        self.end()
