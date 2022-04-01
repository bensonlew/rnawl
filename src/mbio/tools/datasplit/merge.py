# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""使用pear对质控之后的数据进行merge"""
import os
import re
import subprocess
import errno
import shutil
from biocluster.config import Config
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class MergeAgent(Agent):
    def __init__(self, parent=None):
        super(MergeAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'},  # 样本拆分信息表
            {'name': 'chomp_path', 'type': "string"}  # 经过质控之后序列的目录
        ]
        self.add_option(options)

    def check_option(self):
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        if not self.option('chomp_path'):
            raise OptionError("参数chomp_path不能为空")
        return True

    def set_resource(self):
        self._cpu = 32
        self._memory = ''


class MergeTool(Tool):
    def __init__(self, config):
        super(MergeTool, self).__init__(config)
        self._version = 1.0
        self.option('sample_info').get_info()
        self.pear_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/seq/pear-0.9.8/bin/pear")

    def make_ess_dir(self):
        merge_dir = os.path.join(self.work_dir, "merge")
        dir_list = [merge_dir]
        for name in dir_list:
            try:
                os.makedirs(name)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(name):
                    pass
                else:
                    raise OSError("创建目录失败")

    def pear(self):
        """
        用软件pear对质量控制后的fastq进行merge
        """
        merge_dir = os.path.join(self.work_dir, "merge")
        shLogPath = os.path.join(self.work_dir, "sh.log")
        with open(shLogPath, "ab") as a:
            a.write("PEAR".center(79, "#") + "\n")
        cmd_list = list()
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                file_r1 = os.path.join(self.option('chomp_path'), p['sample_id'] + ".chomped_r1.fastq")
                file_r2 = os.path.join(self.option('chomp_path'), p['sample_id'] + ".chomped_r2.fastq")
                pearstr = (self.pear_path + "  -p 1.0 -j 16 -f " + file_r1 + " -r " + file_r2 + " -o " +
                           merge_dir + "/pear_" + p['sample_id'] + "> " +
                           merge_dir + "/" + p['sample_id'] + ".pear.log")
                with open(shLogPath, "ab") as a:
                    a.write(pearstr + "\n")
                command = subprocess.Popen(pearstr, shell=True)
                cmd_list.append(command)
        for mycmd in cmd_list:
            self.logger.info("开始运行pear")
            mycmd.communicate()
        for mycmd in cmd_list:
            if mycmd.returncode == 0:
                self.logger.info("pear运行完成")
            else:
                self.set_error("pear运行出错")
                raise Exception("pear运行出错")

    def pear_stat(self):
        """
        对pear输出的log文件进行统计,输出到sta.xlst当中
        """
        merge_dir = os.path.join(self.work_dir, "merge")
        stat_dir = os.path.join(self.work_dir, "output")
        stat_name = os.path.join(stat_dir, "stat.xls")
        with open(stat_name, 'a') as a:
            a.write("#library_id\tlibrary_name\tpair_merged\tmerge_rate\n")
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                file_name = os.path.join(merge_dir, p['sample_id'] + ".pear.log")
                with open(file_name, 'r') as r, open(stat_name, 'a') as a:
                    for line in r:
                        if re.search(r'Assembled\sreads\s\.\.', line):
                            total_reads = re.search(r'/\s(.+)\s\(', line).group(1)
                            total_reads = re.sub(r',', '', total_reads)
                            rate = re.search(r'\((.+)\)', line).group(1)
                            a.write(p['sample_id'] + "\t" + p["library_name"] + "\t" + total_reads + "\t" + rate + "\n")
            else:
                with open(stat_name, 'a') as a:
                    a.write(p['sample_id'] + "\t" + p["library_name"] + "\t" + "NA" + "\t" + "NA" + "\n")

    def rename(self):
        merge_dir = os.path.join(self.work_dir, "merge")
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                p_id = p["sample_id"]
                file_ = os.path.join(merge_dir, "pear_" + p_id + ".assembled.fastq")
                newFile = os.path.join(merge_dir, p_id + ".merged.fastq")
                shutil.move(file_, newFile)

    def run(self):
        super(MergeTool, self).run()
        self.make_ess_dir()
        self.pear()
        self.pear_stat()
        self.rename()
        self.end()
