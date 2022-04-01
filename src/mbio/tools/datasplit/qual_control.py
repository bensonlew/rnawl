# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""用trimmomatic对raw valid reads进行质量过滤"""
import re
import os
import errno
import subprocess
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.config import Config


class QualControlAgent(Agent):
    def __init__(self, parent=None):
        super(QualControlAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'},  # 样本拆分信息表
            {'name': 'var2end_path', 'type': "string"}  # 对一次拆分结果滤去嵌合和barcode确实后的reads
        ]
        self.add_option(options)

    def check_option(self):
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        if not self.option('var2end_path'):
            raise OptionError("参数var2end_path不能为空")
        return True

    def set_resource(self):
        self._cpu = 32
        self._memory = ''


class QualControlTool(Tool):
    def __init__(self, config):
        super(QualControlTool, self).__init__(config)
        self._version = 1.0
        self.option('sample_info').get_info()
        self.javaPath = os.path.join(Config().SOFTWARE_DIR, "program/sun_jdk1.8.0/bin/java")
        self.trimmomaticPath = os.path.join(Config().SOFTWARE_DIR, "bioinfo/seq/trimmomatic-0.36/trimmomatic-0.36.jar")

    def run_trim(self):
        cmd_list = list()
        trimDir = os.path.join(self.work_dir, "trimmo")
        try:
            os.makedirs(trimDir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(trimDir):
                pass
            else:
                raise OSError("创建目录失败")
        shLogPath = os.path.join(self.work_dir, "sh.log")
        with open(shLogPath, "ab") as a:
            a.write("trimmomatic".center(79, "#") + "\n")
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                leading = str(p["config"]["leading"])
                trailing = str(p["config"]["trailing"])
                slidingWindow = "{}:{}".format(p["config"]['slidingWindow_windowSize'], p["config"]["slidingWindow_requiredQuality"])
                minLen = str(p["config"]["minLen"])
                file_log = os.path.join(self.work_dir, "output", p['sample_id'] + ".trimmo.log")
                file_r1 = os.path.join(self.option('var2end_path'), p['sample_id'] + ".valid_r1.fastq")
                file_r2 = os.path.join(self.option('var2end_path'), p['sample_id'] + ".valid_r2.fastq")
                file_trim1 = os.path.join(trimDir, p['sample_id'] + ".trimmo_r1.fastq")
                file_unpair1 = os.path.join(trimDir, p['sample_id'] + ".trimmo_unpair_r1.fastq")
                file_trim2 = os.path.join(trimDir, p['sample_id'] + ".trimmo_r2.fastq")
                file_unpair2 = os.path.join(trimDir, p['sample_id'] + ".trimmo_unpair_r2.fastq")
                trimStr = (self.javaPath + " -jar " + self.trimmomaticPath + " PE -phred33 " + file_r1
                           + " " + file_r2 + " " + file_trim1 + " " + file_unpair1 + " " + file_trim2
                           + " " + file_unpair2 + " " + "LEADING:" + leading + " TRAILING:" + trailing
                           + " SLIDINGWINDOW:" + slidingWindow + " MINLEN:" + minLen + " 2>" + file_log)
                with open(shLogPath, "ab") as a:
                    a.write(trimStr + "\n")
                command = subprocess.Popen(trimStr, shell=True)
                cmd_list.append(command)
        for mycmd in cmd_list:
            self.logger.info("开始运行trimmomatic")
            mycmd.communicate()
        for mycmd in cmd_list:
            if mycmd.returncode == 0:
                self.logger.info("trimmomatic运行完成")
            else:
                self.set_error("trimmomatic运行出错")
                raise Exception("trimmomatic运行出错")

    def trim_stat(self):
        statPath = os.path.join(self.work_dir, "output", "stat.xls")
        with open(statPath, "wb") as w:
            w.write("#library_id\tlibrary_name\ttrimSurvived\tTrimSurvivedRate\n")
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                file_log = os.path.join(self.work_dir, "output", p['sample_id'] + ".trimmo.log")
                with open(file_log, 'rb') as r, open(statPath, 'ab') as a:
                    for line in r:
                        if re.search('Both\sSurviving:', line):
                            trimSurvived = re.search('Both\sSurviving:\s*(\d+)\s*\(.+\)\s*Forward\s*Only', line).group(1)
                            SurvivedRate = re.search('Both\sSurviving:\s*(\d+)\s*\((.+)\)\s*Forward\s*Only', line).group(2)
                            str_ = (p['sample_id'] + "\t" + p["library_name"] + "\t" + str(trimSurvived)
                                    + "\t" + str(SurvivedRate) + "\n")
                            a.write(str_)

    def run(self):
        super(QualControlTool, self).run()
        self.run_trim()
        self.trim_stat()
        self.end()
