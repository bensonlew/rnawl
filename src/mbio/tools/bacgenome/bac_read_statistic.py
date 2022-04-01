#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import subprocess
import os,re



class BacReadStatisticAgent(Agent):
    """
    用于做微生物基因组qc结果统计
    version 1.0
    author: 高豪
    last_modify: 2018.03.22
    """

    def __init__(self, parent):
        super(BacReadStatisticAgent, self).__init__(parent)
        options = [
            {"name": "sample_name", "type": "string"},  # 样品名称
            {"name": "raw_list", "type": "infile", "format": "bacgenome.list_file"},  #rawdata list
            {"name": "clean_list", "type": "infile", "format": "bacgenome.list_file"},  #cleandata list
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('sample_name'):
            raise OptionError("请设置样品名参数!", code="31400501")
        if not self.option('raw_list').is_set:
            raise OptionError("请提供rawdata list_file参数", code="31400502")
        if not self.option('clean_list').is_set:
            raise OptionError("请提供cleandata list_file参数", code="31400503")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '25G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./qc_stat.xls", "xls", "质控信息统计表"]
        ])
        super(BacReadStatisticAgent, self).end()


class BacReadStatisticTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(BacReadStatisticTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.shell_path = "/program/sh"
        self.stat_path = self.config.PACKAGE_DIR + "/bacgenome/fastq_statistics.pl"
        self.java_path = self.config.SOFTWARE_DIR +  "/program/sun_jdk1.8.0/bin/java"
        self.FastqStat_path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/QC/FastqTotalHighQualityBase.jar"
        self.bash_path = self.config.PACKAGE_DIR + "/bacgenome/q20.sh"
        self.bash2_path = self.config.PACKAGE_DIR + "/bacgenome/q30.sh"
        self.sample_path = defaultdict(list)

    def run_q20(self):
        for sample in self.sample_path:
            self.logger.info(self.sample_path[sample])
            stat_file = self.work_dir + '/' + sample + '.' + "q20.stats"
            cmd = "{} {} {} {} {} {} {} {} {} {}".format(self.shell_path, self.bash_path, self.java_path, self.FastqStat_path,
                                                      self.sample_path[sample][0], self.sample_path[sample][1],
                                                      self.sample_path[sample][2], self.sample_path[sample][3],
                                                      self.sample_path[sample][4], stat_file)
            self.logger.info(cmd)
            self.logger.info("开始运行q20统计")
            path = sample.lower() + "run_q20"
            command = self.add_command(path, cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行%s run_q20完成" %sample)
            else:
                self.set_error("运行%s run_q20运行出错!" , variables=(sample), code="31400501")
                return False


    def run_q30(self):
        for sample in self.sample_path:
            self.logger.info(self.sample_path[sample])
            stat_file = self.work_dir + '/' + sample + '.' + "q30.stats"
            cmd = "{} {} {} {} {} {} {} {} {} {}".format(self.shell_path, self.bash2_path, self.java_path, self.FastqStat_path,
                                                      self.sample_path[sample][0], self.sample_path[sample][1],
                                                      self.sample_path[sample][2], self.sample_path[sample][3],
                                                      self.sample_path[sample][4], stat_file)
            self.logger.info(cmd)
            self.logger.info("开始运行q30统计")
            path =sample.lower() + "run_q30"
            command = self.add_command(path, cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行%s run_q30完成" %sample)
            else:
                self.set_error("运行%s run_q30运行出错!" , variables=(sample), code="31400502")
                return False

    def fastq_stat(self):
        self.logger.info(self.sample_path)
        for sample in self.sample_path:
            m =re.search(r'([0-9]+)',sample)
            insert =m.group()
            self.logger.info(insert)
            q20 =self.work_dir + '/' + sample + '.' + "q20.stats"
            q30 = self.work_dir + '/' + sample + '.' + "q30.stats"
            stat_file = os.path.join(self.work_dir, "qc_stat.xls")
            cmd = "{} {} {} {} {} {} {} {} {} {}".format(self.perl_path, self.stat_path, self.option("sample_name"),insert, q20,q30,self.sample_path[sample][0],self.sample_path[sample][2],
    self.sample_path[sample][4], stat_file)
            self.logger.info(cmd)
            self.logger.info("开始运行汇总统计")
            path = sample.lower() + "summary"
            command = self.add_command(path, cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行%s完成" % path)
            else:
                self.set_error("运行%s 运行出错!" , variables=( path), code="31400503")
                return False


    def set_output(self):
        self.logger.info("set output")
        if os.path.exists(self.output_dir+"/qc_stat.xls"):
            os.remove(self.output_dir+"/qc_stat.xls")
        os.link(self.work_dir+"/qc_stat.xls", self.output_dir+"/qc_stat.xls")
        self.logger.info("done")
        self.end()

    def run(self):
        """
        运行
        """
        super(BacReadStatisticTool, self).run()
        self.get_info()
        self.run_q20()
        self.run_q30()
        self.fastq_stat()
        self.set_output()

    def get_info(self):
        base_dir = os.path.dirname(self.option("raw_list").prop['path'])
        base2_dir = os.path.dirname(self.option("clean_list").prop['path'])
        with open(self.option('raw_list').prop['path'])as fr:
            for line in fr:
                tmp = line.strip().split('\t')
                if tmp[1] in self.sample_path.keys():
                    if tmp[2] == 'l':
                        self.sample_path[tmp[1]].insert(0, base_dir + '/' + tmp[0])
                    else:
                        self.sample_path[tmp[1]].append(base_dir + '/' + tmp[0])
                else:
                    self.sample_path[tmp[1]].append(base_dir + '/' + tmp[0])
        with open(self.option('clean_list').prop['path'])as f:
            tmp_dic = {}
            for line in f:
                tmp = line.strip().split('\t')
                if tmp[1] not in tmp_dic.keys():
                    tmp_dic[tmp[1]] = {tmp[2]: tmp[0]}
                else:
                    tmp_dic[tmp[1]][tmp[2]]=tmp[0]
            for s in tmp_dic.keys():
                for k in ['l','r','s']:
                    if  k in tmp_dic[s]:
                        self.sample_path[s].append(base2_dir + '/' + tmp_dic[s][k])
                    else:
                        if k == 's':
                            os.system('touch null.unpaired.fq')
                            self.sample_path[s].append('null.unpaired.fq')




