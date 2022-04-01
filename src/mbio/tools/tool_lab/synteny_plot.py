#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import os
# import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import shutil
import subprocess
import re


class SyntenyPlotAgent(Agent):
    def __init__(self, parent):
        super(SyntenyPlotAgent, self).__init__(parent)
        options = [
            {"name": "result", "type": "infile", "format": "tool_lab.simple"},
            {"name": "ref", "type": "string"},
            {"name": "samples", "type": "string"},
            {"name": "seq_dir", "type": "infile", "format": "sequence.fasta_dir"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检测
        """
        if not self.option("result").is_set and not self.option("seq_dir").is_set:
            raise OptionError("请传入mumer结果和序列文件夹！")
        if not self.option("ref") and not self.option("samples"):
            raise OptionError("请传入样本名！")

    def set_resource(self):
        """
        设置所需资源
        """
        self.cpu = 4
        self.memory = '20G'

    def end(self):
        super(SyntenyPlotAgent, self).end()


class SyntenyPlotTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(SyntenyPlotTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.mcscanx = software_dir + '/bioinfo/tool_lab/MCScanX/MCScanX-master/'
        self.java_path = self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/bin/"

    def get_gff(self):
        """
        得到简化的gff
        """
        with open(self.option("seq_dir").prop["path"] + "/" + self.option("samples") + ".fna","r") as f,\
                open(self.option("seq_dir").prop["path"] + "/" + self.option("ref") + ".fna", "r") as v,\
                open(self.work_dir+ "/tmp_dir/all.gff","w") as t:
            data1 = f.read()
            num = 1
            for i in data1.split(">"):
                if i.strip():
                    seq_len = 0
                    for x in i.strip().split("\n")[1:]:
                        seq_len += len(x.strip())
                    t.write(self.option("samples") + "\t" + i.strip().split("\n")[0].split(" ")[0] + "_query" + "\t" + str(num) + "\t" + str(num + seq_len) + "\n")
                    num += (seq_len+1)

            data2 = v.read()
            num = 1
            for i in data2.split(">"):
                if i.strip():
                    seq_len = 0
                    for x in i.strip().split("\n")[1:]:
                        seq_len += len(x.strip())
                    t.write(
                        self.option("ref") + "\t" + i.strip().split("\n")[0].split(" ")[0] + "_ref" + "\t" + str(num) + "\t" + str(num + seq_len) + "\n")
                    num += seq_len

    def get_result(self):
        with open(self.option("result").prop["path"],"r") as f,open(self.work_dir+ "/tmp_dir/all.blast","w") as t:
            data = f.readlines()
            for i in data[1:]:
                tmp = i.strip().split("\t")
                if tmp[10] == self.option("samples") and tmp[9] == self.option("ref"):
                    t.write(tmp[8]+"_query"+"\t"+tmp[7]+"_ref"+"\t"+tmp[6]+"\t"+tmp[4]+"\t"+"0"+"\t"+"0"+"\t"+tmp[2]+"\t"+tmp[3]+"\t"+tmp[0]+"\t"+tmp[1]+"\t"+"0"+"\t"+"500"+"\n")

    def plot(self):
        cmd1 = '{}  {}'.format("/bioinfo/tool_lab/MCScanX/MCScanX-master/MCScanX",self.work_dir+ "/tmp_dir/all")
        command = self.add_command('mcscanx', cmd1).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('MCScanX运行成功')
        else:
            self.set_error('MCScanX运行出错')

    def plot2(self):
        with open(self.work_dir+ "/tmp_dir/control1","w") as t1:
            t1.write("800"+"\n"+"800"+"\n"+self.option("samples")+"\n"+self.option("ref")+"\n")

        with open(self.work_dir+ "/tmp_dir/control2","w") as t2:
            t2.write("800" + "\n" + self.option("samples")+","+self.option("ref") + "\n")

        os.system("dos2unix {}".format(self.work_dir+ "/tmp_dir/control1"))
        os.system("dos2unix {}".format(self.work_dir + "/tmp_dir/control2"))
        os.chdir(self.mcscanx+"downstream_analyses")
        cmd2 = '{}java -Djava.awt.headless=true dot_plotter  -g {} -s {} -c {} -o {}'.format("/program/sun_jdk1.8.0/bin/", self.work_dir+ "/tmp_dir/all.gff",self.work_dir+ "/tmp_dir/all.collinearity",
                                                                               self.work_dir+ "/tmp_dir/control1",self.work_dir+ "/tmp_dir/dot.PNG")
        command = self.add_command('dot_plotter', cmd2).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('dot_plotter运行成功')
        else:
            self.set_error('dot_plotter运行出错')

    def plot3(self):
        os.chdir(self.mcscanx + "downstream_analyses")
        cmd3 = '{}java -Djava.awt.headless=true dual_synteny_plotter -g {} -s {} -c {} -o {}'.format("/program/sun_jdk1.8.0/bin/"
                                                                      , self.work_dir + "/tmp_dir/all.gff",
                                                                      self.work_dir + "/tmp_dir/all.collinearity",
                                                                      self.work_dir + "/tmp_dir/control1",
                                                                      self.work_dir + "/tmp_dir/dual_synteny.PNG")
        command = self.add_command('dual_synteny_plotter', cmd3).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('dual_synteny_plotter运行成功')
        else:
            self.set_error('dual_synteny_plotter运行出错')

    def plot4(self):
        os.chdir(self.mcscanx + "downstream_analyses")
        cmd4 = '{}java -Djava.awt.headless=true circle_plotter  -g {} -s {} -c {} -o {}'.format("/program/sun_jdk1.8.0/bin/"
                                                                      , self.work_dir + "/tmp_dir/all.gff",
                                                                      self.work_dir + "/tmp_dir/all.collinearity",
                                                                      self.work_dir + "/tmp_dir/control2",
                                                                      self.work_dir + "/tmp_dir/circle.PNG")
        command = self.add_command('circle_plotter', cmd4).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('circle_plotter运行成功')
        else:
            self.set_error('circle_plotter运行出错')



    def set_output(self):
        """
        将结果文件链接至output
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        os.link(self.work_dir+ '/tmp_dir/dot.PNG', self.output_dir + '/dot.PNG')
        os.link(self.work_dir + '/tmp_dir/dual_synteny.PNG', self.output_dir + '/dual_synteny.PNG')
        os.link(self.work_dir + '/tmp_dir/circle.PNG', self.output_dir + '/circle.PNG')
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(SyntenyPlotTool, self).run()
        if os.path.exists(self.work_dir + "/tmp_dir"):
            shutil.rmtree(self.work_dir + "/tmp_dir")
        os.mkdir(self.work_dir + "/tmp_dir")
        self.get_gff()
        self.get_result()
        self.plot()
        self.plot2()
        self.plot3()
        self.plot4()
        self.set_output()
        self.end()