# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'WangWenjie'

import os
import time
import sys
import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd


class DeCarrierAgent(Agent):
    """
    去载体,用于U1直接去载体或U2拼接后去载体
    """
    def __init__(self, parent):
        super(DeCarrierAgent,self).__init__(parent)
        options = [
            {"name": "input_path", "type": "infile", "format": "denovo_rna_v2.common_dir"}, # 输入去载体结果文件夹
            {"name": "operation_type", "type": "string", "default": "U1"}  # 去载体类型
        ]
        self.add_option(options)       

    def check_options(self):
        """
        参数检查
        """
        if not self.option("input_path"):
            raise OptionError("没有找到input_path")
        if self.option("operation_type") not in ["U1", "U2"]:
            raise OptionError("operation_type只能是U1/U2")
    
    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 8
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(DeCarrierAgent, self).end()
       
class DeCarrierTool(Tool):  
    def __init__(self, config):
        super(DeCarrierTool, self).__init__(config)
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.software_dir = self.config.SOFTWARE_DIR
        self.python = "miniconda2/bin/python"
        self.phred = 'bioinfo/de_carrier/genome/bin/phred'
        # self.phredpar_dat = self.config.SOFTWARE_DIR+'bioinfo/de_carrier/genome/lib/phredpar.dat'
        self.set_environ(PHRED_PARAMETER_FILE=self.config.SOFTWARE_DIR + '/bioinfo/de_carrier/genome/lib/phredpar.dat')
        self.cross_match = 'bioinfo/de_carrier/genome/bin/cross_match'
        self.pairassembly = self.config.SOFTWARE_DIR+'/bioinfo/de_carrier/genome/bin/pairassembly'
        self.contig_crossmatch = self.config.SOFTWARE_DIR+'/bioinfo/de_carrier/genome/bin/contig-crossmatch'
        self.sample_seq = self.config.SOFTWARE_DIR+'/bioinfo/de_carrier/genome/lib/pMD18-T.seq'

    def run_U1(self):
        """
        U1:直接去载体
        """
        cmd = "{} -id  {} -sa 1Y.fa -qa raw.fa.screen.qual -trim_fasta -trim_alt '' -trim_cutoff 0.01".format(self.phred, self.option("input_path").prop["path"])
        out_sa = os.path.join(self.work_dir, "1Y.fa")
        out_qa = os.path.join(self.work_dir, "raw.fa.screen.qual")
        self.logger.info(cmd)
        command = self.add_command("phred", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行phred完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        cmd1 = "{} {} {} -minmatch 12 -penalty -2 -minscore 20 -screen".format(self.cross_match, out_sa, self.sample_seq)
        out_u1 = os.path.join(self.work_dir, "1Y.fa.screen")
        self.logger.info(cmd1)
        command1 = self.add_command("cross_match", cmd1).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行cross_match完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")

    def run_U2(self):
        """
        U2:拼接后去载体
        """
        cmd2 = os.system('cd {} && {}'.format(self.option("input_path").prop["path"], self.pairassembly))
        self.wait()
        if cmd2 == 0:
            self.logger.info("运行pairassembly成功，拼接完成")
        else:
            self.logger.info("运行pairassembly失败")
        cmd3 = os.system('cd {} && cp -r ace {}'.format(self.option("input_path").prop["path"], self.work_dir))
        self.wait()
        if cmd3 == 0:
            self.logger.info("ace文件生成成功")
        else:
            self.logger.info("ce文件生成失败")
        cmd4 = os.system('cd {} && cp -r contigs {}'.format(self.option("input_path").prop["path"], self.work_dir))
        self.wait()
        if cmd4 == 0:
            self.logger.info("contigs文件生成成功")
        else:
            self.logger.info("contigs文件生成失败")
        cmd5 = os.system('cd {} && {} {}'.format(self.work_dir+"/contigs", self.contig_crossmatch, self.sample_seq))
        self.wait()
        if cmd5 == 0:
            self.logger.info("拼接后contigs_fasta文件生成完成")
        else:
            self.logger.info("拼接后contigs_fasta文件生成失败")
        out_contigs_fasta = os.path.join(self.work_dir+"/contigs", "all.Contig.fasta")
        cmd6 = "{} {} {} -minmatch 12 -penalty -2 -minscore 20 -screen".format(self.cross_match, out_contigs_fasta, self.sample_seq)
        self.logger.info(cmd6)
        command6 = self.add_command("crossmatch", cmd6).run()
        self.wait(command6)
        if command6.return_code == 0:
            self.logger.info("拼接后去载体完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")

    def set_output(self):
        # localtime = time.localtime(time.time()) # 获取当前时间
        # time = time.strftime('%Y%m%d',time.localtime(time.time())) # 把获取的时间转换成"年月日格式”
        if(self.option('operation_type') == "U1"): 
            os.link(os.path.join(self.work_dir, "1Y.fa.screen"),
                        os.path.join(self.output_dir, "1Y.fa.screen"))
        elif(self.option('operation_type') == "U2"):
            os.link(os.path.join(self.work_dir+"/contigs", "all.Contig.fasta.screen"),
                        os.path.join(self.output_dir, "all.Contig.fasta.screen"))

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        if self.option("operation_type") == "U1":
            self.run_U1()  
        else:
            self.run_U2()
        self.set_output()
        self.end()
        super(DeCarrierTool, self).run()