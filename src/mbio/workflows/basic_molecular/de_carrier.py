# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'WangWenjie'

import os
import re
import time
import shutil
import json
import unittest
import glob
import datetime
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class DeCarrierWorkflow(Workflow):
    """
    去载体流程：
    U1：序列小于700，直接去载体
    U2：序列大于等于700，拼接后去载体
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DeCarrierWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"main_id", "type":"string"},
            # {"name": "input_path", "type": "infile", "format": "denovo_rna_v2.common_dir"}, 
            # {"name": "operation_type", "type": "string", "default": "U1"},  # 去载体类型
            {"name":"list_file","type":"infile","format":"denovo_rna_v2.common"},
            # {"name":"u2_list","type":"infile","format":"denovo_rna_v2.common"},
            {"name":"update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.de_carrier = self.add_tool("basic_molecular.de_carrier")

    def check_option(self):
        """
        参数检查
        """
        if not self.option("list_file"):
            raise OptionError("没有任何样品输入")
        # if self.option("operation_type") not in ["U1", "U2"]:
        #     raise OptionError("operation_type只能是U1/U2")
            
    def run_decarrier(self):
        """
        输入参数运行de_carrier小工具
        """
        self.decarrier_tools = []
        if self.u1num > 0 :
            self.de_carrier.set_options({
                "input_path":self.u1dir,
                "operation_type":"U1"
                })
            # self.de_carrier.on("end", self.set_output, "de_carrier")
            self.decarrier_tools.append(self.de_carrier)
            # self.de_carrier.run()
        if self.u2num > 0 :
            self.de_carrier.set_options({
                "input_path":self.u2dir,
                "operation_type":"U2"
                })
            # self.de_carrier.on("end", self.set_output, "de_carrier")
            self.decarrier_tools.append(self.de_carrier)
        if len(self.decarrier_tools) == 1:
            self.decarrier_tools[0].on('end',self.set_output)
        else:
            self.on_rely(self.decarrier_tools,self.set_output)
        for tool in self.decarrier_tools:
            tool.run()

    def set_output(self):
        # obj = event['bind_object']
        # if event['data'] == 'de_carrier':
        #     self.linkdir(obj.output_dir, 'de_carrier')
        #     self.end()
        now_time = datetime.datetime.now().strftime("%Y%m%d")
        self.output_name = ""
        if self.u1num + self.u2num == 1:
            if self.u1num == 1:
                target_file = os.path.join(self.decarrier_tools[0].output_dir,"1Y.fa.screen")
                link_file = os.path.join(self.output_dir,"{}-{}--{}.screen".format(self.info,self.sn,now_time))
                self.output_name ="{}-{}--{}.screen".format(self.info,self.sn,now_time)
                if os.path.exists(link_file):
                    os.remove(link_file)
                os.link(target_file,link_file)
            elif self.u2num == 1:
                target_file = os.path.join(self.decarrier_tools[0].output_dir,"all.Contig.fasta.screen")
                link_file = os.path.join(self.output_dir,"{}-{}--{}.screen".format(self.info,self.sn,now_time))
                self.output_name ="{}-{}--{}.screen".format(self.info,self.sn,now_time)
                if os.path.exists(link_file):
                    os.remove(link_file)
                os.link(target_file,link_file)
            self.logger.info('finish set_output at {}'.format(
                self.__class__.__name__))
        else:
            self.output_name = "{}-all--{}.screen".format(self.info,now_time)
            with open(os.path.join(self.output_dir, "{}-all--{}.screen".format(self.info,now_time)), 'w') as rs:
                for tool in self.decarrier_tools:
                    for source in glob.glob(os.path.join(tool.output_dir, '*')):
                        self.logger.info("{}".format(source))
                        with open(source,"r") as sor:
                            while 1:
                                line = sor.readline()
                                if not line:
                                    break
                                rs.write(line)
            self.logger.info('finish set_output at {}'.format(
                self.__class__.__name__))
        self.set_db()

    # def linkdir(self, dirpath, dirname):
    #     """
    #     link一个文件夹下的所有文件到本的output目录
    #     """
    #     allfiles = os.listdir(dirpath)
    #     newdir = os.path.join(self.output_dir, dirname)
    #     if not os.path.exists(newdir):
    #         os.mkdir(newdir)
    #     oldfiles = [os.path.join(dirpath, i) for i in allfiles]
    #     newfiles = [os.path.join(newdir, i) for i in allfiles]
    #     for newfile in newfiles:
    #         if os.path.exists(newfile):
    #             if os.path.isfile(newfile):
    #                 os.remove(newfile)
    #             else:
    #                 os.system('rm -r %s' % newfile)
    #                 # self.logger.info('rm -r %s' % newfile)
    #     for i in range(len(allfiles)):
    #         if os.path.isfile(oldfiles[i]):
    #             os.link(oldfiles[i], newfiles[i])
    #         elif os.path.isdir(oldfiles[i]):
    #             # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
    #             os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        """
        保存结果到mongo数据库
        """
        self.logger.info("保存结果到mongo数据库")
        s3_upload_path = os.path.join(self._sheet.output, '{}'.format(self.output_name))
        # ot_dir = os.path.join(self.bac_identools.output_dir, "{}".format(self.result_name))
        main_id = self.option("main_id")
        de_carrier_api = self.api.api("basic_molecular.de_carrier")
        de_carrier_api.update_decarrier_record(self.option("main_id"), s3_upload_path)
        self.end()

    def get_list(self):
        list_file = self.option("list_file").prop["path"]
        self.tools = []
        self.u1num =0
        self.u2num =0
        self.u1dir = os.path.join(self.work_dir,"u1")
        self.u2dir = os.path.join(self.work_dir,"u2")
        os.mkdir(self.u1dir)
        os.mkdir(self.u2dir)
        with open(list_file,"r") as lf:
            while 1:
                line = lf.readline()
                options = {}
                if not line:
                    break
                fd = line.rstrip().split("\t")
                self.sn, anlysis_type, ab1_s3_path,seq_s3_path,self.info = fd[0],fd[1],fd[2],fd[3],fd[4]
                if anlysis_type == "U1":
                    self.u1num += 1
                    options = {
                        "input_file": ab1_s3_path,
                        "out_path": self.u1dir

                    }
                elif anlysis_type == "U2":
                    self.u2num += 1
                    options = {
                        "input_file": ab1_s3_path,
                        "out_path": self.u2dir

                    }
                file_download = self.add_tool("basic_molecular.download_file")
                file_download.set_options(options)
                self.tools.append(file_download)
                if anlysis_type == "U1":
                    # self.u1num += 1
                    options = {
                        "input_file": seq_s3_path,
                        "out_path": self.u1dir

                    }
                elif anlysis_type == "U2":
                    # self.u2num += 1
                    options = {
                        "input_file": seq_s3_path,
                        "out_path": self.u2dir

                    }
                file_download = self.add_tool("basic_molecular.download_file")
                file_download.set_options(options)
                self.tools.append(file_download)
        if len(self.tools) == 1:
            self.tools[0].on('end',self.run_decarrier)
        else:
            self.on_rely(self.tools,self.run_decarrier)
        for tool in self.tools:
            tool.run()

    def run(self):
        self.get_list()
        super(DeCarrierWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(DeCarrierWorkflow, self).end()