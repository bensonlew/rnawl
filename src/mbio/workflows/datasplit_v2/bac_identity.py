# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
import unittest
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class BacIdentityWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BacIdentityWorkflow, self).__init__(wsheet_object)
        options = [
            # {"name":"list_xls","type":"infile","format": "denovo_rna_v2.common"},
            # {"name": "sample_info","type":"infile","format":"denovo_rna_v2.common"},
            # {"name": "list_url", "type": "string"},
            # {"name": "sample_info_url", "type": "string"},
            {"name":"main_id", "type":"string"},
            {"name":"update_info", "type": "string"},
            {"name": "operation_type", "type": "string", "default": "merge"},
            {"name": "fasta_dir", "type": "infile", "format": "denovo_rna_v2.common_dir"},
            {"name": "list_url", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "sample_info_url","type":"infile","format":"denovo_rna_v2.common"},
            {"name": "new_pictrue_name","type":"infile","format":"denovo_rna_v2.common_dir"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.bac_identools = self.add_tool("tool_lab.bac_identitytonb")
        # self.bac_datapre = self.add_tool("tool_lab.bac_datapre")
        self.zip_tool = self.add_tool("tool_lab.zipresult")

    # def check_option(self):
        # if not self.option("list_url"):
        #     raise OptionError("请设置list_url")
        # if not self.option("sample_info_url"):
        #     raise OptionError("请设置sample_info_url")
    def check_option(self):
        if not self.option("fasta_dir"):
            raise OptionError("请设置fasta_dir")
        if not self.option("list_url"):
            raise OptionError("请设置list_url")
        if not self.option("sample_info_url"):
            raise OptionError("请设置sample_info_url")
        if not self.option("new_pictrue_name"):
            raise OptionError("没有找到new_pictrue_name")
# def mk_fasta_dir(self):
    #     """
    #     从list表里得到样本编号，然后从样本编号里找到需要的seq文件和ab1文件
    #     """
    #     self.list_file = self.option("lsit_xls").prop["path"]
    #     fn_list = []
    #     if os.path.exists(os.path.join(self.work_dir, "fasta_dir")):
    #         shutil.rmtree(os.path.join(self.work_dir, "fasta_dir"))
    #     os.mkdir(os.path.join(self.work_dir, "fasta_dir"))
    #     with open(self.list_file,"r") as lf:
    #         while 1:
    #             line = lf.readline()
    #             if not line:
    #                 break
    #             file_name = line.split("\t")[1]
    #             if file_name not in fn_list:
    #                 fn_list.append(file_name)
    #     for root,fdir,file in os.walk(self.all_fasta_dir):
    #         if file.split("-")[0] in fn_list:
    #             try:
    #                 shutil.copy(os.path.join(root,file),os.path.join(self.work_dir, "fasta_dir"))
    #             except Exception as e:
    #                 self.logger.info("输入原始数据失败{}".format(e))

    # def run_bac_datapre(self):
    #     """
    #     下载list文件
    #     检查原始数据命名和格式
    #     """
    #     self.bac_datapre.set_options({
    #     self.bac_datapre.set_options({
    #         "list_url":self.option("list_url"),
    #         "sample_info_url":self.option("sample_info_url")
    #     })
    #     self.bac_datapre = self.add_tool("tool_lab.bac_datapre")
    #     self.bac_datapre.set_options(options)
    #     self.bac_datapre.on("end",self.set_output,"datapre")
    #     self.bac_datapre.run()  

    def run_bac_indent(self):
        """
        输入参数运行bac_indnet小工具
        """
        self.bac_identools.set_options({
            "fasta_dir": self.option("fasta_dir"),
            "list_xls": self.option("list_url"),
            "sample_info": self.option("sample_info_url"),
            "picture_path": self.option("new_pictrue_name")
            })
        self.bac_identools.on("end",self.run_zip_result,"bac_identity")
        self.bac_identools.run()

    def run_zip_result(self):
        """
        压缩结果文件
        """
        # with open(os.path.join(self.bac_datapre.output_dir, "sample.detail.xls"),"r") as sinfo:
        with open(os.path.join(self.option("sample_info_url").prop["path"]),"r") as sinfo:
            sinfo.readline()
            line = sinfo.readline()
            temp= line.split("\t")
            temp= line.split("\t")
            client = temp[5].replace(' ','_')
            self.dir_name = "{}-{}-{}-{}-{}sample".format(temp[0],temp[1],client,temp[4],temp[6])
            # self.result_name = "{}-{}-{}-{}-{}sample".format(temp[0],temp[1],client,temp[4],temp[6])
            self.result_name = "result"
        self.zip_tool.set_options({
            "result_file": os.path.join(self.bac_identools.output_dir,self.result_name),
            "output_name": self.result_name
        })
        self.zip_tool.on("end",self.set_output,"zipresult")
        self.zip_tool.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'zipresult':
            self.linkdir(obj.output_dir, 'zipresult')
        s3_upload_dir = os.path.join(self._sheet.output, 'zipresult/{}.zip'.format(self.result_name))
        self.set_db(s3_upload_dir)

    def linkdir(self, dirpath, dirname):
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
    
    def set_db(self,s3_upload_dir):
        self.logger.info("开始导表")
        # with open(os.path.join(self.bac_datapre.output_dir, "sample.detail.xls"),"r") as sinfo:
        #     sinfo.readline()
        #     line = sinfo.readline()
        #     temp= line.split("\t")
        #     dir_name = "{}-{}-{}-{}-{}sample".format(temp[0],temp[1],temp[5],temp[4],temp[6])
        #     result_name = "{}-{}-{}-{}".format(temp[0],temp[1],temp[5],temp[4])
        ot_dir = os.path.join(self.bac_identools.output_dir, "{}".format(self.result_name))
        main_id = self.option("main_id")
        # s1 = os.path.basename(self.option('sample1').prop["path"]).split('.')[0]
        # s2 = self.option('sample2').split('/')[-1].split('.')[0]
        api_bac_ind = self.api.api("datasplit.bac_identity")
        # api_ancestorage.add_ancestorage_detail(self.option('main_id'), os.path.join(self.output_dir, "ancestor_age/{}_ancestor_age.txt".format(s1)))
        api_bac_ind.update_bac_db_record(main_id, s3_upload_dir,self.dir_name)
        api_bac_ind.add_qc_report(main_id, os.path.join(ot_dir, "sanger.report"))
        api_bac_ind.add_blastnt_tax(main_id , os.path.join(ot_dir, "{}/NCBI/blastM8.xls".format(self.result_name)))
        self.logger.info("导表结束")
        self.end()

    def run(self):
        # self.mk_fasta_dir()
        # self.run_bac_datapre()
        self.run_bac_indent()
        super(BacIdentityWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BacIdentityWorkflow, self).end()

# class TestFunction(unittest.TestCase):
#     """
#     测试脚本用
#     """

#     def test(self):
#         import random
#         from mbio.workflows.single import SingleWorkflow
#         from biocluster.wsheet import Sheet
#         data = {
#             'id': 'bac_datapre_' + str(random.randint(1, 10000)),
#             "type": "workflow",
#             "name": "tool_lab.bac_identity",
#             "options": {
#                  "list_url":"http://ngs.tngs.com/xlsx/list.xls",
#                  "sample_info_url":"http://ngs.tngs.com/xlsx/sample.detail.xls"
#             }
#         }
#         wsheet = Sheet(data=data)
#         wf = SingleWorkflow(wsheet)
#         wf.run()


# if __name__ == '__main__':
#     unittest.main()

if __name__ == '__main__':
    import random
    from biocluster.wsheet import Sheet
    data = {
            'id': 'bac_datapre_' + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.bac_identity",
            "options": {
                 "main_id":"5e9e6a6017b2bf2049a81be3",
                 "list_url":"http://ngs.tngs.com/xlsx/list.xls",
                 "sample_info_url":"http://ngs.tngs.com/xlsx/sample.detail.xls"
            }
        }
    
    wsheet =Sheet(data=data)
    wf = BacIdentityWorkflow(wsheet)
    wf.run()
