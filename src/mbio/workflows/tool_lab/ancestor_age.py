# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class AncestorAgeWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AncestorAgeWorkflow, self).__init__(wsheet_object)
        options = [
            # {"name":"out_bam", "type":"infile"},#存放外群的bam文件
            # {"name":"out_bam", "type":"string"},
            {"name":"sample1","type":"string"},#用作共祖时间计算的样本1，要提供单倍型
            # {"name":"c_bam", "type":"infile"},#存放用于比较的样本的bam文件
            # {"name":"c_bam", "type":"string"},
            {"name":"csample1","type":"string"},#用作共祖时间计算的比较样本1
            {"name":"csample2","type":"string"},#用作共祖时间计算的比较样本2
            {"name":"csample3","type":"string"},#用作共祖时间计算的比较样本3
            {"name":"csample4","type":"string"},#用作共祖时间计算的比较样本4
            {"name":"csample5","type":"string"},#用作共祖时间计算的比较样本5
            {"name":"outgroup1","type":"string"},#用作外群的bam文件
            {"name":"outgroup2","type":"string"},
            {"name":"outgroup3","type":"string"},
            {"name":"outgroup4","type":"string"},
            {"name":"outgroup5","type":"string"},
            {"name":"num_c_samples","type":"int"},
            {"name":"num_outgroup","type":"int"},
            {"name":"haplogroup", "type": "string"},#样本一的单倍型
            {"name":"main_id", "type":"string"},
            {"name":"update_info", "type": "string"},
        ]

        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.bam2vcf = self.add_tool("tool_lab.bam2vcf")
        self.ancestor_age = self.add_tool("tool_lab.ancestor_age_v2")

    def check_option(self):
        # if not self.option("bam_dir"):
        #     raise OptionError("必须输入bam文件夹")
        if not self.option("sample1"):
            raise OptionError("请选择需要计算共祖时间的样品1")
        # if not self.option("sample2"):
        #     raise OptionError("请选择另一个需要计算共祖时间的样品2")
        if not self.option("haplogroup"):
            raise OptionError("请输入样品1所属的单倍型")
        if not self.option("outgroup1"):
            raise OptionError("请至少输入一个外群bam文件")
        if not self.option("csample1"):
            raise OptionError("请至少输入一个比较样本bam文件")

    # def mk_bam_string(self):
    #     self.outbam = "{}".format(self.option("outgroup1").prop["path"])
    #     for i in range(2, 6):
    #         self.outbam = self.outbam + ",{}".format(self.option("outgroup{}".format(i)).prop["path"])
    #     self.cbam = "{}".format(self.option("outgroup1").prop["path"])
    #     for i in range(2, 6):
    #         self.outbam = self.outbam + ",{}".format(self.option("outgroup{}".format(i)).prop["path"])
        
    def mk_bam_dir(self):
        '''
        把比较样本csample和外群outgroup，分别放入c_bam和out_bam文件夹
        '''
        if not os.path.exists(self.work_dir + "/temp"):
            os.mkdir(self.work_dir + "/temp")
        if not os.path.exists(self.work_dir + "/temp/out_bam"):
            os.mkdir(self.work_dir + "/temp/out_bam")
        if not os.path.exists(self.work_dir + "/temp/c_bam"):
            os.mkdir(self.work_dir + "/temp/c_bam")
        
        for i in range(1, 6):
            if self.option("outgroup{}".format(i)):
                outgroup_bam = ""
                if self.option("outgroup{}".format(i)).startswith("YGB"):
                    outgroup_bam = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/kehu/{}.hg38.chrY.dedup.sort.bam".format(self.option("outgroup{}".format(i)))
                    if not os.path.exists(outgroup_bam):
                        outgroup_bam = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/kehu/{}.hg38.chrY.dedup.bam".format(self.option("outgroup{}".format(i)))
                else:
                    outgroup_bam = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/Yoogene/{}.hg38.chrY.dedup.sort.bam".format(self.option("outgroup{}".format(i)))
                    if not os.path.exists(outgroup_bam):
                        outgroup_bam = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/Yoogene/{}.hg38.chrY.dedup.bam".format(self.option("outgroup{}".format(i)))
                try:
                    shutil.copy(outgroup_bam,os.path.join(self.work_dir, "temp/out_bam"))
                except Exception as e:
                    self.set_error("获取原始文件发生错误：{}".format(e))
        for i in range(1, 6):
            if self.option("csample{}".format(i)):
                csample = ""
                if self.option("csample{}".format(i)).startswith("YGB"):
                    csample = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/kehu/{}.hg38.chrY.dedup.sort.bam".format(self.option("csample{}".format(i)))
                    if not os.path.exists(outgroup_bam):
                        csample = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/kehu/{}.hg38.chrY.dedup.bam".format(self.option("csample{}".format(i)))
                else:
                    csample = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/Yoogene/{}.hg38.chrY.dedup.sort.bam".format(self.option("csample{}".format(i)))
                    if not os.path.exists(outgroup_bam):
                        csample = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/Yoogene/{}.hg38.chrY.dedup.bam".format(self.option("csample{}".format(i)))
                try:
                    shutil.copy(csample,os.path.join(self.work_dir, "temp/c_bam"))
                except Exception as e:
                    self.set_error("获取原始文件发生错误：{}".format(e))
        # if not os.path.exists("{}/temp/c_bam".format(self.work_dir)):
        #     self.set_error("文件{}没有传输成功".format(os.path.basename(self.option("csample1").prop["path"])))
        # else:
        #     self.logger.info("文件{}传输成功".format(os.path.basename(self.option("csample1").prop["path"])))
        # if not os.path.exists("{}/temp/c_bam".format(self.work_dir)):
        #     self.set_error("文件{}没有传输成功".format(os.path.basename(self.option("outgroup1").prop["path"])))
        # else:
        #     self.logger.info("文件{}传输成功".format(os.path.basename(self.option("outgroup1").prop["path"])))
    
    # def download_from_s3_(self):
    #     """
    #     从对象存储中下载文件到指定路径
    #     :return:
    #     """
    #     if not os.path.exists(self.work_dir + "/temp"):
    #         os.mkdir(self.work_dir + "/temp")
    #     if not os.path.exists(self.work_dir + "/temp/out_bam"):
    #         os.mkdir(self.work_dir + "/temp/out_bam")
    #     if not os.path.exists(self.work_dir + "/temp/c_bam"):
    #         os.mkdir(self.work_dir + "/temp/c_bam")
    #     self.logger.info("开始下载对象存储中的文件！")
    #     transfer = MultiFileTransfer()
    #     for sample in self.option("out_bam").split(","):
    #         transfer.add_download(sample,"{}/temp/out_bam".format(self.work_dir))
    #         if not exists(sample):
    #             self.set_error("文件%s不存在！"%sample)
    #     for sample2 in self.option("c_bam").split(","):
    #         transfer.add_download(sample,"{}/temp/c_bam".format(self.work_dir))
    #         if not exists(sample2):
    #             self.set_error("文件%s不存在！"%sample2)    
    #     transfer.perform()
    #     self.logger.info("下载对象存储中的文件成功！")

    def run_bam2vcf(self):
        self.sample1 = ""
        if self.option("sample1").startswith("YGB"):
            self.sample1 = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/kehu/{}.hg38.chrY.dedup.sort.bam".format(self.option("sample1"))
            if not os.path.exists(self.sample1):
                self.sample1 = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/kehu/{}.hg38.chrY.dedup.bam".format(self.option("sample1"))
        else:
            self.sample1 = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/Yoogene/{}.hg38.chrY.dedup.sort.bam".format(self.option("sample1"))
            if not os.path.exists(self.sample1):
                self.sample1 = "/mnt/lustre/users/sanger-dev/ysource/workspace/data/bam/Yoogene/{}.hg38.chrY.dedup.bam".format(self.option("sample1"))
        # sample1_path = os.path.join()
        self.bam2vcf.set_options({
            "og_bam_dir": os.path.join(self.work_dir, "temp/out_bam"),#存放外群的bam文件的文件夹
            "sample1": self.sample1,#用作共祖时间计算的样本1，要提供单倍型
            "c_bam_dir": os.path.join(self.work_dir, "temp/c_bam"),#存放用于比较的样本的bam文件的文件夹
            # "bam_dir": self.option("bam_dir"),
            # "sample1": self.option("sample1"),
            # "sample2": self.option("sample2"),
            # "outgroup1":self.option("outgroup1"),
            # "outgroup2":self.option("outgroup2"),
            # "outgroup3":self.option("outgroup3"),
            # "outgroup4":self.option("outgroup4"),
            # "outgroup5":self.option("outgroup5"),
        })
        self.bam2vcf.on('end', self.run_ancestorage, "bam2vcf")
        
        self.bam2vcf.run()
    
    def run_ancestorage(self):
        self.logger.info('x'*30)
        self.logger.info('{}/{}'.format(self.bam2vcf.output_dir,'samples-filter-8.4.vcf'))
        self.logger.info('{}/{}'.format(os.path.join(self.work_dir, 'Bam2vcf/'), 'bam1.list'))
        self.ancestor_age.set_options({
            "vcf": self.bam2vcf.output_dir + '/samples-filter-8.4.vcf',
            "bam_list": self.bam2vcf.output_dir + '/bam1.list',
            "haplogroup": self.option("haplogroup"),
            "sample1": self.sample1,
            "out_group": self.option("num_c_samples") + 2,
        })
        self.ancestor_age.on('end', self.set_output, "ancestor_age")
        self.ancestor_age.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'bam2vcf':
            self.linkdir(obj.output_dir, 'bam2vcf')
        if event['data'] == 'ancestor_age':
            self.linkdir(obj.output_dir, 'ancestor_age')
        self.set_db()

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

    def set_db(self):
        self.logger.info("开始导表")
        s1 = self.option('sample1')
        # s2 = self.option('sample2').split('/')[-1].split('.')[0]
        api_ancestorage = self.api.api("tool_lab.ancestor_age")
        api_ancestorage.add_ancestorage_detail(self.option('main_id'), os.path.join(self.output_dir, "ancestor_age/{}_ancestor_age.txt".format(s1)))
        self.logger.info("导表结束")
        self.end()

    def run(self):
        # self.download_from_s3_()
        self.mk_bam_dir()
        self.run_bam2vcf()
        super(AncestorAgeWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AncestorAgeWorkflow, self).end()





