# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180425

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os


class AssemblyModule(Module):
    """
    用于对数据组装中捆绑bam2fastq, make_config, ngsqc_stat, soap_denovo
    """
    def __init__(self, work_id):
        super(AssemblyModule, self).__init__(work_id)
        options = [
            {"name": "bam_file", "type": "string"},  # bam文件
            {"name": "sample_id", "type": "string"},  # 页面传进来的样本id
            {"name": "pos", "type": "string"},  # 页面中的选择区域 chr1,chr1::100,000,chr1:1:,chr1:10000:20000
            {"name": "unmapping", "type": "string", "default": "true"}  # 页面选择添加的时候为true，否则为false
        ]
        self.add_option(options)
        self.bam2fastq = self.add_tool("wgs.bam2fastq")
        self.make_config = self.add_tool("wgs.make_config")
        self.assembly_stat = self.add_tool("wgs.assembly_stat")
        self.ngsqc_stat_tools = []
        self.soap_denovo_tools = []
        # self.anno = self.add_tool("wgs.anno_analysis")  # 用于后面注释

    def check_options(self):
        if not self.option("bam_file"):
            raise OptionError("缺少bam_file参数", code="24500201")
        if not self.option("sample_id"):
            raise OptionError("缺少sample_id参数", code="24500202")
        if not self.option("pos"):
            raise OptionError("pos", code="24500203")
        if not self.option("unmapping"):
            raise OptionError("缺少unmapping参数", code="24500204")
        return True

    def bam2fastq_run(self):
        self.bam2fastq.set_options({
            "bam_file": self.option("bam_file"),
            "sample_id": self.option("sample_id"),
            "pos": self.option("pos"),
            "unmapping": self.option("unmapping")
        })
        self.bam2fastq.on('end', self.set_output, 'bam2fastq')
        self.bam2fastq.on('end', self.make_config_run)
        self.bam2fastq.run()

    def make_config_run(self):
        self.make_config.set_options({
            "fastq_dir": self.bam2fastq.output_dir
        })
        self.make_config.on("end", self.set_output, 'config')
        self.make_config.on("end", self.ngsqc_stat_run)
        self.make_config.run()

    def ngsqc_stat_run(self):
        with open(os.path.join(self.make_config.output_dir, "fastq_list.txt"), 'r') as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split('\t')
                if os.path.getsize(temp[1]) < 3000:  # 提取的fastq大小小于3kb不计算
                    continue
                ngsqc_stat = self.add_tool("wgs.ngsqc_stat")
                ngsqc_stat.set_options({
                    "fastq_l": temp[1],
                    "fastq_r": temp[1],
                    "sample_name": temp[0]
                })
                self.ngsqc_stat_tools.append(ngsqc_stat)
            for j in range(len(self.ngsqc_stat_tools)):
                self.ngsqc_stat_tools[j].on("end", self.set_output, 'qc_stat')
            if self.ngsqc_stat_tools:
                if len(self.ngsqc_stat_tools) > 1:
                    self.on_rely(self.ngsqc_stat_tools, self.soap_denovo_run)
                elif len(self.ngsqc_stat_tools) == 1:
                    self.ngsqc_stat_tools[0].on('end', self.soap_denovo_run)
                for tool in self.ngsqc_stat_tools:
                    gevent.sleep(1)
                    tool.run()
            else:
                self.set_error("根据筛选条件从bam中没有提取出相应的fastq文件，修改参数重新运行！", code="24500201")

    def soap_denovo_run(self):
        for file_ in os.listdir(self.make_config.output_dir):
            if file_ == "fastq_list.txt":
                continue
            else:
                soap_denovo = self.add_tool("wgs.soap_denovo")
                soap_denovo.set_options({
                    "config_file": os.path.join(self.make_config.output_dir, file_)
                })
                self.soap_denovo_tools.append(soap_denovo)
        for j in range(len(self.soap_denovo_tools)):
            self.soap_denovo_tools[j].on("end", self.set_output, 'soap_denovo')
        if self.soap_denovo_tools:
            if len(self.soap_denovo_tools) > 1:
                self.on_rely(self.soap_denovo_tools, self.assembly_stat_run)
            elif len(self.soap_denovo_tools) == 1:
                self.soap_denovo_tools[0].on('end', self.assembly_stat_run)
            for tool in self.soap_denovo_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("soap_denovo_tools列表为空！", code="24500207")

    def assembly_stat_run(self):
        self.assembly_stat.set_options({
            "denovo_scafseq":
                self.soap_denovo_tools[0].output_dir + "/{}.denovo.scafSeq".format(self.option("sample_id")),
            "denovo_scafstatistics":
                self.soap_denovo_tools[0].output_dir + "/{}.denovo.scafStatistics".format(self.option("sample_id"))
        })
        self.assembly_stat.on("end", self.set_output, "assembly_stat")
        self.assembly_stat.on("end", self.end)
        self.assembly_stat.run()

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

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'bam2fastq':
            self.linkdir(obj.output_dir, 'bam2fastq')
        elif event['data'] == 'config':
            self.linkdir(obj.output_dir, 'config')
        elif event['data'] == 'assembly_stat':
            self.linkdir(obj.output_dir, 'assembly_stat')
        elif event['data'] == 'qc_stat':
            self.linkdir(obj.output_dir, 'qc_stat')
        elif event['data'] == 'soap_denovo':
            self.linkdir(obj.output_dir, 'soap_denovo')
        elif event['data'] == 'anno':
            self.linkdir(obj.output_dir, 'anno')
        else:
            pass

    def run(self):
        super(AssemblyModule, self).run()
        self.bam2fastq_run()

    def end(self):
        super(AssemblyModule, self).end()
