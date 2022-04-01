# -*- coding: utf-8 -*-
# __author__ = 'wentian'
# last_modify:20180517

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import json


class TransgeneMappingModule(Module):
    """
    用于对转基因中mapping的地方
    """
    def __init__(self, work_id):
        super(TransgeneMappingModule, self).__init__(work_id)
        options = [
            {"name": "samples", "type": 'string'},
            {"name": "clean_path", "type": 'string'},
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)
        self.bwa_mem_tools = []
        self.samtools_view_tools = []
        self.samtools_merge_tools = []

    def check_options(self):
        if not self.option("samples"):
            raise OptionError("缺少samples参数", code="24501601")
        if not self.option("clean_path"):
            raise OptionError("缺少clean_path参数", code="24501602")
        if not self.option("ref_fa"):
            raise OptionError("缺少ref_fa参数", code="24501603")
        return True

    def bwa_mem_run(self):
        n = 1
        for i in json.loads(self.option("samples")):
            fastq_list_path = os.path.join(self.work_dir, "fastq.list")
            with open(fastq_list_path, "r") as fr:
                lines = fr.readlines()
                for line in lines[1:]:
                    data_split = line.strip().split("\t")
                    compare_samples = data_split[0]
                    compare_r1r2 = data_split[2]
                    if compare_samples == i:
                        fastq_l_path = compare_r1r2.strip().split(",")[0]
                        fastq_r_path = compare_r1r2.strip().split(",")[1]
                        bwa_mem = self.add_tool("wgs.bwa_mem")
                        bwa_mem.set_options({
                                 "fastq_l": fastq_l_path,
                                 "fastq_r": fastq_r_path,
                                 "sample_name": data_split[1],  # 带批次
                                 "ref_fa": self.option("ref_fa"),
                                 "num": str(n)
                             })
                        n += 1
                        self.bwa_mem_tools.append(bwa_mem)
        for j in range(len(self.bwa_mem_tools)):
            self.bwa_mem_tools[j].on("end", self.set_output, 'bwa_mem')
        if self.bwa_mem_tools:
            if len(self.bwa_mem_tools) > 1:
                self.on_rely(self.bwa_mem_tools, self.samtools_view_run)
            elif len(self.bwa_mem_tools) == 1:
                self.bwa_mem_tools[0].on('end', self.samtools_view_run)
        else:
            self.set_error("bwa_mem_tools列表为空！", code="24501609")
        for tool in self.bwa_mem_tools:
            gevent.sleep(1)
            tool.run()

    def samtools_view_run(self):
        fastq_list_path = os.path.join(self.work_dir, "fastq.list")
        with open(fastq_list_path, "r") as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                data_split = line.strip().split("\t")
                compare_samples = data_split[0]
                batch_name = data_split[1]
                for i in json.loads(self.option("samples")):
                    if i == compare_samples:
                        samtools_view = self.add_tool("wgs.samtools_view")
                        samtools_view.set_options({
                            "sam_file": self.output_dir + "/sam_dir/" + batch_name + ".sam"
                            if len(self.bwa_mem_tools) > 1 else self.bwa_mem_tools[0].output_dir + "/{}.sam".format(
                                batch_name)
                        })
                        self.samtools_view_tools.append(samtools_view)
        for j in range(len(self.samtools_view_tools)):
            self.samtools_view_tools[j].on("end", self.set_output, 'samtools_view')
        if self.samtools_view_tools:
            if len(self.samtools_view_tools) > 1:
                self.on_rely(self.samtools_view_tools, self.samtools_merge_run)
            elif len(self.samtools_view_tools) == 1:
                self.samtools_view_tools[0].on('end', self.samtools_merge_run)
        else:
            self.set_error("samtools_view_tools列表为空！", code="24501610")
        for tool in self.samtools_view_tools:
            gevent.sleep(1)
            tool.run()

    def samtools_merge_run(self):
        for i in json.loads(self.option("samples")):
            batch_list = []
            fastq_list_path = os.path.join(self.work_dir, "fastq.list")
            with open(fastq_list_path, "r") as fr:
                lines = fr.readlines()
                for line in lines[1:]:
                    data_split = line.strip().split("\t")
                    compare_samples = data_split[0]
                    batch_name = data_split[1] + ".bam"
                    batch_path = self.output_dir + "/bam_dir/" + batch_name
                    if compare_samples == i:
                        batch_list.append(batch_path)
            bam_list = ";".join(batch_list)
            samtools_merge = self.add_tool("wgs.samtools_merge")
            samtools_merge.set_options({
                "bam_list": bam_list,
                "specimen_id": i,
            })
            self.samtools_merge_tools.append(samtools_merge)
        for j in range(len(self.samtools_merge_tools)):
            self.samtools_merge_tools[j].on("end", self.set_output, 'samtools_merge')
        if self.samtools_merge_tools:
            if len(self.samtools_merge_tools) > 1:
                self.on_rely(self.samtools_merge_tools, self.end)
            elif len(self.samtools_merge_tools) == 1:
                self.samtools_merge_tools[0].on('end', self.end)
        else:
            self.set_error("samtools_merge_tools列表为空！", code="24501611")
        for tool in self.samtools_merge_tools:
            gevent.sleep(1)
            tool.run()

    def fastq_list(self):
        out = os.path.join(self.work_dir, "fastq.list")
        path_list = os.listdir(self.option("clean_path"))
        list_header = "sample_name" + "\t" + "batch" + "\t" + "path" + "\n"
        r1_list = []
        with open(out, "w") as f0:
            f0.write(list_header)
            for path in path_list:
                if re.match(".*\.1\.fastq\.gz$", path):
                    r1_list.append(path)
            for i in json.loads(self.option("samples")):
                for r1 in r1_list:
                    print r1_list
                    path_name = r1.strip().split(".")[0]
                    name_list = path_name.strip().split("-")
                    if name_list[0] == i:
                        r2_path = path_name + ".clean.2.fastq.gz"
                        print r2_path
                        if r2_path in path_list:
                            r1 = path_name + ".clean.1.fastq.gz"
                            r2 = path_name + ".clean.2.fastq.gz"
                            r1r2 = os.path.join(self.option("clean_path"), r1) + "," + os.path.join(self.option("clean_path"), r2)
                            write_content = name_list[0] + "\t" + path_name + "\t" + r1r2 + "\n"
                            f0.write(write_content)
                        else:
                            self.set_error("没有r2文件", code="24501612")

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
        if event['data'] == 'bwa_mem':
            self.linkdir(obj.output_dir, 'sam_dir')
        if event['data'] == 'samtools_view':
            self.linkdir(obj.output_dir, 'bam_dir')
        elif event['data'] == 'samtools_merge':
            self.linkdir(obj.output_dir, 'samtools_merge')
        else:
            pass

    def run(self):
        super(TransgeneMappingModule, self).run()
        self.fastq_list()
        # self.samtools_view_run()
        self.bwa_mem_run()

    def end(self):
        super(TransgeneMappingModule, self).end()
