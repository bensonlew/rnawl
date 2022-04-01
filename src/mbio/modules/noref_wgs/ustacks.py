# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last_modify:20181227

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import json


class UstacksModule(Module):
    """
    ustacks module,其中包含一个tool，即ustacks
    输入文件是一个fastq的列表
    样本名     样本fastq文件路径
    1080    /mnt/ilustre/users/tong.wang/Project/MJ20180611054-lihaiwen/var20180821/02.uniform/1080.fastq.gz
    """
    def __init__(self, work_id):
        super(UstacksModule, self).__init__(work_id)
        options = [
            {"name": "fastq_list", "type": 'infile', 'format': "noref_wgs.fastq_list"}
        ]
        self.add_option(options)
        self.ustacks_tools = []
        self.fastq_files = []

    def check_options(self):
        if not self.option("fastq_list").is_set:
            raise OptionError("please input sample's fastq list file!", code="25500303")

    def ustacks_run(self):
        i = 1
        for fastq in self.fastq_files:
            ustacks = self.add_tool("noref_wgs.ustacks")
            ustacks.set_options({
                "fastq": fastq,
                "id": i,
                "max_distance": 6,
                "min_depth": 2,
                "num_threads": 16
            })
            i += 1
            self.ustacks_tools.append(ustacks)
        for j in range(len(self.ustacks_tools)):
            self.ustacks_tools[j].on("end", self.set_output, "ustacks")
        if self.ustacks_tools:
            if len(self.ustacks_tools) > 1:
                self.on_rely(self.ustacks_tools, self.end)
            elif len(self.ustacks_tools) == 1:
                self.ustacks_tools[0].on('end', self.end)
        else:
            self.set_error("ustacks_tools列表为空！", code="25500303")
        for tool in self.ustacks_tools:
            gevent.sleep(1)
            tool.run()

    def set_sample_list(self):
        """
        将fastq list中按照样本进行排序，当然源头也要按照样本进行sort排序
        :return:
        """
        samples = []
        fastq_list = {}
        with open(self.option("fastq_list").prop['path'], "r") as r:
            data = r.readlines()
            for line in data:
                temp = line.split("\t")
                samples.append(temp[0])
                fastq_list[temp[0]] = temp[1]
        samples.sort()
        for m in samples:
            self.fastq_files.append(fastq_list[m])

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'ustacks':
            self.linkdir(obj.output_dir)
        else:
            pass

    def linkdir(self, dirpath, dirname=None):
        allfiles = os.listdir(dirpath)
        if dirname:
            newdir = os.path.join(self.output_dir, dirname)
        else:
            newdir = self.output_dir
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
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        super(UstacksModule, self).run()
        self.set_sample_list()
        self.ustacks_run()

    def end(self):
        self.make_ustack_list()
        super(UstacksModule, self).end()

    def make_ustack_list(self):
        with open('{}/ustacks.list'.format(self.output_dir), 'w') as w:
            for m in os.listdir(self.output_dir):
                result = re.match('(.*)\.check$', m)
                if result:
                    w.write('{}\t{}/{}\n'.format(result.group(1), self.output_dir, result.group(1)))

        allfiles = os.listdir(self.output_dir)
        total_data = os.path.join(self.output_dir, "tags.stat")
        with open(total_data, "w") as w:
            w.write("#sampleID\ttotal tags\ttotal depth\taverage depth\taverage length\tdep5\tdep10\n")
            for file in allfiles:
                file_path = os.path.join(self.output_dir, file)
                if file.endswith(".tags.stat"):
                    with open(file_path, 'r')as fr:
                        lines = fr.readlines()
                        for line in lines:
                            if not re.match('#.*', line):
                                w.write(line)

