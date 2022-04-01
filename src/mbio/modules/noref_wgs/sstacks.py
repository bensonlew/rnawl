# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modify:201812

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import json


class SstacksModule(Module):
    """
    sstacks module,其中包含一个tool，即sstacks
    传入文件分别为第三步03.ustacks和第四步04.cstacks结果文件夹。
    """
    def __init__(self, work_id):
        super(SstacksModule, self).__init__(work_id)
        options = [
            {"name": "sample_loci_dir", "type": 'string'},  # 格式需要修改
            {"name": "catalog_dir", "type": 'string'},  # 格式需要修改
        ]
        self.add_option(options)
        self.sstacks_tools = []


    def sstacks_run(self):
        """
        这里获取样本有两种方法，一种是通过遍历获取该文件夹下所有样本的名字，另一种则是直接读取ustacks.list文件，这里采取这种方法
        20190102 第二种方法取决于是否有list的文件，在某些条件下没有，所以还是自己来生成
                """
        sample_list = []
        # file_path = self.option('sample_loci_dir') + "/" + "ustacks.list"
        # self.logger.info(file_path)
        # self.logger.info("----------------------------------")
        # print
        # with open(file_path) as f:
        #     lines = f.readlines()
        #     for line in lines:
        #         sample = line.strip().split("\t")[0]
        #         sample_list.append(sample)
        temp = os.listdir(self.option("sample_loci_dir"))
        for i in temp:
            if i.strip().split(".")[0] not in sample_list:
                if not i.startswith('ustacks.list') and not i.startswith('tags.stat'):
                    sample_list.append(i.strip().split(".")[0])
        for j in sample_list:
            self.logger.info("*************************")
            self.logger.info(j)
            sstacks = self.add_tool("noref_wgs.sstacks")
            sstacks.set_options(
                {
                    "catlog_dir": self.option("catalog_dir"),
                    "sample_loci_dir": self.option("sample_loci_dir"),
                    "sample_name": j
                }
            )
            self.sstacks_tools.append(sstacks)
        for j in range(len(self.sstacks_tools)):
            self.sstacks_tools[j].on("end", self.set_output, "sstacks_dir")
        if self.sstacks_tools:
            if len(self.sstacks_tools) > 1:
                self.on_rely(self.sstacks_tools, self.end)
            elif len(self.sstacks_tools) == 1:
                self.sstacks_tools[0].on('end', self.end)
        else:
            self.set_error("sstacks_tools列表为空！", code="25500503")
        for tool in self.sstacks_tools:
            gevent.sleep(1)
            tool.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'sstacks_dir':
            self.linkdir(obj.output_dir, 'sstacks_dir')
        else:
            pass

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

        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        super(SstacksModule, self).run()
        self.sstacks_run()

    def end(self):
        super(SstacksModule, self).end()
