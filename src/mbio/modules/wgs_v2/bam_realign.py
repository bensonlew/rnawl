# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modify:20190220

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import json


class BamRealignModule(Module):
    """
     bam_realign module，用于将bam文件分别以tool的形式投出去。
    """
    def __init__(self, work_id):
        super(BamRealignModule, self).__init__(work_id)
        options = [
            {"name": "bam_list", "type": 'infile', "format": "wgs_v2.bam_list"},
            {"name": "fa_file", "type": "infile", "format": "sequence.fasta"},
            {"name": 'bam_realign_list', "type": "infile", 'format': "wgs_v2.bam_list"}
        ]
        self.add_option(options)
        self.bam_realign_tools = []

    def check_options(self):
        if not self.option("bam_list"):
            raise OptionError("请输入bam_list文件")
        return True

    def run_bam_realign(self):
        with open(self.option("bam_list").prop["path"]) as f:
            lines = f.readlines()
            for line in lines:
                sample_name = line.strip().split("\t")[0]
                sample_path = line.strip().split("\t")[1]
                bam_realign = self.add_tool("wgs_v2.bam_realign")
                options = ({
                    "bam_file": sample_path,
                    "fa_file": self.option("fa_file").prop["path"],
                    "name": sample_name + ".realign.bam"
                })
                bam_realign.set_options(options)
                self.bam_realign_tools.append(bam_realign)
            for j in range(len(self.bam_realign_tools)):
                self.bam_realign_tools[j].on("end", self.set_output, 'bam_realign_dir')
            if self.bam_realign_tools:
                if len(self.bam_realign_tools) > 1:
                    self.on_rely(self.bam_realign_tools, self.end)
                elif len(self.bam_realign_tools) == 1:
                    self.bam_realign_tools[0].on('end', self.end)
            else:
                raise Exception("bam_realign_tools列表为空！")
            for tool in self.bam_realign_tools:
                gevent.sleep(1)
                tool.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'bam_realign_dir':
            self.linkdir(obj.output_dir, 'bam_realign_dir')
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
        super(BamRealignModule, self).run()
        self.run_bam_realign()

    def end(self):
        self.set_realign_list()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(BamRealignModule, self).end()

    def set_realign_list(self):
        if len(self.bam_realign_tools) > 1:
            dir_path = self.output_dir + "/bam_realign_dir"
        else:
            # dir_path = self.bam_realign_tools[0] + "/bam_realign_dir"
            dir_path = self.bam_realign_tools[0].output_dir  # modified by hd 20190904
        with open(self.output_dir + "/bam.list", 'w') as w:
            for m in os.listdir(dir_path):
                n = re.match('(.*)\.realign\.bam$', m)
                if n:
                    w.write('{}\t{}\n'.format(n.group(1), dir_path + "/{}".format(m)))
        self.option("bam_realign_list").set_path(self.output_dir + "/bam.list")


