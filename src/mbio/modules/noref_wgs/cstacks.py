# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# modified 20200221

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
import os, re
import json
import time, datetime
import shutil
import gevent
from collections import defaultdict


class CstacksModule(Module):
    """
    无参变异检测Cstack部分流程，因为以前的版本由于样本比较多的时候，cstack的时候会超内存，所以现在增加新的流程.
    将所有的样本分成15分，分别进行cstack，然后将结果在合并成一个结果，用于后面的分析
    add by hd @20200220
    第一部分是将所有的样本分成20份
    第二部分是tool cstackv2,并发分别对20份结果进行cstack聚类
    第三部分是tool cstackmerge，将20份结果进行merge，将生成的20份cstacks结果进行合并成一个
    catalog.alleles.tsv.gz  catalog.tags.tsv.gz catalog.snps.tsv.gz 
    """
    def __init__(self, work_id):
        super(CstacksModule, self).__init__(work_id)
        options = [
            {"name": "group_list", "type": "string"},
            {"name": "ustacks_path", "type": "string"},  # 03.ustacks结果路径 ./Demo/01.stacks/03.ustacks/
            {"name": 'split_num', 'type': 'int', "default": 15},  # 将样本分成多少份进行cstack
        ]
        self.add_option(options)
        self.sample_run_list = defaultdict(list)
        self.ustacks_tools = []
        self.cstackmerge = self.add_tool("noref_wgs.cstacks_merge")

    def check_options(self):
        if not self.option("ustacks_path"):
            raise OptionError("请设置ustacks_path", code="200103")

    def set_sample_run_list(self):
        """
        进行样本的分组，如150个样本分成15份，每份是10个样本进行cstack聚类
        {'0':[]}
        add by hd @ 20200221
        """
        ustack_list = os.path.join(self.option("ustacks_path"), "ustacks.list")
        self.make_group(ustack_list)
        split_num = int(self.option('split_num'))
        with open(ustack_list, 'r') as r:
            n = 0
            for line in r:
                temp = line.strip().split('\t')
                tip_num = str(n % split_num)
                self.sample_run_list[tip_num].append(temp[1])
                n += 1
        self.logger.info(self.sample_run_list)
        # for keys in self.sample_run_list.keys():
        #     self.mkdirs(os.path.join(self.work_dir, keys))

    def make_group(self, ustack_list):
        """
        制造分组方案，这里默认只分一组数据
        """
        groups = {}
        if self.option("group_list"):
            with open(self.option("group_list"), "r")as fr:
                for line in fr:
                    if re.match('#.*', line):
                        continue
                    tmp = line.strip().split("\t")
                    if len(tmp) < 2:
                        self.set_error("分组文件第一列应为样本名，第二列应为分组名！", code="35500208")
                    groups[tmp[0]] = tmp[1]

        with open(ustack_list, 'r') as r, open(self.work_dir + "/group.list", 'w') as w, \
                open(self.work_dir + "/sample.list", 'w') as w1:
            for line in r:
                temp = line.strip().split('\t')
                if temp[0] in groups.keys():
                    group_name = groups[temp[0]]
                else:
                    group_name = str(1)
                w.write(temp[0] + '\t' + group_name + "\n")
                w1.write(temp[0] + '\n')

    def mkdirs(self, filepath):
        if os.path.exists(filepath):
            os.system('rm -rf {}'.format(filepath))
            os.mkdir(filepath)
        else:
            os.mkdir(filepath)

    def run_cstacks_1(self):
        """
        执行cstack聚类，分别生成20个文件夹，每个文件夹中10个样本进行聚类
        add by hd @ 20200221
        """
        for keys in self.sample_run_list.keys():
            ustacks = self.add_tool("noref_wgs.cstacksv2")
            ustacks.set_options({
                "ustacksfiles": json.dumps(self.sample_run_list[keys]),
                "outnum": str(keys),
            })
            self.ustacks_tools.append(ustacks)
        for j in range(len(self.ustacks_tools)):
            self.ustacks_tools[j].on("end", self.set_output, "ustacks")
        if self.ustacks_tools:
            if len(self.ustacks_tools) > 1:
                self.on_rely(self.ustacks_tools, self.run_cstacks_2)
            elif len(self.ustacks_tools) == 1:
                self.ustacks_tools[0].on('end', self.run_cstacks_2)
        else:
            self.set_error("ustacks_tools列表为空！", code="25500303")
        for tool in self.ustacks_tools:
            gevent.sleep(1)
            tool.run()
        pass

    def run_cstacks_2(self):
        """
        将run_cstacks_1生成的20份cstacks结果进行合并成一个catalog.alleles.tsv.gz  catalog.tags.tsv.gz catalog.snps.tsv.gz
        add by hd @ 20200221
        """
        if len(self.ustacks_tools) > 1:
            catalog_path = self.output_dir + "/cstackssplit/"
        else:
            catalog_path = self.ustacks_tools[0].output_dir
        self.cstackmerge.set_options({
            "catalog_path": catalog_path
        })
        self.cstackmerge.on('end', self.set_output, "ustacksmerge")
        self.cstackmerge.on('end', self.end)
        self.cstackmerge.run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
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
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'ustacks':
            self.linkdir(obj.output_dir, 'cstackssplit')
        elif event['data'] == 'ustacksmerge':
            self.linkdir(obj.output_dir, 'cstacksmerge')
        else:
            pass

    def run(self):
        super(CstacksModule, self).run()
        self.set_sample_run_list()
        self.run_cstacks_1()

    def end(self):
        super(CstacksModule, self).end()
