
# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last_modify:20180825

import os
import re
import gevent
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class GwasMvpModule(Module):
    """
    GWAS关联分析的module:gwas_analysis.py接口调用这个module
    包MVP的R的tool（Tool:trit_mvp.py）,
    可以让多性状多台机器并行
    """
    def __init__(self, work_id):
        super(GwasMvpModule, self).__init__(work_id)
        options = [
            {"name": "trit_list_path", "type": "string"},    # Tool: vcf2trit产生的list
            {"name": "pop_hapmap_path", "type": "string"},    # Tool: vcf2hapmap产生
        ]
        self.add_option(options)
        self.trit_mvp_tools = []
        self.end_times = 0

    def check_options(self):
        if not self.option("trit_list_path"):
            raise OptionError("必须有trit_list_path")
        if not self.option("pop_hapmap_path"):
            raise OptionError("必须有pop_hapmap_path")

    def run_trit_mvp(self):
        """
        并行投递trit_mvp.py的tool
        # 出很多结果的csv文件，需要链接到一个路径下，然后存这个路径到主表里
        """
        if not os.path.exists(self.option('trit_list_path')):
            self.set_error("trit_list文件：{}不存在".format(self.option('trit_list')))
        self.trit_l = []
        with open(self.option('trit_list_path'), 'r') as f:
            trit_list = f.readlines()
            self.trit_l = trit_list
            for i in trit_list:
                list_ = i.strip().split('\t')
                trit_name = list_[0]
                trt_path = list_[1]
                self.trit_mvp = self.add_tool("dna_evolution.trit_mvp")
                self.trit_mvp.set_options({
                    "trait_name": trit_name,
                    "trait_trt": trt_path,
                    "pop_hapmap": self.option('pop_hapmap_path')
                })
                self.logger.info(trit_name)
                self.trit_mvp_tools.append(self.trit_mvp)
        # self.logger.info(self.trit_mvp_tools)
        for j in range(len(self.trit_mvp_tools)):
            self.trit_mvp_tools[j].on("end", self.set_output, 'trit_mvp')
        if self.trit_mvp_tools:
            if len(self.trit_mvp_tools) > 1:
                self.on_rely(self.trit_mvp_tools, self.set_output)
            elif len(self.trit_mvp_tools) == 1:
                self.trit_mvp_tools[0].on('end', self.set_output, 'trit_mvp')
        else:
            raise Exception("trit_mvp_tools列表为空！")
        for tool in self.trit_mvp_tools:
            gevent.sleep(1)
            tool.run()

    def set_output(self, event):
        self.logger.info("---------")
        obj = event['bind_object']
        if event['data'] == 'trit_mvp':
            self.linkdir(obj.output_dir, self.output_dir)
        self.end_times += 1
        if len(self.trit_mvp_tools) == self.end_times:
            self.end()

    def end(self):
        super(GwasMvpModule, self).end()

    def run(self):
        super(GwasMvpModule, self).run()
        self.run_trit_mvp()

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
