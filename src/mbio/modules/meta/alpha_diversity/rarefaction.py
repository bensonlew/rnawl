# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'
# last_modify:2020.11.11

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.bacgenome.common import sum_stat
import unittest

class RarefactionModule(Module):
    def __init__(self, work_id):
        super(RarefactionModule, self).__init__(work_id)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table,meta.otu.tax_summary_dir"},  # 输入文件
            {"name": "indices", "type": "string", "default": "sobs,shannon"},  # 指数类型
            {"name": "freq", "type": "int", "default": 100},  # 取样频数
            {"name": "level", "type": "string", "default": "otu"}  # level水平
        ]
        self.add_option(options)
        self.rarefaction_list = []

    def check_options(self):
        if not self.option("otu_table").is_set:
            raise OptionError("参数otu_table不能为空", code="O32705201")
        return True

    def run_rarefaction(self):
        self.otu_path = self.option("otu_table").prop['path']
        options = {
            'otu_table': self.option('otu_table'),
            'indices': self.option('indices'),
            'freq': self.option('freq')
        }
        indices_list = self.option("indices").split(",")
        for f in indices_list:
            rarefaction = self.add_tool("meta.alpha_diversity.rarefaction")
            rarefaction.set_options({
                'otu_table': self.otu_path,
                'indices': f,
                'freq': self.option('freq'),
                'level': self.option('level')
            })
            self.rarefaction_list.append(rarefaction)
        if len(self.rarefaction_list) > 1:
            self.on_rely(self.rarefaction_list, self.set_output)
        else:
            self.rarefaction_list[0].on('end', self.set_output)
        for tool in self.rarefaction_list:
            tool.run()

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
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))

    def set_output(self):
        for tool in self.rarefaction_list:
            if os.path.exists(tool.output_dir):
                try:
                    file_name = os.listdir(tool.output_dir)[0]
                    if (os.path.exists(tool.output_dir + '/' + file_name)) and (len(os.listdir(tool.output_dir + '/' + file_name)) > 0):
                        self.linkdir(tool.output_dir + '/' + file_name, file_name)
                except:
                    self.logger.info("{} hava none result".format(tool.output_dir))
        self.end()

    def run(self):
        super(RarefactionModule, self).run()
        if self.option("indices") != "":
            self.run_rarefaction()
        else:
            self.end()

    def end(self):
        super(RarefactionModule, self).end()
