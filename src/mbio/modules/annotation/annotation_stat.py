# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'

from biocluster.module import Module
import os
import shutil,gevent
from biocluster.core.exceptions import OptionError
#import unittest
from mbio.packages.align.blast.xml2table import xml2table


class AnnotationStatModule(Module):
    def __init__(self, work_id):
        super(AnnotationStatModule, self).__init__(work_id)
        options = [
            {"name": "nr_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # nr比对结果文件夹
            {"name": "nr_method", "type": "string", "default": "best_hit"},  # 新增NR不同筛选结果方法,best_hit,lca,deunclassified
        ]
        self.add_option(options)
        self.cat_file = self.add_tool("sequence.cat_table")
        self.ncbi_tools = []
        self.number = 1

    def check_options(self):
        if not self.option("nr_xml_dir").is_set:
            raise OptionError("必须设置参数nr_xml_dir")
        return True

    def run_anno(self):
        self.anno_file = self.work_dir + "/anno_file"
        if os.path.exists(self.anno_file):
            pass
        else:
            os.mkdir(self.anno_file)
            xml_file = os.listdir(self.option('nr_xml_dir').prop['path'])
            for i in xml_file:
                file_path = os.path.join(self.option('nr_xml_dir').prop['path'], i)
                xml2table(file_path, self.anno_file + '/temp_blastable_{}.xls'.format(self.number))
                self.number +=1
        self.run_cat_file()

    def run_cat_file(self):
        """
        合并注释的结果文件
        """
        opts = {
            'fa_dir': self.anno_file,
            "prefix": self.option('nr_method'),
        }
        self.cat_file.set_options(opts)
        self.cat_file.on("end", self.set_output)
        self.cat_file.run()

    def set_output(self):
        """
        设置结果文件夹
        """
        if os.path.exists(self.output_dir + '/{}.xls'.format(self.option("nr_method"))):
            os.remove(self.output_dir + '/{}.xls'.format(self.option("nr_method")))
        os.link(self.cat_file.output_dir + '/{}.xls'.format(self.option("nr_method")), self.output_dir + '/{}.xls'.format(self.option("nr_method")))
        self.end()

    def run(self):
        super(AnnotationStatModule, self).run()
        self.run_anno()

    def end(self):
        super(AnnotationStatModule, self).end()