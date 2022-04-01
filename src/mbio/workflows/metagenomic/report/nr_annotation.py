# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

"""群落组成分析workflow"""
import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
import unittest,gevent
from comm_table import CommTableWorkflow


class NrAnnotationWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NrAnnotationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "query", "type": "infile", "format":"sequence.fasta"},
            {"name": "upload", "type": "string", "default": "best_hit"},

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.nr = self.add_module('align.meta_diamond')
        self.anno = self.add_module('annotation.annotation_stat')
        self.anno_nr_lca = self.add_module('annotation.annotation_stat')
        self.anno_nr_deunclass = self.add_module('annotation.annotation_stat')
        #self.ncbi_taxon = self.add_tool('annotation.mg_ncbi_taxon')
        self.anno_stat = []
        self.merge = []

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query")
        return True

    def run_nr(self):
        opts = {
            'query': self.option('query'),
            'query_type': "prot",
            'database': 'nr',
            'lines': '50000',
            "target_num": 5,
        }
        self.nr.set_options(opts)
        self.nr.on("end", self.run_anno_stat)
        self.nr.run()

    def run_anno(self):
        opts = {
            'nr_xml_dir': self.nr.option('outxml_dir'),
            "nr_method": "best_hit",
        }
        self.anno.set_options(opts)
        self.anno_stat.append(self.anno)

    def run_anno_nr_lca(self):
        opts = {
            'nr_xml_dir': self.nr.option('outxml_dir'),
            "nr_method": "lca",
        }
        self.anno_nr_lca.set_options(opts)
        self.anno_stat.append(self.anno_nr_lca)

    def run_anno_nr_deunclass(self):
        opts = {
            'nr_xml_dir': self.nr.option('outxml_dir'),
            "nr_method": "deunclassied",
        }
        self.anno_nr_deunclass.set_options(opts)
        self.anno_stat.append(self.anno_nr_deunclass)

    def run_anno_stat(self):
        self.run_anno()
        self.run_anno_nr_lca()
        self.run_anno_nr_deunclass()
        for module in self.anno_stat:
            module.run()
            gevent.sleep(0)
        self.on_rely(self.anno_stat, self.set_output)

    def set_output(self):
        for tool in self.anno_stat:
            for file in os.listdir(tool.output_dir):
                nr_align_path = tool.output_dir + "/" + file
                nr_align_basename = os.path.basename(nr_align_path)
                if nr_align_basename == 'best_hit.xls':
                    nr_dir = self.output_dir + "/nr"
                    if os.path.exists(nr_dir):
                        shutil.rmtree(nr_dir)
                    os.mkdir(nr_dir)
                    os.link(tool.output_dir+"/best_hit.xls", nr_dir + "/nr_align_table.xls")
                elif nr_align_basename == 'lca.xls':
                    nr_lca = self.output_dir + "/nr_lca"
                    if os.path.exists(nr_lca):
                        shutil.rmtree(nr_lca)
                    os.mkdir(nr_lca)
                    os.link(tool.output_dir+"/lca.xls", nr_lca + "/nr_align_table.xls")
                elif nr_align_basename == 'deunclassied.xls':
                    nr_deunclassied = self.output_dir + "/nr_deunclassied"
                    if os.path.exists(nr_deunclassied):
                        shutil.rmtree(nr_deunclassied)
                    os.mkdir(nr_deunclassied)
                    os.link(tool.output_dir+"/deunclassied.xls", nr_deunclassied + "/nr_align_table.xls")
                else:
                    pass
        self.end()

    def send_file(self):
        """
        上传结果文件到指定的位置
        """
        self.sheet.output = self.option('upload')
        dir_up = self.output_dir
        repaths = [
            ["nr/nr_align_table.xls", "", "物种序列比对结果", 0, "120037"],
            ["nr_lca/nr_align_table.xls", "", "物种序列比对结果", 0, "120037"],
            ["nr_deunclassied/nr_align_table.xls", "", "物种序列比对结果", 0, "120037"]
        ]
        sdir = self.add_upload_dir(dir_up)
        sdir.add_relpath_rules(repaths)

    def end(self):
        self.send_file()
        super(NrAnnotationWorkflow, self).end()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_nr()
        super(NrAnnotationWorkflow, self).run()
