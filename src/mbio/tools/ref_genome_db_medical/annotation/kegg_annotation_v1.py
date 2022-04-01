# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import xml.etree.ElementTree as ET
import subprocess


class KeggAnnotationAgent(Agent):
    """
    to perform KEGG annotation
    author:wangbixuan
    last_modified:20160729
    """

    def __init__(self, parent):
        super(KeggAnnotationAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "kegg_table", "type": "outfile", "format": "ref_rna_v2.kegg_table"},
        ]
        self.add_option(options)
        self.step.add_steps('kegg_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.kegg_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.kegg_annotation.finish()
        self.step.update()

    def check_options(self):
        if not self.option("blastout").is_set:
            raise OptionError("必须提供BLAST结果文件")
        else:
            pass

    def set_resource(self):
        self._cpu = 10
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./kegg_table.xls", "xls", "KEGG annotation table"],
            ["./pathway_table.xls", "xls", "Sorted pathway table"],
            ["./kegg_taxonomy.xls", "xls", "KEGG taxonomy summary"]
        ])
        result_dir.add_regexp_rules([
            [r"pathways/ko\d+", 'pdf', '标红pathway图']
        ])
        super(KeggAnnotationAgent, self).end()


class KeggAnnotationTool(Tool):

    def __init__(self, config):
        super(KeggAnnotationTool, self).__init__(config)
        self._version = "2.0"

    def run(self):
        super(KeggAnnotationTool, self).run()
        self.run_path_search()

    def run_path_search(self):
        cmd = '{}/program/Python/bin/python {}/bioinfo/annotation/scripts/pathSearch.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        cmd += ' %s %s %s %s' % (self.option('blastout').prop['path'], self.work_dir +
                                 '/kegg_table.xls', self.work_dir + '/pid.txt', self.work_dir + '/pathway_table.xls')
        self.logger.info("运行pathSearch.py")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            if os.path.exists(self.output_dir + '/kegg_table.xls'):
                os.remove(self.output_dir + '/kegg_table.xls')
            if os.path.exists(self.output_dir + '/pathway_table.xls'):
                os.remove(self.output_dir + '/pathway_table.xls')
            os.link(self.work_dir + '/kegg_table.xls',
                    self.output_dir + '/kegg_table.xls')
            self.option('kegg_table', self.output_dir + '/kegg_table.xls')
            os.link(self.work_dir + '/pathway_table.xls',
                    self.output_dir + '/pathway_table.xls')
        except subprocess.CalledProcessError:
            self.set_error("运行pathSearch.py出错")
        self.run_kegg_layer()

    def run_kegg_layer(self):
        cmd1 = '{}/program/Python/bin/python {}/bioinfo/annotation/scripts/keggLayer.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        cmd1 += ' %s %s %s' % (self.work_dir + '/pathway_table.xls', self.work_dir +
                               '/kegg_layer.xls', self.work_dir + '/kegg_taxonomy.xls')
        self.logger.info("运行keggLayer.py")
        self.logger.info(cmd1)
        try:
            subprocess.check_output(cmd1, shell=True)
            self.logger.info("运行keggLayer.py完成")
            if os.path.exists(self.output_dir + '/kegg_layer.xls'):
                os.remove(self.output_dir + '/kegg_layer.xls')
            if os.path.exists(self.output_dir + '/kegg_taxonomy.xls'):
                os.remove(self.output_dir + '/kegg_taxonomy.xls')
            os.link(self.work_dir + '/kegg_layer.xls',
                    self.output_dir + '/kegg_layer.xls')
            os.link(self.work_dir + '/kegg_taxonomy.xls',
                    self.output_dir + '/kegg_taxonomy.xls')
        except subprocess.CalledProcessError:
            self.set_error("运行keggLayer.py出错")
        self.run_get_pic()

    def run_get_pic(self):
        cmd2 = '{}/program/Python/bin/python {}/bioinfo/annotation/scripts/getPic.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        cmd2 += ' %s %s' % (self.work_dir + '/pid.txt',
                            self.work_dir + '/pathways')
        self.logger.info("运行getPic.py")
        self.logger.info(cmd2)
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info("运行getPic.py完成")
            if not os.path.exists(self.output_dir + '/pathways'):
                os.makedirs(self.output_dir + '/pathways')
            f = open(self.work_dir + '/pid.txt').read().split('\n')
            pids = []
            for record in f:
                if record != '':
                    pid = record.split('\t')[0]
                    pids.append(pid)
            for p in pids:
                if p != 'ko00312' and p != 'ko00351':
                    if os.path.exists(self.output_dir + '/pathways/' + p + '.pdf'):
                        os.remove(self.output_dir + '/pathways/' + p + '.pdf')
                    os.link(self.work_dir + '/pathways' + '/' + p + '.pdf',
                            self.output_dir + '/pathways/' + p + '.pdf')
        except subprocess.CalledProcessError:
            self.set_error("运行getPic.py出错")
        self.end()
