# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os,re
from biocluster.core.exceptions import OptionError
import subprocess
import pandas as pd
import math


class MgGoAnnotationAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: gaohao
    last_modified: 20181114
    """

    def __init__(self, parent):
        super(MgGoAnnotationAgent, self).__init__(parent)
        options = [
            {"name": "blast2go_annot", "type": "infile", "format": "annotation.go.blast2go_annot"},
            {"name": "go2level_out", "type": "outfile", "format": "annotation.go.level2"},
        ]
        self.add_option(options)
        self.step.add_steps('go_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 40  # 每次重运行增加内存40G by qingchen.zhang @ 20190527

    def step_start(self):
        self.step.go_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.go_annotation.finish()
        self.step.update()

    def check_options(self):
        if self.option("blast2go_annot").is_set:
            pass
        else:
            raise OptionError("Must provide blast2go_annot result file !", code="31205501")

    def set_resource(self):
        self._cpu = 5
        count = 0
        for index, line in enumerate(open(self.option('blast2go_annot').prop['path'],'r')):
            count += 1
        print count
        total = int(count) / (50000)
        total = math.ceil(total)
        self._memory = '{}G'.format(int(total) + 10)

    def end(self):
        super(MgGoAnnotationAgent, self).end()


class MgGoAnnotationTool(Tool):

    def __init__(self, config):
        super(MgGoAnnotationTool, self).__init__(config)
        self.annot_go = self.config.PACKAGE_DIR + "/prok_rna/goAnnot.py"

    def run(self):
        super(MgGoAnnotationTool, self).run()
        self.run_gomerge()
        self.run_annotation()
        self.run_go_gene_anno()
        self.set_output()
        self.end()

    def run_gomerge(self):
        cmd1 = 'miniconda2/bin/python {}/bioinfo/annotation/scripts/goMerge.py'.format(
            self.config.SOFTWARE_DIR)
        cmd1 += ' %s %s' % (self.option("blast2go_annot").prop['path'], 'GO.list')
        self.logger.info("运行mergeGO.py")
        self.logger.info(cmd1)
        #try:
            #subprocess.check_output(cmd1, shell=True)
        #except subprocess.CalledProcessError:
            #self.set_error('running mergeGO.py error', code="31205501")
        command = self.add_command("run_gomerge", cmd1)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行mergeGO.py")
        else:
            self.set_error('running mergeGO.py error', code="31205501")

    def run_annotation(self):
        self.run_filter(self.work_dir + '/GO.list',self.work_dir + '/GO_2.list')
        cmd2 = 'miniconda2/bin/python {} {}'.format(self.annot_go, self.work_dir + '/GO_2.list')
        self.logger.info("运行goAnnot.py")
        self.logger.info(cmd2)
        command = self.add_command("run_annotation", cmd2, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行run_annotation结束")
        elif command.return_code in [-9, 1]:  # change return_code by qingchen.zhang @ 20190610
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by qingchen.zhang @ 20190610
        else:
            self.set_error("运行run_annotation出错")

    def run_go_gene_anno(self):
        cmd3 = 'miniconda2/bin/python {}/annotation/go/go_anno.py'.format(
            self.config.PACKAGE_DIR)
        cmd3 += ' -i %s -o %s' % (self.work_dir + '/go1234level_statistics.xls',self.work_dir + '/all.go.annotation.xls')
        self.logger.info("运行go_anno.py")
        command = self.add_command("run_go_gene_anno", cmd3, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行run_go_gene_anno结束")
        elif command.return_code in [-9, 1]:  # change return_code by qingchen.zhang @ 20190610
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by qingchen.zhang @ 20190610
        else:
            self.set_error("运行run_go_gene_anno出错")

    def set_output(self):
        for i in ["go1234level_statistics.xls","all.go.annotation.xls"]:
            if os.path.exists(self.output_dir + '/' + i):
                os.remove(self.output_dir + '/' + i)
            os.link(self.work_dir + '/' + i, self.output_dir + '/' + i)

    def run_filter(self,file,output):
        table =pd.read_table(file,sep="\t",header=None)
        ll = re.compile('_1$')
        table[0] = table[0].replace(ll, '')
        table.to_csv(output, sep="\t", header=None, index=0)

