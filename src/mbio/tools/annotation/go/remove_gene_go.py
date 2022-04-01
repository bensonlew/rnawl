# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess


class RemoveGeneGoAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: gaohao
    last_modified: 20181114
    """

    def __init__(self, parent):
        super(RemoveGeneGoAgent, self).__init__(parent)
        options = [
            {"name": "go1234level_out", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_anno", "type": "infile", "format": "sequence.profile_table"},  # 基因注释file
            {"name": "gogene_out", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.step.add_steps('go_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.go_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.go_annotation.finish()
        self.step.update()

    def check_options(self):
        if not self.option("go1234level_out").is_set:
            raise OptionError("Must provide go1234level_out result file !", code="31204801")
        if not self.option("gene_anno").is_set:
            raise OptionError("Must provide gene_anno result file !", code="31204802")

    def set_resource(self):
        self._cpu = 2
        self._memory = '30G'

    def end(self):
        super(RemoveGeneGoAgent, self).end()


class RemoveGeneGoTool(Tool):

    def __init__(self, config):
        super(RemoveGeneGoTool, self).__init__(config)


    def run(self):
        super(RemoveGeneGoTool, self).run()
        self.run_remove_gene()
        self.set_output()
        self.end()

    def run_remove_gene(self):
        cmd1 = '/program/Python/bin/python {}/annotation/go/go_remove_gene.py'.format(self.config.PACKAGE_DIR)
        cmd1 += ' -a %s -i %s -o %s' % (self.option("gene_anno").prop['path'],self.option("go1234level_out").prop['path'], self.work_dir + '/go1234level_statistics_new.xls')
        self.logger.info("运行run_remove_gene.py")
        command = self.add_command('run_go_anno_all', cmd1).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_remove_gene succeed")
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("run_remove_gene failed", code="31204801")

    def set_output(self):
        for i in ["go1234level_statistics_new.xls"]:
            if os.path.exists(self.output_dir + '/' + i):
                os.remove(self.output_dir + '/' + i)
            os.link(self.work_dir + '/' + i, self.output_dir + '/' + i)
        self.option('gogene_out',self.output_dir + '/go1234level_statistics_new.xls')