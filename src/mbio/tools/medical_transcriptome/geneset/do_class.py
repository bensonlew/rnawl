# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan, qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest
import pandas as pd

class DoClassAgent(Agent):
    '''
    last_modify: 2019.04.16
    '''
    def __init__(self, parent):
        super(DoClassAgent, self).__init__(parent)
        options = [
            {'name': 'do_ids', 'type': 'string', 'default': None},
            {'name': 'geneset_names', 'type': 'string', 'default': None},
            {'name': 'database', 'type': 'string', 'default': 'do'},
            {'name': 'source', 'type': 'string'},
            {"name": "do_version", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.step.add_steps('do_annot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.do_annot.start()
        self.step.update()

    def step_end(self):
        self.step.do_annot.finish()
        self.step.update()

    def check_options(self):
        if len(self.option("geneset_names").split(";")) != len(self.option("do_ids").split(";")):
            self.set_error("geneset_names 数量和 do_ids数量不一致")



    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(DoClassAgent, self).end()

class DoClassTool(Tool):
    def __init__(self, config):
        super(DoClassTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.do_annotation_py = os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/annotation/do_class.py')


    def run(self):
        super(DoClassTool, self).run()
        for name, go_ids in zip(self.option("geneset_names").split(","), self.option("do_ids").split(",")):
            self.run_do_annotation(go_ids, name)
        self.set_output()
        self.end()

    def run_do_annotation(self, do_file, name):
        cmd = '{} {} {} {}'.format(
            self.python,
            self.do_annotation_py,
            do_file,
            name + ".tsv"
        )
        cmd_name = 'run_do_class_{}'.format(name.lower())
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run %s, abord', variables=(n), code="33710402")

    def set_output(self):

        if len(self.option("geneset_names").split(",")) == 1:

            name1 = self.option("geneset_names")
            table1 = pd.read_table(name1 + ".tsv", header=0, sep="\t")
            table1.rename(columns={"DO1 (Lev1)": "Term Type",
                                   "DO Term (Lev2)": "DO Name",
                                   "DO ID (Lev2)": "DO ID",
                                   "Seq Number": name1 + " numbers",
                                   "Percent": name1 + " percent",
                                   "Seq List": name1 + " seqs"}, inplace=True)
            a = table1
            a[name1 + " numbers"]=a[name1 + " numbers"].fillna(0).map(int)
            a[name1 + " percent"] =  a[name1 + " percent"].fillna(0)
            a[name1 + " seqs"] =  a[name1 + " seqs"].fillna("")
            b = a.sort_values(["Term Type", "DO Name"])
            b.to_csv(self.output_dir + "/do_level2_class.xls", sep="\t", index=False)

        elif len(self.option("geneset_names").split(",")) >= 2:
            ## merge 2 class table
            name_list = self.option("geneset_names").split(",")

            name1 = name_list[0]
            table1 = pd.read_table(name1 + ".tsv", header=0, sep="\t")
            table1.rename(columns={"DO1 (Lev1)": "Term Type",
                                   "DO Term (Lev2)": "DO Name",
                                   "DO ID (Lev2)": "DO ID",
                                   "Seq Number": name1 + " numbers",
                                   "Percent": name1 + " percent",
                                   "Seq List": name1 + " seqs"}, inplace=True)

            a = table1
            for name2 in name_list[1:]:
                table2 = pd.read_table(name2 + ".tsv", header=0, sep="\t")
                table2.rename(columns={"DO1 (Lev1)": "Term Type",
                                       "DO Term (Lev2)": "DO Name",
                                       "DO ID (Lev2)": "DO ID",
                                       "Seq Number": name2 + " numbers",
                                       "Percent": name2 + " percent",
                                       "Seq List": name2 + " seqs"}, inplace=True)
                a = pd.merge(a, table2, how='outer')

            for name in name_list:
                a[name + " numbers"]=a[name + " numbers"].fillna(0).map(int)
                a[name + " percent"] =  a[name + " percent"].fillna(0)
                a[name + " seqs"] =  a[name + " seqs"].fillna("")
            b = a.sort_values(["Term Type", "DO Name"])
            b.to_csv(self.output_dir + "/do_level2_class.xls", sep="\t", index=False)
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        '''
        data = {
            'id': 'do_class_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'medical_transcriptome.geneset.do_class',
            'instant': False,
            'options': {
                'do_ids': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class/do/x;/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class/do/x2',
                'geneset_names': 'set1;set2',
            }
        }
        '''
        data = {
            'id': 'do_class_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'medical_transcriptome.geneset.do_class',
            'instant': False,
            'options': {
                'do_ids': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class/do/x',
                'geneset_names': 'set1',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
