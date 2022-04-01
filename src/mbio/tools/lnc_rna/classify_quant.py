# -*- coding: utf-8 -*-
# import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

import os
import unittest
#__author__ = 'fuwenyao'


class ClassifyQuantAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(ClassifyQuantAgent, self).__init__(parent)
        options = [
            {'name': 'new_lnc_rna', 'type': 'string'},#rna鉴定结果，新lnc_rna list
            {'name': 'known_lnc_rna', 'type': 'string'},#rna鉴定结果，已知lnc_rna list
            {'name': 'new_mrna', 'type': 'string'},#rna鉴定结果，新mrna list
            {'name': 'known_mrna', 'type': 'string'},#rna鉴定结果，已知mrna list
            {'name': 'quant_results', 'type': 'string'},#表达量结果目录
            {'name': 't2g', 'type': 'string'},#trans_id to gene_id，from 组装结果目录
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('15')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(ClassifyQuantAgent,self).end()

class ClassifyQuantTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(ClassifyQuantTool, self).__init__(config)
        self.python_path = 'program/Python/bin/python'
        self.QcStatIntergration = self.config.PACKAGE_DIR + "/lnc_rna/get_classify_quant.py"

    def run_ClassifyQuant(self):
        quant_files=glob.glob(r"{}/*".format(self.option("quant_results")))
        for n, quant_file in enumerate(quant_files):
           if quant_file.endswith("matrix"):
             if os.path.basename(quant_file).startswith("transcript"):
                cmd = '{} {} '.format(self.python_path, self.QcStatIntergration)
                cmd += '-{} {} '.format("new", self.option("new_lnc_rna"))
                cmd += '-{} {} '.format("known", self.option("known_lnc_rna"))
                cmd += '-{} {} '.format("newm", self.option("new_mrna"))
                cmd += '-{} {} '.format("knownm", self.option("known_mrna"))
                cmd += '-{} {} '.format("quant",quant_file)
                cmd += '-{} {} '.format("out",os.path.basename(quant_file))
                cmd_name = 'classify_quant_{}'.format(n)
                command = self.add_command(cmd_name, cmd)
                command.run()
                self.wait()
                if command.return_code == 0:
                  self.logger.info("{} Finished successfully".format(cmd_name))
                elif command.return_code is None:
                  self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
                  command.rerun()
                  self.wait()
                  if command.return_code is 0:
                      self.logger.info("{} Finished successfully".format(cmd_name))
                  else:
                      self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
                else:
                  self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
             else:
                cmd = '{} {} '.format(self.python_path, self.QcStatIntergration)
                cmd += '-{} {} '.format("new", self.option("new_lnc_rna"))
                cmd += '-{} {} '.format("known", self.option("known_lnc_rna"))
                cmd += '-{} {} '.format("quant", quant_file)
                cmd += '-{} {} '.format("newm", self.option("new_mrna"))
                cmd += '-{} {} '.format("knownm", self.option("known_mrna"))
                cmd += '-{} {} '.format("out", os.path.basename(quant_file))
                cmd += '-{} {} '.format("t2g",self.option("t2g"))
                cmd_name = 'classify_quant_{}'.format(n)
                command = self.add_command(cmd_name, cmd)
                command.run()
                self.wait()
                if command.return_code == 0:
                    self.logger.info("{} Finished successfully".format(cmd_name))
                elif command.return_code is None:
                    self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
                    command.rerun()
                    self.wait()
                    if command.return_code is 0:
                        self.logger.info("{} Finished successfully".format(cmd_name))
                    else:
                        self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
                else:
                    self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

    def set_output(self):
        all_files = os.listdir(self.work_dir)
        for each in all_files:
            if each.endswith('matrix'):
                fname = os.path.basename(each)
                each = os.path.join(self.work_dir, fname)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

    def run(self):
        super(ClassifyQuantTool, self).run()
        self.run_ClassifyQuant()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        lncrna_dir='/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/others/test_for_quant'
        quant_dir='/mnt/ilustre/users/sanger-dev/workspace/20190325/Single_Quant+3302/Quant/output'
        t2g='/mnt/ilustre/users/sanger-dev/workspace/20190322/Single_assemble_3643_9041/Assemble/output/NewTranscripts/trans2gene'
        data = {
            "id": "classify_quant" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "lnc_rna.classify_quant",
            "instant": False,
            "options": dict(
                known_lnc_rna=lncrna_dir + "/" + "known_lncrna_ids.list",
                new_lnc_rna=lncrna_dir + "/" + "new_lncrna_ids.list",
                known_mrna=lncrna_dir + "/" + "known_mrna_ids.list",
                new_mrna=lncrna_dir + "/" + "new_mrna_ids.list",
                quant_results=quant_dir,
                t2g=t2g
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()