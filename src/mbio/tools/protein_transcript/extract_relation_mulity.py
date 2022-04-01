# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import unittest


class ExtractRelationMulityAgent(Agent):
    def __init__(self, parent):
        super(ExtractRelationMulityAgent, self).__init__(parent)
        options = [
            {"name": "pep", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
            {"name": "relation_file", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
            {"name": "outlist", "type": "outfile", "format": "itraq_and_tmt.common"},  # 输出的有关联关系的list文件
            ]
        self.add_option(options)
        self.step.add_steps('relation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.relation.start()
        self.step.update()

    def step_end(self):
        self.step.relation.finish()
        self.step.update()

    def check_options(self):
        if not os.path.exists(self.option('relation_file').prop['path']):
            raise OptionError("没有正确传入biomart文件")
        if not os.path.exists(self.option('pep').prop['path']):
            raise OptionError("没有正确传入蛋白序列文件")
        return True

    def set_resource(self):
        self._cpu = 15
        self._memory = '40G'

    def end(self):
        super(ExtractRelationMulityAgent, self).end()


class ExtractRelationMulityTool(Tool):
    def __init__(self, config):
        super(ExtractRelationMulityTool, self).__init__(config)
        self.python_path = '/program/Python/bin/python'
        self.extract = self.config.PACKAGE_DIR + "/protein_transcript/extract_relation_mulity.py"

    def extract_relation(self):
        cmd = '{} {} '.format(self.python_path, self.extract)
        cmd += '{} {} '.format(self.option("pep").prop['path'], self.option("relation_file").prop['path'])
        cmd_name = 'extract_relation_mulity'
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
                self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004101")
        else:
            self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004102")

    def set_output(self):
        all_files = os.listdir(self.work_dir)
        all_files = [self.work_dir + '/' + each for each in all_files ]
        for each in all_files:
            if each.endswith('g2p.list'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)
                self.option('outlist').set_path(link)

    def run(self):
        super(ExtractRelationMulityTool, self).run()
        self.extract_relation()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87'
        data = {
            "id": "Extract_relation_" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "protein_transcript.extract_relation",
            "instant": False,
            "options": dict(
                pep = test_dir + "/" + "cds/Homo_sapiens.GRCh37.pep.fa",
                relation_file = test_dir + "/" + "biomart/Homo_sapiens.GRCh37.biomart",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
