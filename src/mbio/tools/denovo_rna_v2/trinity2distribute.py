# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
__author__ = 'liubinxu'


class Trinity2DistributeAgent(Agent):
    """
    Trinity
    """
    def __init__(self, parent):
        super(Trinity2DistributeAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'distribute_cmd', 'format': 'denovo_rna_v2.common'},
            {'default': 10, 'type': 'int', 'name': 'cpu'},
            {'default': 0, 'type': 'float', 'name': 'mis_rate'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('20')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        super(Trinity2DistributeAgent, self).end()


class Trinity2DistributeTool(Tool):
    """
    Trinity
    """
    def __init__(self, config):
        super(Trinity2DistributeTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.parafly = software_dir + '/bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/trinity-plugins/ParaFly-0.1.0/bin/ParaFly'

    def run_parafly(self):
        cmd = '{} '.format(self.parafly)
        cmd += '-{} {} '.format("C", self.option("distribute_cmd").prop['path'])
        cmd += '-{} {} '.format("cpu", self.option("cpu"))
        cmd += '-{} {} '.format("failed_cmds", self.option("distribute_cmd").prop['path'] )
        cmd_name = 'trinity2distribute'
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
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32006401")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32006402")

    def set_output(self):
        pass

    def run(self):
        super(Trinity2DistributeTool, self).run()
        self.run_parafly()
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
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20180102/Single_de_tr_data2.2'
        data = {
            "id": "Trinity2Distribute" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.trinity2distribute",
            "instant": False,
            "options": dict(
                distribute_cmd=test_dir + "/" + "trinity_cmd_distribute10",
                cpu="20",
                mis_rate="0",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
