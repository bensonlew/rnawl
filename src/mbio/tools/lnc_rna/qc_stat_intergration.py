import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
#__author__ = 'fuwenyao'


class QcStatIntergrationAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(QcStatIntergrationAgent, self).__init__(parent)
        options = [
            {'name': 'stat_before', 'type': 'string'},
            {'name': 'stat_after', 'type': 'string',},
            {'name': 'out', 'type': 'string'},
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
        super(QcStatIntergrationAgent,self).end()

class QcStatIntergrationTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(QcStatIntergrationTool, self).__init__(config)
        self.python_path = 'miniconda2/bin/python'
        self.QcStatIntergration = self.config.PACKAGE_DIR + "/lnc_rna/get_qc_stat_intergration.py"

    def run_QcStatIntergration(self):
        cmd = '{} {} '.format(self.python_path, self.QcStatIntergration)
        cmd += '-{} {} '.format("before", self.option("stat_before"))
        cmd += '-{} {} '.format("after", self.option("stat_after"))
        cmd += '-{} {} '.format("out", self.option("out"))
        cmd_name = 'qc_stat_intergration'
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
            if each.endswith('results'):
                fname = os.path.basename(each)
                each = os.path.join(self.work_dir, fname)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

    def run(self):
        super(QcStatIntergrationTool, self).run()
        self.run_QcStatIntergration()
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
        test_dir_before='/mnt/ilustre/users/sanger-dev/workspace/20190226/Single_HiseqReadsStat_5971/HiseqReadsStat/output'
        test_dir_after='/mnt/ilustre/users/sanger-dev/workspace/20190227/Single_HiseqReadsStat_7969/HiseqReadsStat/output'
        data = {
            "id": "qc_stat_intergration" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "lnc_rna.qc_stat_intergration",
            "instant": False,
            "options": dict(
                stat_before=test_dir_before + "/" + "stat_results",
                stat_after=test_dir_after + "/" + "stat_results",
                out="qc_stat_results"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()