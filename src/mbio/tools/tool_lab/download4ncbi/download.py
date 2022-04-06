# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
__author__ = 'zjx'


class DownloadAgent(Agent):
    """
    Expression clustering analysis
    """
    def __init__(self, parent):
        super(DownloadAgent, self).__init__(parent)
        options = [
            {'name': 'id', 'type': 'string'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('10')

    def end(self):
        super(DownloadAgent, self).end()


class DownloadTool(Tool):
    """
    Expression clustering analysis
    """
    def __init__(self, config):
        super(DownloadTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.download_path = self.config.PACKAGE_DIR + '/tool_lab/download4ncbi.py'
        self.miniconda = software_dir + '/bioinfo/ref_rna_v3/gene_fusion/miniconda3/bin'
        self.set_environ(PATH=self.miniconda)

    def id_list(self):
        with open(os.path.join(self.work_dir, 'ID_list'), 'w') as i:
            i.write(self.option('id') + '\n')

    def download(self):
        cmd = '{} {}'.format(self.python_path, self.download_path)
        cmd += ' -ID_list {}'.format(os.path.join(self.work_dir, 'ID_list'))
        cmd += ' -output {}'.format(self.output_dir)
        cmd += ' -sorfware {}'.format(self.config.SOFTWARE_DIR)
        cmd_name = 'download4ncbi'
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
                self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35002902")
        else:
            self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35002903")

    def set_output(self):
        all_files = glob.glob(self.option("output") + '/PCA.xls')
        all_files += glob.glob(self.option("output") + '/Explained_variance_ratio.xls')
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(DownloadTool, self).run()
        self.id_list()
        self.download()
        # self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "download4ncbi" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "prok_rna.download4ncbi.download",
            "instant": False,
            "options": dict(
                id='GCF_004214875.1',
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
