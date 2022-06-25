# coding=utf-8
# __author__ = 'linmeng.liu'
# last_modifiy = modified 2018.0816

import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'who'


class ProcrustesAgent(Agent):
    """
    procrustes description
    use qiime transform_coordinate_matrices.py
    """
    def __init__(self, parent):
        super(ProcrustesAgent, self).__init__(parent)
        options = [
            {'name': 'coord_ref', 'type': 'infile', 'default': 'None', 'format': 'metabolome.coord_matrices', 'required':True}, # 只能有一个
            {'name': 'coord_query', 'type': 'infile', 'default': 'None', 'format': 'metabolome.coord_matrices','required':True}, # 如有多个可逗号隔开，每个均会与ref进行计算
            {'name': 'random', 'type': 'int', 'default': '1000', 'format': 'None'},
            {'name': 'out', 'type': 'string', 'default': self.work_dir, 'format': 'None'}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('5')

    def end(self):
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Procrustes"]
            ])

        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
        ])
        """
        super(ProcrustesAgent, self).end()


class ProcrustesTool(Tool):
    """
    procrustes description
    """
    def __init__(self, config):
        super(ProcrustesTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = '/miniconda2/bin/python'
        self.procrustes = software_dir + '/miniconda2/bin/transform_coordinate_matrices.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run_procrustes(self):
        cmd = '{} {} '.format(self.python_path, self.procrustes)
        cmd += '-i {},{} '.format(self.option("coord_ref").prop["path"], self.option("coord_query").prop["path"])
        cmd += '-o {} '.format(self.option("out"))
        cmd += '-r {} '.format(self.option("random"))
        cmd_name = 'procrustes'
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
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32707301")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32707302")

    def set_output(self):
        target_files = glob.glob(self.option("out") + '/*transformed_reference.txt')
        target_files += glob.glob(self.option("out") + '/*transformed_q*.txt')
        target_files += glob.glob(self.option("out") + '/procrustes_results.txt')
        for each in target_files:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(ProcrustesTool, self).run()
        self.run_procrustes()
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
        data = {
            "id": "Procrustes" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "metabolome.metabset.procrustes",
            "instant": False,
            "options": dict(
                coord_ref="/mnt/ilustre/users/sanger-dev/sg-users/liulinmeng/metabolome/package/procrustes/procrustes_test/pcoa1.txt",
                coord_query="/mnt/ilustre/users/sanger-dev/sg-users/liulinmeng/metabolome/package/procrustes/procrustes_test/pcoa2.txt",
                random="1000",
                #out="/mnt/ilustre/users/sanger-dev/sg-users/liulinmeng/metabolome/package/procrustes/procrustes_test/ProcrustesToolTest",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
