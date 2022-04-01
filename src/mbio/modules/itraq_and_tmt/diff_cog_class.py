# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
# last_modify:2019.06.24

from biocluster.module import Module
import os
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import unittest


class DiffCogClassModule(Module):
    def __init__(self, work_id):
        super(DiffCogClassModule, self).__init__(work_id)
        options = [
            {'name': 'diff_path', 'type': 'infile', 'format': 'itraq_and_tmt.common_dir'},
            {'name': 'cog_stat', 'type': 'infile', 'format': 'itraq_and_tmt.common'},
        ]
        self.add_option(options)
        self.tools = list()

    def run_cog_class(self):
        diff_files = glob.glob(os.path.join(self.option('diff_path').prop['path'], '*_vs_*'))
        for diff in diff_files:
            options = dict(
                diff_file = diff,
                cog_stat = self.option('cog_stat')
            )
            class_tool = self.add_tool("itraq_and_tmt.export_cog_class")
            class_tool.set_options(options)
            self.tools.append(class_tool)

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        return True

    def set_output(self):
        shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        for cog_class in self.tools:
            out_dir = os.path.join(self.output_dir, os.path.basename(cog_class.option('diff_file').prop['path'].split('_diff.xls')[0]))
            os.mkdir(out_dir)
            for file in os.listdir(cog_class.output_dir):
                source = os.path.join(cog_class.output_dir,file)
                target = os.path.join(out_dir,file)
                os.link(source,target)
        self.end()

    def run(self):
        super(DiffCogClassModule, self).run()
        self.run_cog_class()
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def end(self):
        super(DiffCogClassModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "diff_cog_class_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "itraq_and_tmt.diff_cog_class",
            "instant": False,
            "options": dict(
                diff_path="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/Diff/output",
                cog_stat="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/ProteinAnnotation/output/cog/cog_summary.xls",
            )
        }
        # data['options']['method'] = 'rsem'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        # #
        # data['id'] += '1'
        # data['options']['method'] = 'salmon'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        #
        data['id'] += '_fyt'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()