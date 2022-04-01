# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:2017.9.29

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import unittest


class ProbioAnnoModule(Module):
    def __init__(self, work_id):
        super(ProbioAnnoModule, self).__init__(work_id)
        options = [
            {"name": "nr_gene_anno", "type": "infile", "format": "sequence.profile_table", "required": True},
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table", "required": True},
        ]
        self.add_option(options)
        self.anno_tool = self.add_tool("annotation.probio_anno")
        self.abun_tool = self.add_tool("annotation.probio_abun")

    def check_options(self):
        '''
        '''
        return True

    def probio_anno(self):
        self.anno_tool.set_options({
            "nr_gene_anno": self.option("nr_gene_anno")
        })
        self.anno_tool.on('end', self.probio_anno_stat)
        self.anno_tool.run()

    def probio_anno_stat(self):
        self.abun_tool.set_options({
            "probio_anno": self.anno_tool.option('probio_anno'),
            'reads_profile_table': self.option('reads_profile_table'),
        })
        self.abun_tool.on('end', self.set_output)
        self.abun_tool.run()

    def set_output(self):
        self.link_dir(self.anno_tool.output_dir)
        self.link_dir(self.abun_tool.output_dir)
        self.end()

    def run(self):
        super(ProbioAnnoModule, self).run()
        self.probio_anno()

    def end(self):
        super(ProbioAnnoModule, self).end()

    def link_dir(self, dir_path):
        files = os.listdir(dir_path)
        for each in files:
            oldfile = os.path.join(dir_path, each)
            newfile = os.path.join(self.output_dir, each)
            if os.path.exists(newfile):
                os.remove(newfile)
            os.link(oldfile, newfile)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "probio",
            "type": "module",
            "name": "annotation.probio_anno",
            "instant": True,
            "options": dict(
                nr_gene_anno="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/measurement/metag_v2/nr/gene_nr_anno.xls",
                reads_profile_table="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/annotation/gene_profile.reads_number.total.txt"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
