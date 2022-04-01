# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'gdq'


class FilterAnnotAgent(Agent):
    """
    filter xml result
    """
    def __init__(self, parent):
        super(FilterAnnotAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'blast', 'format': 'denovo_rna_v2.blast_xml'},
            {'type': 'infile', 'name': 'hmm', 'format': 'denovo_rna_v2.common'},
            {'default': 'blast', 'type': 'string', 'name': 'method'},
            {'default': 0.001, 'type': 'float', 'name': 'evalue'},
            {'default': 0, 'type': 'float', 'name': 'identity'},
            {'default': 0, 'type': 'float', 'name': 'similarity'},
            {'default': "", 'type': 'string', 'name': 'exclude_taxon'},
            {'type': 'outfile', 'name': 'outxml', 'format': 'denovo_rna_v2.blast_xml'},
            {'type': 'outfile', 'name': 'outtable', 'format': 'denovo_rna_v2.common'},
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
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(FilterAnnotAgent, self).end()


class FilterAnnotTool(Tool):
    """
    filter xml result
    """
    def __init__(self, config):
        super(FilterAnnotTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
 
    def run_filterxml(self):
        basename = os.path.basename(self.option('blast').prop['path']) + ".filter.xml"
        self.option('blast').filter_blast_xml(basename, self.option('evalue'), self.option('identity'), self.option('similarity'))
 
    def set_output(self):
        pass
        '''Example:
        diff_files = glob.glob(self.option("output") + '/*_vs_*.xls')
        diff_list = glob.glob(self.option("output") + '/*.DE.list')
        diff_summary = glob.glob(self.option("output") + '/*summary.xls')
        all_files = diff_files + diff_list + diff_summary
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
        '''

    def run(self):
        super(FilterAnnotTool, self).run()
        self.run_filterxml()
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
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/test_data2/annotation/output/blast_xml/'
        data = {
            "id": "Blastfilter" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.filterannot",
            "instant": False,
            "options": dict(
                blast=test_dir + "Trinity_vs_swissprot.xml",
                method="blast",
                evalue=0.000001,
                identity=90,
                similarity=90
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
