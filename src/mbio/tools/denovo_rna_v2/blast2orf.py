# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'liubinxu'


class Blast2orfAgent(Agent):
    """
    predict orf from blast result
    """
    def __init__(self, parent):
        super(Blast2orfAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'fasta', 'format': 'denovo_rna_v2.fasta'},
            {'type': 'infile', 'name': 'blast_nr_xml', 'format': 'ref_rna_v2.blast_xml'},
            {'type': 'infile', 'name': 'blast_swissprot_xml', 'format': 'ref_rna_v2.blast_xml'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('10')

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
        super(Blast2orfAgent, self).end()


class Blast2orfTool(Tool):
    """
    predict orf from blast result
    """
    def __init__(self, config):
        super(Blast2orfTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        package_dir = self.config.PACKAGE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.blast2orf = package_dir + '/denovo_rna_v2/blast2orf.py'


    def run_blast2orf(self):
        cmd = '{} {} '.format(self.python_path, self.blast2orf)
        cmd += '{} '.format(self.option("blast_nr_xml").prop['path'] + "," + self.option("blast_swissprot_xml").prop['path'])
        cmd += '{} '.format(self.option("fasta").prop['path'])

        cmd_name = 'blast2orf'
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
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32007701")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32007702")

    def set_output(self):
        pass
        all_files = glob.glob(self.work_dir + '/blast_hit*')
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(Blast2orfTool, self).run()
        self.run_blast2orf()
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
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo_rna_v2'
        data = {
            "id": "Blast2orf" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.blast2orf",
            "instant": False,
            "options": dict(
                fasta=test_dir + "/" + "assemble_2/assemble_raw.fasta",
                blast_nr_xml=test_dir + "/" + "map_db/nr/blast.xml",
                blast_swissprot_xml=test_dir + "/" + "map_db/swissprot/blast.xml",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
