# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class PlantNameTranslationWorkflow(Workflow):
    """
    Used for convert RNA-seq expression units convert, count/CPM/RPM/TPM/FPKM/RPKM are implemented.

    RPM/CPM: Reads/Counts of exon model per million mapped reads
    RPM/CPM=Total exon reads/ Mapped reads(Millions)

    RPKM/FPKM: Reads/Fragments Per Kilobase of exon model per Million mapped reads
    RPKM/FPKM=Total exon reads/[Mapped reads(Millions)*Exon length(Kb)]

    TPM is like RPKM and FPKM, except the order of operation is switched.
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PlantNameTranslationWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'plant_list', 'type': 'string', 'default': None},
            {'name': 'name_type', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.plant_name_translation")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(PlantNameTranslationWorkflow, self).run()

    def run_tool(self):
        opts = {
            'name_type': self.option('name_type'),
            'plant_list': self.option('plant_list')
        }
        #
        #
        # if self.option('chinese_name'):
        #     opts = {
        #         'name_type': 'chinese',
        #         'plant_list': self.option('chinese_name')
        #     }
        # elif self.option('latin_name'):
        #     opts = {
        #         'name_type': 'latin',
        #         'plant_list': self.option('latin_name')
        #     }
        # elif self.option('english_name'):
        #     opts = {
        #         'name_type': 'english',
        #         'plant_list': self.option('english_name')
        #     }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):

        """
         保存结果标准化数据到mongo数据库中
        """
        name_trans = self.api.api('tool_lab.translation')
        name_trans.add_translation(self.option('main_id'), self.tool.option('translation_file').path)

        self.set_output()


    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "植物名字互译目录", 0],
        ])
        super(PlantNameTranslationWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.exp_norm import ExpNormWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "exp_norm" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.exp_norm",
            "options": dict(
                exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/known_seqs_count.matrix",
                convert_type="DESeq2",
                # float_num=4,
            )
        }
        wsheet = Sheet(data=data)
        wf =ExpNormWorkflow(wsheet)
        wf.sheet.id = 'exp_norm'
        wf.sheet.project_sn = 'exp_norm'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
