# -*- coding: utf-8 -*-
# __author__ = 'ysh'

import os
import re
import types
from biocluster.workflow import Workflow
from mainapp.models.mongo.bacgenome import Bacgenome
from bson import ObjectId
import datetime
from biocluster.file import download, exists


class TestApiWorkflow(Workflow):
    """

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TestApiWorkflow, self).__init__(wsheet_object)
        options = [
            # {"name": "specimens", "type": "string"},  #
            # {"name": "main_id", "type": "string"},
            #{"name": "task_id", "type": "string"},
            {"name":"test","type":"string"}

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.set_db()
        super(TestApiWorkflow, self).run()



    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.prephage = self.api.api('bacgenome.prephage')
        prephage_id = self.prephage.add_prephage("prephage")

        self.output_dir = '/mnt/ilustre/users/sanger-dev/workspace/20190506/Bacgenome_tsg_34115/output'
        sample = 'U00096.3'

        # if os.path.exists(self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + ".stat.xls"):
        #     self.prephage.add_prephage_stat(prephage_id,
        #                                         self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + ".stat.xls",
        #                                         sample)
        # if os.path.exists(
        #                                                 self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + "_prephage_summary.xls"):
        #     self.prephage.add_prephage_detail(prephage_id,
        #                                       self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + "_prephage_summary.xls",
        #                                       sample)
        # if os.path.exists(
        #                                                 self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + "_prephage_detail.xls"):
        #     self.prephage.add_prephage_gene(prephage_id,
        #                                     self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + "_prephage_detail.xls",
        #                                     sample)

        self.two_component_api= self.api.api("bacgenome.common_api")
        two_component_id = self.two_component_api.add_main('anno_regulator', name='two_component')
        mongo_key1 = 'gene_id,type,pfam_id,domain,domain_desc,,,,,,,location,'
        two_com1 = self.output_dir + '/'+sample + '/annotation/Two_component/'+ sample+'.senser_regulator.xls'
        two_com2 = self.output_dir + '/' + sample + '/annotation/Two_component/'+sample+ '.senser_regulator.stat'
        self.two_component_api.add_main_detail(two_com1,'anno_regulator_detail', two_component_id, mongo_key1, has_head =True,
                                               main_name='regulator_id',other_dic={'specimen_id':sample})
        mongo_key2 = 'senser,regulator,hybrid'
        self.two_component_api.add_main_detail(two_com2,'anno_regulator_stat', two_component_id, mongo_key2, has_head =True,
                                               main_name='regulator_id',other_dic={'specimen_id':sample})


        self.end()



    def end(self):
        super(TestApiWorkflow, self).end()

if __name__ == "__main__":

    from biocluster.wsheet import Sheet

    data = {
        'name': 'api_export',
        'id': 'tsg_34115',
        'type': 'workflow',
        'options': {
            "test" : 't'
            }
        }

    wsheet = Sheet(data = data)
    wf = TestApiWorkflow(wsheet)
    wf.run()

