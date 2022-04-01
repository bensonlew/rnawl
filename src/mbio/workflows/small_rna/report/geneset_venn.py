# -*- coding: utf-8 -*-
# __author__ = "konghualei, qinjincheng"

from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
import datetime
import shutil
import re,os, time
from biocluster.workflow import Workflow
import time

class GenesetVennWorkflow(Workflow):
    """
    Update a document in  sg_geneset_venn by specific main_id with inserting {'workflow': 'used'}
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.start_listener()
        self.fire("start")
        time.sleep(3)
        self.set_db()

    def end(self):
        super(GenesetVennWorkflow, self).end()

    def set_db(self):
        inst = self.api.api("small_rna.geneset_venn")
        inst.insert_key_value_pair(main_id=self.option("main_id"))
        self.end()
