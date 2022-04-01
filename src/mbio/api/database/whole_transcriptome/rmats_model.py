# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.whole_transcriptome.api_base import ApiBase
from biocluster.api.database.base import report_check
import os
from bson.objectid import ObjectId
import unittest

class RmatsModel(ApiBase):
    def __init__(self, bind_object):
        super(RmatsModel, self).__init__(bind_object)

    @report_check
    def add_rmats_model(self, output_dir, s3_output, main_id):
        insert_dict = {
            'pdf_files': [i for i in os.listdir(output_dir) if i.endswith('.pdf')],
            'png_files': [i for i in os.listdir(output_dir) if i.endswith('.png')],
            'graph_dir': s3_output,
            'main_id': ObjectId(main_id),
            'status': 'end'
        }
        self.update_db_record('splicing_rmats_model', main_id, insert_dict=insert_dict)
