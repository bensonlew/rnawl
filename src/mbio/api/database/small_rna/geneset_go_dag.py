# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.small_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
import sys
import types
from bson.objectid import ObjectId

class GenesetGoDag(ApiBase):
    def __init__(self, bind_object):
        super(GenesetGoDag, self).__init__(bind_object)

    @report_check
    def update_geneset_go_dag(self, go_enrich_id, output_dir):
        self.bind_object.logger.debug('incoming go_enrich_id in {} is {}'.format(sys._getframe().f_code.co_name,
                                                                                 go_enrich_id))
        if go_enrich_id == None:
            self.bind_object.set_error('incoming go_enrich_id is None, abord')
        elif isinstance(go_enrich_id, types.StringTypes):
            go_enrich_id = ObjectId(go_enrich_id)
        elif isinstance(go_enrich_id, ObjectId):
            go_enrich_id = go_enrich_id

        self.update_db_record(table_name='sg_geneset_go_dag',
                              record_id=go_enrich_id,
                              output_dir=output_dir,
                              main_id=go_enrich_id,
                              status='end')
