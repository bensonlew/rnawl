# -*- coding: utf-8 -*-
# __author__ = 'shijin, qinjincheng'

from mainapp.controllers.project.small_rna_controller import SmallRnaController
from mainapp.libs.signature import check_sig
import web
import json
import types
from bson.objectid import ObjectId

class GenesetDeleteAction(SmallRnaController):
    def __init__(self):
        super(GenesetDeleteAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        input_data = web.input()

        ret = self.check_params(input_data)
        if ret != True:
            return ret

        ret = self.delete_geneset(input_data.geneset_id, input_data.task_id)
        if ret == True:
            task_info = {'success': True, 'info': 'succeed in deleting geneset'}
            return json.dumps(task_info)
        else:
            task_info = {'success': True, 'info': 'fail to delete geneset'}
            return json.dumps(task_info)

    def check_params(self, data):
        expected_args = ['task_id', 'geneset_id']
        for arg in expected_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Lack argument: {}'.format(arg)}
                return json.dumps(info)
        return True

    def delete_geneset(self, geneset_id, task_id):
        if isinstance(geneset_id, types.StringTypes):
            geneset_id = ObjectId(geneset_id)
        elif isinstance(geneset_id, ObjectId):
            pass
        else:
            raise Exception('geneset_id must be str or ObjectId!')
        sg_geneset_info = self.db['sg_geneset_info']
        sg_status = self.db['sg_status']
        results = sg_geneset_info.find({'geneset_id': geneset_id})
        for result in results:
            col_name = result['col_name']
            col_id = result['col_id']
            col = self.db[col_name]
            try:
                col_result = col.find_one({'_id': col_id})
                if col_result['task_id'] == task_id:
                    col_result['params'] = ''
                    col.update({'_id': col_id}, {'$set': col_result}, upsert=True)
            except:
                print 'can not find document of which _id is {} in {}'.format(col_id, col_name)
            try:
                col_result = sg_status.find_one({'table_id': col_id})
                col_result['status'] = 'deleted'
                sg_status.update({'table_id': col_id}, {'$set': col_result})
            except:
                print 'can not find document of which table_id is {} in sg_status'.format(col_id)

        sg_geneset = self.db['sg_geneset']
        sg_geneset_detail = self.db['sg_geneset_detail']
        result = sg_geneset.find_one({'main_id': geneset_id, 'task_id': task_id})
        if result:
            sg_geneset_detail.delete_many({'geneset_id': result['_id']})
            sg_geneset.remove({'main_id': geneset_id, 'task_id': task_id})
        return True
