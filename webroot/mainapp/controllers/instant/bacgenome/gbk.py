#- coding: utf-8 -*-
# __author__ = 'gaohao'
import web
import json,re
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class GbkAction(BacgenomeController):
    def __init__(self):
        super(GbkAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['gbk_id', 'specimen_id', 'gbk_detail_id', 'submit_location', 'task_type','task_id','gbk_path','client']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.gbk'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        params = self.bacgenome.get_gbk(ObjectId(data.gbk_detail_id),ObjectId(data.gbk_id))
        path = self.bacgenome.get_gbk_path(ObjectId(data.gbk_detail_id))
        path1 = ''
        if re.search(r'://',path):
            if data.gbk_path in ['Scaffold', 'scaffold']:
                target_file =path + data.specimen_id + '.gbk'
            else:
                target_file = path + data.specimen_id + '_' + data.gbk_path + '.gbk'
        else:
            if data.client == 'client03':
                path1 = '/mnt/ilustre/tsanger-data'
            if data.client == 'client01':
                path1 = '/mnt/ilustre/data'
            if data.gbk_path in ['Scaffold', 'scaffold']:
                target_file = path1 + "/rerewrweset/" + path + '/' + data.specimen_id + '.gbk'
            else:
                target_file = path1 + "/rerewrweset/" + path + '/' + data.specimen_id + '_' + data.gbk_path + '.gbk'
        main_table_name = data.specimen_id + '_' + 'Gbk_' + data.gbk_path + '_' + datetime.datetime.now().strftime(
            "%Y%m%d_%H%M%S%f")[:-3]
        options = {'des': params['des'],
                   'source': params['source'],
                   'organism': params['organism'],
                   'author': params['author'],
                   'title': params['title'],
                   'journal': params['journal'],
                   'gbk_detail_id': data.gbk_detail_id,
                   'gbk_id': data.gbk_id,
                   'task_id':data.task_id,
                   'specimen_id':data.specimen_id,
                   'seq_type':data.gbk_path,
                   'gbk_file':target_file,
                   }
        self.set_sheet_data(name=task_name,
                            options=options,
                            module_type=module_type,
                            main_table_name="Gbk/" + main_table_name,
                            project_sn=project_sn,
                            task_id=data.task_id,params='analysis:gbk')
        task_info = super(GbkAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(data.gbk_id), 'name': main_table_name}}
        return json.dumps(task_info)