# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class CircosAction(BacgenomeController):
    def __init__(self):
        super(CircosAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'specimen_id', 'submit_location', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bacgenome.report.circos'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        if not hasattr(data, "seq_type"):
            seq_type = "Circular"
        else:
            seq_type = data.seq_type
        params = {
            'task_id': data.task_id,
            'specimen_id': data.specimen_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'seq_type':seq_type
        }
        lab = ''

        if hasattr(data, 'para_def'):  #guanqing 20190410
            params['para_def'] = data.para_def   #'Scaffold1,200,400;Scaffold2,600,700'
            lab += ';para_def'

        if hasattr(data, 'para1'):
            params['para1'] = data.para1
            lab += ';' + data.para1

        if hasattr(data, 'para2'):
            params['para2'] = data.para2
            lab += ';' + data.para2
        if hasattr(data, 'para3'):
            params['para3'] = data.para3
            lab += ';' + data.para3
        if hasattr(data, 'para4'):
            params['para4'] = data.para4
            lab += ';' + data.para4
        if hasattr(data, 'para5'):
            params['para5'] = data.para5
            lab += ';' + data.para5
        if lab !='':
            lab = lab[1:]
        if hasattr(data, 'location'):
            params['location'] = data.location
            pre_file,color_version = self.bacgenome.get_pre_file(data.task_id, data.specimen_id, data.location,return_color_version=True)
        else:
            pre_file,color_version = self.bacgenome.get_pre_file(data.task_id, data.specimen_id, return_color_version=True)
        if re.search(r'/$',pre_file):
            pre_file =pre_file
        else:
            pre_file = pre_file + '/'
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        if hasattr(data, 'location'):
            main_table_name = data.specimen_id + '_' + 'Circos_' + data.location + '_' + datetime.datetime.now().strftime(
                "%Y%m%d_%H%M%S%f")[:-3]
        else:
            main_table_name = data.specimen_id + '_' + 'Circos_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[
                                                                   :-3]
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('version','3.0'),   # by zzg 从params移出来
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params),
            ('lab', lab.replace(';',','))
        ]
        main_id = self.bacgenome.insert_main_table('circos', mongo_data)
        update_info[str(main_id)] = 'circos'
        options = {'task_id': data.task_id,
                   'specimen_id': data.specimen_id,
                   # 'old_species':data.specimen_id,
                   'update_info': json.dumps(update_info),
                   'params': params,
                   'main_id': str(main_id),
                   "pre_file": pre_file,
                   'main_table_data': SON(mongo_data),
                   'labs':lab,
                   'seq_type':seq_type,
                   'color_version' : str(color_version)
                   }
        to_file = []

        if hasattr(data, 'para_def'):   #zouguanqing 20190410
            options['para_def'] = data.para_def  #'Scaffold1,200,400;Scaffold2,600,700'
            to_file.append('bacgenome.export_def_circle(para_def)')
        if hasattr(data, 'para1'):
            options['para1'] = data.para1
            to_file.append('bacgenome.export_prephage_circos(para1)')
        if hasattr(data, 'para2'):
            options['para2'] = data.para2
            to_file.append('bacgenome.export_prephage_circos(para2)')
        if hasattr(data, 'para3'):
            options['para3'] = data.para3
            to_file.append('bacgenome.export_prephage_circos(para3)')
        if hasattr(data, 'para4'):
            options['para4'] = data.para4
            to_file.append('bacgenome.export_prephage_circos(para4)')
        if hasattr(data, 'para5'):
            options['para5'] = data.para5
            to_file.append('bacgenome.export_prephage_circos(para5)')
        if hasattr(data, 'location'):
            options['location'] = data.location
        else:
            options['location'] = "Scaffold"
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="Circos/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params,
                            to_file=to_file)
        task_info = super(CircosAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)

