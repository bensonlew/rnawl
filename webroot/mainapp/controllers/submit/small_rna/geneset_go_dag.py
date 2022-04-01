# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

from mainapp.controllers.project.small_rna_controller import SmallRnaController
from mainapp.libs.signature import check_sig
import web
import json
import datetime
import unittest
import os
from bson.objectid import ObjectId

class GenesetGoDagAction(SmallRnaController):
    def __init__(self):
        super(GenesetGoDagAction,self).__init__(instant=False)

    @check_sig
    def POST(self):
        ret = self.check_params()
        if ret is not True:
            return ret
        self.pack_params()

        ret = self.check_main_table()
        if ret:
            self.delete_main_table()
        else:
            self.create_main_table()

        self.prepare_workflow()
        task_info = super(GenesetGoDagAction, self).POST()
        task_info['content'] = {'ids': {'id': str(self.main_id), 'name': self.name}}

        self.small_rna.insert_geneset_info(self.geneset_id, 'sg_geneset_go_dag', str(self.main_id))

        return json.dumps(task_info)

    def check_params(self):
        self.input_data = web.input()
        self.expected_args = ['task_id', 'submit_location', 'task_type']
        # go_enrich_id ~ sg_geneset_go_enrich.main_id
        self.expected_args.extend(['go_enrich_id'])
        # 获取geneset_id
        go_enrich_id = self.input_data.go_enrich_id
        result = self.db['sg_geneset_go_enrich'].find_one({'main_id': ObjectId(go_enrich_id)})
        if result:
            params = json.loads(result['params'])
            self.geneset_id = params['geneset_id']
        if hasattr(self.input_data, 'go_list'):
            # go_list ~ sg_annotation_go_detail.goid_2 (multiple inputs would be joined by comma)
            self.expected_args.extend(['go_list'])
        elif hasattr(self.input_data, 'significant_diff') and hasattr(self.input_data, 'significant_value') and hasattr(self.input_data, 'top_num'):
            # significant_diff in ['pvalue', 'padjust']
            # significant_value ~ 0.05
            # top_num ~ 10
            self.expected_args.extend(['significant_diff', 'significant_value', 'top_num'])
        else:
            info = {'success': False, 'info': 'Must specify argument in significant_diff&significant_value&top_num or go_list'}
            return json.dumps(info)
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                info = {'success': False, 'info': 'Lack argument: {}'.format(arg)}
                return json.dumps(info)
        if hasattr(self.input_data, 'significant_value'):
            if float(self.input_data['significant_value']) <= 0 or float(self.input_data['significant_value']) > 1:
                info = {'success': False, 'info': 'The value of Pvalue or Padjust must in (0, 1]'}
                return json.dumps(info)
        if hasattr(self.input_data, 'top_num'):
            if int(self.input_data['top_num']) <= 0:
                info = {'success': False, 'info': 'The number of GO term must greater than 0'}
                return json.dumps(info)
        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        for each in self.expected_args:
            if each == 'task_type':
                params_dict[each] = int(input_data_dict[each])
            else:
                params_dict[each] = input_data_dict[each]
        self.params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

    def check_main_table(self):
        sg_geneset_go_dag = self.db['sg_geneset_go_dag']
        result = sg_geneset_go_dag.find_one({'task_id': self.input_data.task_id})
        if result:
            return True
        else:
            return False

    def delete_main_table(self):
        sg_geneset_go_dag = self.db['sg_geneset_go_dag']
        sg_geneset_go_dag.remove({'task_id': self.input_data.task_id})
        self.create_main_table()

    def create_main_table(self):
        result_info = self.small_rna.get_task_info(self.input_data.task_id)
        self.project_sn = result_info['project_sn']
        time_now = datetime.datetime.now()
        self.name = 'GenesetGoDag_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'Geneset go Dag main table built at controller'
        main_info = {
            # essential keys
            'task_id': self.input_data.task_id,
            'project_sn': self.project_sn,
            'created_ts': created_ts,
            'name': self.name,
            'desc': desc,
            'params': self.params,
            # status of process
            'status': 'start',
        }
        self.main_id = self.small_rna.insert_main_table('sg_geneset_go_dag', main_info)

    def prepare_workflow(self):
        input_data_dict = dict(self.input_data)
        main_info = self.small_rna.get_main_info(self.input_data.go_enrich_id, 'sg_geneset_go_enrich', self.input_data.task_id)
        options = dict()
        for each in self.expected_args:
            if input_data_dict[each] == '':
                options[each] = None
            else:
                options[each] = input_data_dict[each]
        options.update({
            'go_enrich_id': str(self.main_id),
            'update_info': json.dumps({str(self.main_id): 'sg_geneset_go_dag'}),
            'go_enrich_detail': main_info['result_file']
        })
        self.set_sheet_data(name='small_rna.report.geneset_go_dag',
                            options=options,
                            main_table_name=self.name,
                            module_type='workflow',
                            project_sn=self.project_sn,
                            task_id=self.input_data.task_id)

class TestFunction(unittest.TestCase):
    '''
    This is test for the controller. Just run this script to do test.
    '''
    def test_default(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += '-fr no '
        cmd += '-c {} '.format('client03')
        cmd += 's/small_rna/geneset_go_dag '
        cmd += '-b http://192.168.12.102:9090 '
        args = {
            'task_id': 'small_rna',
            'submit_location': 'genesetgodag',
            'task_type': '2',
            'go_enrich_id': '5bf390c0a4e1af62168452ff',
            'significant_diff': 'pvalue',
            'significant_value': '0.05',
            'top_num': '10',
        }
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)

    def test_go_list(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += '-fr no '
        cmd += '-c {} '.format('client03')
        cmd += 's/small_rna/geneset_go_dag '
        cmd += '-b http://192.168.12.102:9090 '
        args = {
            'task_id': 'small_rna',
            'submit_location': 'genesetgodag',
            'task_type': '2',
            'go_enrich_id': '5bf390c0a4e1af62168452ff',
            'go_list': 'GO:0044464;GO:0044424;GO:0005488;GO:0009987;GO:0065007;GO:0043226;GO:0044699;GO:0043227;GO:0043229;GO:0044763;GO:0043231;GO:0044444;GO:0044422;GO:0044446;GO:0005515;GO:0019222;GO:0043167;GO:0031323;GO:0060255;GO:0080090',
        }
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()