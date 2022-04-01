# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from mainapp.controllers.project.small_rna_controller import SmallRnaController
from mainapp.libs.signature import check_sig
import web
import json
import datetime
import unittest
import os
from biocluster.config import Config

class CircosRedrawnAction(SmallRnaController):
    def __init__(self):
        super(CircosRedrawnAction,self).__init__(instant=False)

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        self.pack_params()
        self.create_main_table()
        self.prepare_workflow()
        task_info = super(CircosRedrawnAction, self).POST()
        task_info['content'] = {'ids': {'id': str(self.main_id), 'name': self.name}}
        return json.dumps(task_info)

    def check_params(self):
        self.input_data = web.input()
        self.expected_args = ['task_id', 'submit_location', 'task_type']
        self.expected_args.extend(["chr_num", 'group_id', 'group_dict', 'samples'])
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                info = {'success': False, 'info': 'Lack argument: {}'.format(arg)}
                return json.dumps(info)

        origin_info = self.small_rna.get_main_info_by_record('sg_circos', task_id=self.input_data.task_id, type="origin")
        if origin_info and origin_info['status'] == "start":
            info = {'success': False, 'info': '正在运行， 稍后再试'}
            return json.dumps(info)
        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        for each in self.expected_args:
            if each == 'group_dict':
                params_dict[each] = json.loads(input_data_dict[each])
            elif each == 'task_type':
                params_dict[each] = int(input_data_dict[each])
            else:
                params_dict[each] = input_data_dict[each]
        self.params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

    def create_main_table(self):
        origin_circos = self.small_rna.get_main_info_by_record('sg_circos', task_id=self.input_data.task_id)
        if not origin_circos:
            origin_info = self.small_rna.get_main_info_by_record('sg_mapping', task_id=self.input_data.task_id, type="origin")
        else:
            origin_info = origin_circos

        if hasattr(self.input_data, "samples"):
            sample_list = self.input_data.samples.split(",")
        else:
            sample_list = origin_info['sample_list'],
        result_info = self.small_rna.get_task_info(self.input_data.task_id)

        self.project_sn = result_info['project_sn']
        time_now = datetime.datetime.now()
        self.name = 'SmallRnaCircos_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = '基因组分布圈图'
        main_info = {
            # essential keys
            'task_id': self.input_data.task_id,
            'project_sn': self.project_sn,
            'created_ts': created_ts,

            'name': self.name,
            'desc': desc,
            'params': self.params,
            'sample_list': origin_info['sample_list'],
            'sample_list_circos': sample_list,
            'type': "rerun",
            'status': 'start',
        }
        if not origin_circos:
            self.main_id = self.small_rna.insert_main_table('sg_circos', main_info)
        else:
            self.small_rna.delete_main_table('sg_circos', origin_circos['_id'])
            self.main_id = self.small_rna.insert_main_table('sg_circos', main_info)


    def prepare_workflow(self):
        input_data_dict = dict(self.input_data)
        task_info_dict = self.small_rna.get_task_info(self.input_data.task_id)
        if 'genome_id' in task_info_dict:
            genome_info = self.small_rna.get_genome_info_new(task_info_dict['genome_id'])
        else:
            genome_info = self.small_rna.get_genome_info(task_info_dict['organism_name'], task_info_dict['assembly'], task_info_dict['annot_version'])
        origin_info = self.small_rna.get_main_info_by_record('sg_mapping', task_id=self.input_data.task_id, type="origin")
        genome_dir = Config().SOFTWARE_DIR + "/database/Genome_DB_finish/"

        if "result_dir" in origin_info:
            origin_info['result_origin'] = origin_info["result_dir"]

        if origin_info['result_origin'].endswith("/"):
            pass
        else:
            origin_info['result_origin'] = origin_info['result_origin'] + "/"
        options = dict()
        for each in self.expected_args:
            if input_data_dict[each] == '':
                options[each] = None
            else:
                options[each] = input_data_dict[each]
        if 'config' in task_info_dict and os.path.exists(task_info_dict['config']):
            config = task_info_dict['config']
        else:
            config = origin_info['result_origin'] + "qc_file.config"
        options.update({
            'arf': origin_info['result_origin'] + "reads_vs_genome.arf",
            'map_stat': origin_info['result_origin'],
            'gtf': genome_dir + genome_info['gtf'],
            'ref': genome_dir + genome_info['dna_fa'],
            'config': config,
            'main_id': str(self.main_id),
            'chr_num': self.input_data.chr_num,
            'update_info': json.dumps({str(self.main_id): 'sg_circos'}),
        })
        if hasattr(self.input_data, "samples"):
            options.update({
                "samples": self.input_data.samples
            })
        self.set_sheet_data(name='small_rna.report.circos_redrawn',
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
        cmd += 's/small_rna/circos_redrawn '
        cmd += '-b http://192.168.12.102:9090 '
        args = {
            'task_id': 'small_rna',
            'submit_location': 'mapping_circs',
            'task_type': '2',
            'chr_num': '10',
        }
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
