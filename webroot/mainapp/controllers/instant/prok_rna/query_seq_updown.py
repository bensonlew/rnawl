# -*- coding: utf-8 -*-
import os
import web
import json
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mainapp.libs.signature import check_sig
import unittest
import datetime
import re
# __author__ = 'fengyitong'


class QuerySeqUpdownAction(ProkRNAController):
    def __init__(self):
        super(QuerySeqUpdownAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        # check
        for arg in ['upstream', 'downstream', 'task_id', 'seq_id']:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)

        sg_task = self.prok_rna.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        ref_genome = sg_task['ref_genome']

        rock_info = sg_task['rock_index']

        if re.match(r'^\w+://\S+/.+$', rock_info):
            inter_dir = self.create_tmp_dir(data.task_id, "rock_info/")
            rock_info_fna = self.download_from_s3(os.path.join(rock_info, 'genome_fna.db.sqlite3'), inter_dir=inter_dir)
            rock_info_bed = self.download_from_s3(os.path.join(rock_info, 'ptt.bed'), inter_dir=inter_dir)
            rock_info = rock_info_bed.split('ptt.bed')[0]

        ptt_path = os.path.join(rock_info, 'ptt.bed')
        genome_path = os.path.join(rock_info, 'genome_fna.db.sqlite3')

        # save
        mongo_data = dict([('task_id', data.task_id),
                           ('project_sn', project_sn),
                           ('upstream', data.upstream),
                           ('downstream', data.downstream),
                           ('seq_id', data.seq_id),
                           ])

        collection_name = 'sg_query_seq_down'
        main_table_id = self.prok_rna.insert_or_replace_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            # "submit_location": data.submit_location,
            # 'update_info': json.dumps(update_info),
            "main_id": str(main_table_id),
            "seq_id": data.seq_id,
            "upstream": int(data.upstream),
            "downstream": int(data.downstream),
            "ref_genome": ref_genome,
            "ptt_path": ptt_path,
            "genome_path": genome_path
            }
        task_name = 'prok_rna.report.query_seq_updown'
        main_table_name = 'query_seq_updown' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=[], project_sn=project_sn, task_id=data.task_id)
        task_info = super(QuerySeqUpdownAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }}
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/prok_rna/query_seq_updown "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna_workflow",
            seq_id="YE4208",
            # seq_db_path="/mnt/ilustre/users/sanger-dev/workspace/20180831/Single_Srna_8837_fyt/Srna/output/rockhopper",
            upstream='40',
            downstream='50'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
