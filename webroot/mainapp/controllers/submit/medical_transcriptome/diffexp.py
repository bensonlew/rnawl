# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
import os
import unittest
from bson.objectid import ObjectId


class DiffexpAction(MedicalTranscriptomeController):
    def __init__(self):
        super(DiffexpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'control_id',
                       'exp_id', 'level', 'diff_method',
                       'fc', 'kind', 'is_batch']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:   %s" % arg, 'code': 'C2900301', 'variables': variables}
                return json.dumps(info)
        # if not hasattr(data, 'tpm_filter_threshold'):
        #     data.tpm_filter_threshold = 0
        # elif not data.tpm_filter_threshold:
        #     data.tpm_filter_threshold = 0
        exp_info = self.medical_transcriptome.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        group_size = list()
        success = list()
        for key in group_dict:
            group_size.append(len(group_dict[key]))
        group_size.sort()
        if group_size[0] == group_size[-1] == 1:
            if data.diff_method == "DESeq2":
                success.append('If your experiment has no biological replicates, please select DEGseq, edgeR, limma, NOIseq')
        elif group_size[0] == 1 and group_size[1] >= 2:
            if data.diff_method == "DESeq2":
                success.append('If some groups have biological replicates and some groups no, please select edgeR')
        # elif group_size[0] >= 2:
        #     if data.diff_method == "DEGseq":
        #         success.append('If your experiment has biological replicates, please select  DESeq2, edgeR, limma, NOIseq.')
        if success:
            variables = []
            variables.append(success[0])
            info = {'success': False, 'info': '%s' % success[0], 'code': 'C2900302', 'variables': variables}
            return json.dumps(info)

        # create main table record
        level = exp_info['level']
        exp_type, quant_method = exp_info['exp_type'], exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            control_id=data.control_id,
            level=level,
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            fc=str(float(data.fc)),
            filter_method = data.filter_method,
            # stat_type=data.stat_type,
            # stat_cutoff=data.stat_cutoff,
            # quant_method=quant_method,
            diff_method=data.diff_method,
            kind=data.kind,
            # tpm_filter_threshold=str(float(data.tpm_filter_threshold)),
            is_batch=data.is_batch
        )
        # for special args
        if hasattr(data, 'tpm_filter_threshold'):
            tpm_filter_threshold = data.tpm_filter_threshold
            params.update(dict(tpm_filter_threshold=tpm_filter_threshold))
        if hasattr(data, 'stat_type'):
            stat_type = data.stat_type
            params.update(dict(stat_type=stat_type))
        if hasattr(data, 'stat_cutoff'):
            stat_cutoff = data.stat_cutoff
            params.update(dict(stat_cutoff=stat_cutoff))
        if hasattr(data, 'correct_method'):
            correct_method = data.correct_method
            params.update(dict(correct_method=correct_method))
        else:
            correct_method = "BH"
        # if hasattr(data, 'is_batch') and data.is_batch == 'True':
        #     params.update({
        #         "is_batch": data.is_batch
        #     })
        if hasattr(data,'has_batch'):
            params.update({
                "has_batch": data.has_batch
            })
        if hasattr(data, "batch_matrix"):
            params.update({
                "batch_matrix": data.batch_matrix
            })
        if hasattr(data, "prob"):
            params.update({
                "prob": data.prob
            })
        if hasattr(data, 'file_id'):
            params.update({
                'file_id': data.file_id
            })
        # if hasattr(data, "batch_method"):
        #     params.update({
        #         "batch_method": data.batch_method
        #     })
        # if not hasattr(data, 'is_batch'):
        #     data.is_batch = False
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        connect = self.db['sg_exp']
        record = connect.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id)})
        is_rmbe = str(record['is_rmbe']).lower()
        if is_rmbe == 'false':
            exp_id = data.exp_id
        if is_rmbe == 'true':
            exp_id = str(record['batch_main_id'])
        if data.is_batch == 'False':
            name = "DiffExpress" + '_' + level + '_' + quant_method + '_' + data.diff_method + '_'
        else:
            name = "DiffExpress_batch" + '_' + level + '_' + quant_method + '_' + data.diff_method + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        last_record = self.db["sg_diff"].find_one({"task_id": task_id}, sort=[('_id', -1)])
        if last_record:
            last_flag = last_record["flag"]
            current_flag = last_record["flag"] + 1
        else:
            current_flag = 1
        print(str(current_flag))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            version="v1",
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='differential analysis main table',
            level=level,
            group_dict=group_dict,
            group_id=ObjectId(data.group_id),
            exp_id=ObjectId(data.exp_id),
            params=params,
            is_batch=data.is_batch,
            status="start",
            flag =current_flag
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_diff', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        if not hasattr(data, 'tpm_filter_threshold'):
            tpm_filter_threshold = 0
        elif not data.tpm_filter_threshold:
            tpm_filter_threshold = 0

        # prepare option for workflow
        options = {
            "exp_matrix": "{};{};{}".format(exp_id, level, is_rmbe),
            "group_dict": json.dumps(group_dict),
            "diff_main_id": str(main_id),
            "method": data.diff_method,
            "count": self.medical_transcriptome.get_count_path(data.exp_id),
            "fc": float(data.fc),
            "group": json.dumps(group_dict),
            "cmp": data.control_id,
            # "pvalue": float(data.stat_cutoff),
            # "pvalue_padjust": data.stat_type,
            # "padjust_way": correct_method,
            "update_info": json.dumps({str(main_id): "sg_diff"}),
            "type": data.kind,
            "task_id": data.task_id,
            'filter_method':data.filter_method,
            "exp_level": level,
            "exp_type": exp_type,
            "tpm_filter_threshold":tpm_filter_threshold,
            'is_batch': data.is_batch,
            'main_table_data': main_table_data,
        }


        if hasattr(data, 'stat_type'):
            options.update({
                'pvalue_padjust': data.stat_type
            })
        if hasattr(data, 'stat_cutoff'):
            options.update({
                'pvalue':float(data.stat_cutoff)
            })
        if hasattr(data, 'correct_method'):
            options.update({
                'padjust_way': data.correct_method
            })
        if hasattr(data, 'has_batch'):
            options.update({
                "has_batch": data.has_batch
            })
        if hasattr(data, "batch_matrix"):
            options.update({
                "batch_matrix": data.batch_matrix
            })
        if hasattr(data, 'prob'):
            options.update({
                'prob': data.prob
            })
        # if hasattr(data, "batch_method"):
        #     options.update({
        #         "batch_method": data.batch_method
        #     })
        if hasattr(data, 'is_batch') and data.is_batch == 'True' and data.has_batch == 'False':
            options['method'] = 'svaseqlimma'
        # prepare to file
        to_files = ["medical_transcriptome.export_exp_matrix_new(exp_matrix)",
                    "medical_transcriptome.export_group(group)",
                    "medical_transcriptome.export_compare(cmp)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.diffexp_batch'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(DiffexpAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
        # task_info['group_dict'] = group_dict
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
        cmd += "s/medical_transcriptome/diffexp_batch "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="diffexp",
            exp_id='5f8e8ec117b2bf3a1b24bb5c',
            group_id='5f8e8da017b2bf3a1b1dd32a',
            exp_level='T',
            group_dict=r'{"H1": ["H1581_1", "H1581_2", "H1581_3"],"H2": ["H1581_4", "H1581_5", "H1581_6"],"H3": ["H1581_7", "H1581_8", "H1581_9"],"S1": ["SNU16_1", "SNU16_2", "SNU16_3"],"S2": ["SNU16_4", "SNU16_5", "SNU16_6"],"S3": ["SNU16_7", "SNU16_8", "SNU16_9"]}'.replace('"', '\\"'),
            type='all',
            control_id="5f8e8da017b2bf3a1b1dd32b",
            diff_method='DESeq2',
            stat_type='padjust',
            stat_cutoff='0.05',
            correct_method="BH",
            # prob='0.8',
            fc='2',
            tpm_filter_threshold='0',
            is_batch='False',
            # is_rmbe="false"
            # has_batch='True',
            # batch_matrix='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/batch.txt'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
