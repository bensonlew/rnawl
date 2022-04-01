# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from mbio.api.to_file.itraq_tmt import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ProteinsetClusterAction(ItraqTmtController):
    def __init__(self):
        super(ProteinsetClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'use_group',
                       'proteinset_id', 'sct', 'gct',
                       ]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: %s", "variables":[arg], "code" : "C1900901"}
                return json.dumps(info)
        # if data.exp_level != data.proteinset_type:
        #     info = "exp_level is not consistent with proteinset_type!"
        #     return json.dumps(info)
        if data.group_id.lower() == 'all' and data.use_group.lower() != 'no':
            info = {'success': False, 'info': "When group is all, group samples calculation cannot be mean!", "code" : "C1900902"}
            return json.dumps(info)
        exp_info = self.itraq_tmt.get_exp_params_info(data.task_id)

        # print  exp_info
        express_id = str(exp_info['main_id'])
        # print data.task_id
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)

        # print "group_num is {}".format(len(group_dict.keys()))
        # print "group_us is {}".format(data.use_group.lower())
        if len(group_dict.keys()) < 2 and data.use_group.lower() != "no":
            info = {'success': False, 'info': "Using grouped samples calculation, two or more groups should be selected!", "code" : "C1900903"}
            # print "info is {}".format(info)
            # print "return info is {}".format(json.dumps(info))
            return json.dumps(info)


        # create main table record
        # exp_level = exp_info['exp_level']
        # exp_type, quant_method = exp_info['exp_type'], exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            express_id=str(exp_info['main_id']),
            group_id=data.group_id,
            group_dict=group_dict,
            # quant_method=quant_method,
            use_group=data.use_group,
            proteinset_id=data.proteinset_id,
            sct=data.sct,
            gct=data.gct,
        )

        n_clusters = 10
        scm = "complete"
        scd = "correlation"
        gcm = "average"
        gcd = "euclidean"
        if data.sct == 'no':
            pass
        elif data.sct == 'kmeans':
            params.update(dict(scd=data.scd, n_clusters=int(data.n_clusters), ))
            scd = data.scd
            n_clusters = int(data.n_clusters)
        else:
            params.update(dict(scm=data.scm, scd=data.scd))
            scm = data.scm
            scd = data.scd
        if data.gct.startswith('h'):
            params.update(dict(gcm=data.gcm, gcd=data.gcd, n_clusters=int(data.n_clusters)))
            gcm = data.gcm
            gcd = data.gcd
            n_clusters = int(data.n_clusters)
        elif data.gct.startswith('k'):
            params.update(dict(gcd=data.gcd, n_clusters=int(data.n_clusters)))
            gcd = data.gcd
            n_clusters = int(data.n_clusters)
        else:
            pass

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Pset_Cluster" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='proteinset cluster main table',
            params=params,
            status="start"
        )
        main_id = self.itraq_tmt.insert_main_table('sg_proteinset_cluster', main_info)

        print("++++ exp_id protein_set_id")
        print(express_id + ";" + data.proteinset_id)
        # prepare option for workflow
        options = dict(
            exp_matrix=express_id + ";" + data.proteinset_id,
            group_dict=data.group_dict,
            group_id=data.group_id,
            group=data.group_dict,
            cluster_main_id=str(main_id),
            n_clusters=n_clusters,
            use_group=data.use_group,
            scm=scm,
            scd=scd,
            sct=data.sct,
            gcm=gcm,
            gct=data.gct,
            gcd=gcd,
            update_info=json.dumps({str(main_id): "sg_proteinset_cluster"})  # to update sg_status
        )

        # prepare to file
        to_files = ["itraq_tmt.export_proteinset_exp_matrix(exp_matrix)",
                    "itraq_tmt.export_group(group)",]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'itraq_and_tmt.report.proteinset_cluster'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ProteinsetClusterAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # 更新基因集的使用信息
        self.itraq_tmt.insert_proteinset_info(data.proteinset_id, name, str(main_id))
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
        cmd += "s/itraq_and_tmt/proteinset_cluster "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="itraq_tmt",
            task_type="2",
            submit_location="proteinsetcluster",
            group_id="59f3975f28fb4f19526ebee8",
            group_dict=json.dumps({"control_0":["control_0_1","control_0_2","control_0_3"],"BR_0":[ "BR_0_1","BR_0_2","BR_0_3" ], "control_3":["control_3_1","control_3_2","control_3_3"], "BR_3":["BR_3_1","BR_3_2","BR_3_3"]}).replace('"', '\\"'),
            express_id="5aaf2760a4e1af0f99070f45",
            proteinset_id="5a9de1715b8c915efb2a3293",
            quant_method='RSEM',
            use_group='no',
            n_clusters='10',
            sct="hierarchy",
            scm="complete",
            scd="correlation",
            gct="hierarchy",
            gcm="average",
            gcd="euclidean",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
