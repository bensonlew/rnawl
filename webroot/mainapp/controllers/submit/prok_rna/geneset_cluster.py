# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os
from bson.objectid import ObjectId


class GenesetClusterAction(ProkRNAController):
    def __init__(self):
        super(GenesetClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'use_group',
                       'geneset_id', 'sct', 'gct',
                       ]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        if hasattr(data, 'submit_type'):
            submit_type = int(data.submit_type)
        else:
            submit_type = 0
        # if data.exp_level != data.geneset_type:
        #     info = "exp_level is not consistent with geneset_type!"
        #     return json.dumps(info)
        if data.group_id.lower() == 'all' and data.use_group.lower() != 'no':
            info = {'success': False, 'info': "group 设置为all时,不能选择以分组均值分析"}
            return json.dumps(info)

        # print "+++" + data.task_id

        #通过main_id和task_id找到sg_exp表
        exp_info = self.prok_rna.get_exp_params_info(data.exp_id, data.task_id)

        # print  exp_info
        express_id = str(exp_info['main_id'])
        # print data.task_id
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        try:
            group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        except:
            info = {'success': False, 'info': "group_dict参数不正确"}
            return json.dumps(info)

        if data.group_id.lower() != 'all':
            collection = self.db['sg_specimen_group']
            result = collection.find_one({"_id": ObjectId(data.group_id)})
            group_names = result["category_names"]
            group_dict_in = group_dict
            group_dict = OrderedDict()
            for each in group_names:
                if each in group_dict_in.keys():
                    group_dict[each] = group_dict_in[each]

        # print "group_num is {}".format(len(group_dict.keys()))
        # print "group_us is {}".format(data.use_group.lower())
        if len(group_dict.keys()) < 2 and data.use_group.lower() != "no":
            info = {'success': False, 'info': "使用分组样本计算(求中值或中位数),需选择两个或两个以上的组别"}
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
            exp_id=str(exp_info['main_id']),
            group_id=data.group_id,
            group_dict=group_dict,
            # quant_method=quant_method,
            use_group=data.use_group,
            geneset_id=data.geneset_id,
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
        name = "Geneset_Cluster" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='geneset cluster main table',
            params=params,
            status="start",
            submit_type=submit_type
        )
        main_id = self.prok_rna.insert_main_table('sg_geneset_cluster', main_info)

        print("++++ exp_id gene_set_id")
        print(express_id + ";" + data.geneset_id)
        # prepare option for workflow
        options = dict(
            exp_matrix=express_id + ";" + data.geneset_id,
            # group_dict=data.group_dict,
            # group=data.group_dict,
            group_dict=json.dumps(group_dict),
            group=json.dumps(group_dict),
            cluster_main_id=str(main_id),
            n_clusters=n_clusters,
            group_id=data.group_id,
            use_group=data.use_group,
            scm=scm,
            scd=scd,
            sct=data.sct,
            gcm=gcm,
            gct=data.gct,
            gcd=gcd,
            update_info=json.dumps({str(main_id): "sg_geneset_cluster"})  # to update sg_status
        )

        # prepare to file
        to_files = ["prok_rna.export_geneset_exp_matrix(exp_matrix)",
                    "prok_rna.export_group(group)",]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'prok_rna.report.geneset_cluster'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(GenesetClusterAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # 更新基因集的使用信息
        self.prok_rna.insert_geneset_info(data.geneset_id, 'sg_geneset_cluster', str(main_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.prok_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.prok_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/prok_rna/geneset_cluster "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna_srna",
            task_type="2",
            submit_location="genesetcluster",
            group_id="5b600dc677b3f3b113a7bd57",
            group_dict=json.dumps({"WT":["WT_1","WT_2","WT_3"],"rcsBKO":[ "rcsBKO_1","rcsBKO_2","rcsBKO_3"]}).replace('"', '\\"'),
            exp_id="5b76a6f3a4e1af5574afbb09",
            geneset_id="5b70dd58a4e1af13f89b4362",
            quant_method='RSEM',
            use_group='no',
            n_clusters='10',
            sct="hierarchy",
            scm="complete",
            scd="correlation",
            gct="hierarchy",
            gcm="average",
            gcd="euclidean",
            submit_type='0',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
