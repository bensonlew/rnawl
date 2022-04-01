# -*- coding: utf-8 -*-

import web
import json
from mainapp.libs.signature import check_sig
import datetime
from mainapp.controllers.project.small_rna_controller import SmallRnaController
from mbio.api.to_file.small_rna import *
from mainapp.models.mongo.small_rna import *
import os
import unittest
from collections import OrderedDict

class GenesetClusterAction(SmallRnaController):
    def __init__(self):
        super(GenesetClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ['task_id', 'submit_location', 'task_type', 'group_id', 'group_dict',
                      'use_group', 'geneset_id', 'sct', 'gct', 'exp_id']

        # check whether arguments exist
        re_info = {'success': False, 'info': None}

        for arg in basic_args:
            if not hasattr(data, arg):
                re_info['info'] = 'Lack argument: {}'.format(arg)
                return json.dumps(re_info)

        if data.group_id.lower() == 'all' and data.use_group.lower() == 'yes':
            re_info['info'] = 'use_group can not be set to yes while group is set as all'
            return json.dumps(re_info)

        geneset_info = self.small_rna.get_main_info(data.geneset_id, "sg_geneset", data.task_id)
        print "geneset_info is {}".format(geneset_info)
        if geneset_info.get("gene_length", 0) == 0:
            info = {'success': False, 'info': "该基因集为空，不可以进行聚类分析"}
            return json.dumps(info)
        if hasattr(data, "n_clusters") and  geneset_info.get("gene_length", 0) < int(data.n_clusters):
            info = {'success': False, 'info': "基因数量应该大于等于子聚类数目"}
            return json.dumps(info)


        task_id = data.task_id
        exp_info = self.small_rna.get_exp_info(task_id=task_id)
        project_sn = exp_info['project_sn']
        exp_id = str(data.exp_id)
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)

        # exp_type, quant_method = exp_info['exp_type'], exp_info['method']
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            group_id=data.group_id,
            group_dict=group_dict,
            # quant_method=quant_method,
            use_group=data.use_group,
            geneset_id=data.geneset_id,
            # geneset_type=data.geneset_type,
            sct=data.sct,
            gct=data.gct,
            exp_id=exp_id
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
        name = "GeneSet_Cluster" + '_'
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
        )
        main_id = self.small_rna.insert_main_table('sg_geneset_cluster', main_info)

        # prepare option for workflow
        options = dict(
            exp_matrix=exp_id + ";" + data.geneset_id,
            group_dict=data.group_dict,
            group=data.group_dict,
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
            # gene_detail=task_id,
            update_info=json.dumps({str(main_id): "sg_geneset_cluster"})  # to update sg_status
        )

        # prepare to file
        to_files = [
            "small_rna.export_geneset_exp_matrix(exp_matrix)",
            "small_rna.export_group(group)",
            # "small_rna.get_gene_detail(gene_detail)"
        ]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'small_rna.report.geneset_cluster'
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
        self.small_rna.insert_geneset_info(data.geneset_id, "sg_geneset_cluster", str(main_id))

        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.small_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.small_rna.update_group_compare_is_use(data.task_id, data.control_id)

        return json.dumps(task_info)


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """
        def test_this(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/small_rna/geneset_cluster "
            cmd += "-b http://192.168.12.102:9090 "
            args = dict(
                task_id="tsg_33072",  # 按照sg_exp修改
                task_type="2",
                submit_location="genesetcluster",
                group_id="5c1b44fda4e1af70b4bd3e23",  # sg_specimen_group 的 main_id
                group_dict=json.dumps({"GH": ["GH1", "GH2"],
                                       "NFAs": ["NFAs1", "NFAs2"],
                                       "Normal": ["Normal1", "Normal2"],
                                       "PRL": ["PRL1", "PRL2"]}).replace('"', '\\"'),
                # exp_id="5be4f3e4a4e1af79b83858fd",  # 按照sg_exp修改
                # exp_level="A",  # 按照sg_exp修改
                geneset_id="5c1c8233a4e1af7be4955b82",
                # geneset_type="A",  # 按照sg_exp修改
                # quant_method='RSEM',
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

    unittest.main()
