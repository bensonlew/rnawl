#!/usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'zhaozhipeng,20190319'

import datetime
import unittest

import web

from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.lnc_rna import *
import json
from collections import OrderedDict

class GenesetClusterAction(LncRnaController):
    def __init__(self):
        super(GenesetClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'exp_id', 'exp_level', 'use_group',
                       'geneset_id', 'geneset_type', 'sct', 'gct',
                       # 'scm', 'gcm',
                       # 'scd', 'gcd',
                       # 'n_clusters',
                       ]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901101', 'variables': variables}
                return json.dumps(info)

        if data.group_id.lower() == 'all' and data.use_group.lower() == 'yes':
            info = {'success': False, 'info': "group 设置为all时,不能选择以分组均值分析", 'code': 'C2901102', 'variables': ''}
            return json.dumps(info)
        # data.exp_id： 表达主表id
        task_id = data.task_id
        exp_id = data.exp_id
        exp_info = self.lnc_rna.get_exp_params_info(exp_id, task_id)
        project_sn = exp_info["project_sn"]
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)

        # create main table record
        exp_level = exp_info['exp_level']
        # exp_type, quant_method = exp_info['exp_type'], exp_info['method']

        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=exp_id,
            group_id=data.group_id,
            exp_level=exp_level,
            group_dict=group_dict,
            # quant_method=quant_method,
            use_group=data.use_group,
            geneset_id=data.geneset_id,
            geneset_type=data.geneset_type,
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

        #   add geneset name
        geneset_info = self.lnc_rna.get_main_info(data.geneset_id, 'sg_geneset', data.task_id)
        if not geneset_info:
            info = {"success": False, "info": "基因集不存在，请确认参数是否正确！!", 'variables': ''}
            return json.dumps(info)
        geneset_name = geneset_info['name']

        name = geneset_name + "_GeneSet_Cluster" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='geneset cluster main table',
            geneset_id=data.geneset_id,
            geneset_type=data.geneset_type,
            params=params,
            status="start",
            type=exp_level,
        )
        main_id = self.lnc_rna.insert_main_table('sg_geneset_cluster', main_info)

        # prepare option for workflow
        options = dict(
            exp_matrix=exp_id + ";" + data.geneset_id + ";" + data.geneset_type,
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
            gene_detail=task_id,
            update_info=json.dumps({str(main_id): "sg_geneset_cluster"})  # to update sg_status
        )

        # prepare to file src/mbio/api/to_file/lnc_geneset.py
        to_files = ["lnc_geneset.export_geneset_exp_matrix(exp_matrix)",
                    "lnc_rna.export_group(group)",
                    "lnc_rna.get_gene_detail(gene_detail)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'lnc_rna.report.geneset_cluster'
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
        self.lnc_rna.insert_geneset_info(data.geneset_id, name, str(main_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.lnc_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.lnc_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
            cmd += "s/lnc_rna/geneset_cluster "
            cmd += "-b http://192.168.12.101:9090 "
            args = dict(
                task_id="lnc_rna",
                # 提交类型, 2, Int, Conf配置
                task_type="2",
                # 无, genesetcluster, string, Conf配置
                submit_location="genesetcluster",

                # -------------------------- 切换数据, 基因集, 表达量结果表 ------------------------------ #
                # 表达量的id是在基因还是转录本水平上, G/T, string
                exp_level="T",
                # 基因集的id是在基因还是转录本水平上, T/LT/G/LG/, string
                geneset_type="lncRNA,mRNA",
                # str, geneset的主表id
                geneset_id="5c90894a17b2bf407167bc7f,5c90890017b2bf407167bc7d",
                # quant_method为废弃参数
                # quant_method='RSEM',

                # --------------------------- 分组方案, 以分组的均值分析 --------------------------------- #
                # str, 分组方案的主表id
                group_id="5caae8e917b2bf6ddc645eaf",
                # 表示分组的字典, string
                group_dict=json.dumps({"Con": ["Con1", "Con2", "Con3", "Con4", "Con5"],
                                       "Vit": ["Vit1", "Vit2", "Vit3", "Vit4", "Vit5"]}).replace('"', '\\"'),

                # 表达量统计表, str, 表达量的主表id
                exp_id="5cc00b1917b2bf57ab50a705",
                # str, 是否以分组均值分析
                use_group='no',

                # -------------------------------- 聚类分析方法参数 -------------------------------------- #
                # 基因聚类方法, 有三种，hierarchy(对应层级聚类), kmeans, no(表示不聚类）, str， scm是gene cluster type的缩写
                gct="hierarchy",
                # 基因聚类方式, complete, str， scm是gene cluster method的缩写
                gcm="average",
                # 基因距离算法, average, str， scm是gene cluster distance的缩写
                gcd="euclidean",
                # 基因子聚类数目, 10, int
                n_clusters='6',

                # 样本聚类方法, 有两种，hierarchy(对应层级聚类) no(表示不聚类, str， scm是sample cluster type的缩写
                sct="hierarchy",
                # 样本聚类方式, complete, str， scm是sample cluster method的缩写
                scm="complete",
                # 样本距离算法, correlation, str， scm是sample cluster distance的缩写
                scd="correlation",
            )
            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


    unittest.main()
