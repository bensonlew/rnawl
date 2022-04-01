# -*- coding: utf-8 -*-
# __author__ = 'konghualei,zoujiaxun'

import datetime
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.whole_transcriptome import *


class GenesetClusterAction(WholeTranscriptomeController):
    def __init__(self):
        super(GenesetClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['group_id', 'group_dict', 'level', 'use_group',
                       'geneset_id', 'sct', 'gct',
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
                print info
                return json.dumps(info)
        # if data.level != data.geneset_type:
        #     info = "level is not consistent with geneset_type!"
        #     return json.dumps(info)
        if data.group_id.lower() == 'all' and data.use_group.lower() == 'yes':
            info = {'success': False, 'info': "group 设置为all时,不能选择以分组均值分析", 'code': 'C2901102', 'variables': ''}
            return json.dumps(info)
        exp = self.whole_transcriptome.get_exp_id(data.level, data.task_id)
        exp_id = exp['main_id'].__str__()
        exp_info = self.whole_transcriptome.get_exp_params_info(exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)

        # create main table record
        level = exp_info['level']
        # exp_type, quant_method = exp_info['exp_type'], exp_info['method']

        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            # exp_id=exp_id,
            group_id=data.group_id,
            level=level,
            group_dict=group_dict,
            # quant_method=quant_method,
            use_group=data.use_group,
            geneset_id=data.geneset_id,
            # geneset_type=data.geneset_type,
            sct=data.sct,
            gct=data.gct,
        )
        if hasattr(data, "category"):
            params.update({
                "category": data.category
            })  # 当选择转录本时，添加rna_type参数
        # if hasattr(data, "way"):
        #     params.update({
        #         "way": data.way
        #     })     #当选择mRNA或者lncRNA时，添加表达量指标参数（tpm、fpkm）

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
        name = "Cluster" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            version="v1",
            group_dict=json.loads(data.group_dict, object_pairs_hook=OrderedDict),
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='geneset cluster main table',
            params=params,
            status="start",
            level=level
        )
        main_id = self.whole_transcriptome.insert_main_table('geneset_cluster', main_info)

        # prepare option for workflow
        options = dict(
            exp_matrix=exp_id + ";" + data.geneset_id + ";" + data.level,
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
            # gene_detail="|".join([task_id, data.geneset_type]),
            update_info=json.dumps({str(main_id): "geneset_cluster"})  # to update sg_status
        )
        if hasattr(data, "category"):
            options.update({
                "exp_matrix": exp_id + ";" + data.geneset_id + ";" + data.level + ";" + data.category
            })
        # if hasattr(data, "way"):
        #     options.update({
        #         "exp_matrix": data.exp_id + ";" + data.geneset_id + ";" + data.level + ";" +data.rna_type + ";" + data.way
        #     })
        # prepare to file
        to_files = ["whole_transcriptome.geneset_analysis.export_geneset_exp_matrix(exp_matrix)",
                    "whole_transcriptome.geneset_analysis.export_group(group)",
                    # "whole_transcriptome.geneset_analysis.get_geneset_name(geneset_name)"
                    ]
        # "whole_transcriptome.whole_transcriptome.get_gene_detail_new(gene_detail)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'whole_transcriptome.report.geneset_cluster'
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
        self.whole_transcriptome.insert_geneset_info(data.geneset_id, 'geneset_cluster', str(main_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.whole_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.whole_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/whole_transcriptome/geneset_cluster "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_36088",
            task_type="2",
            submit_location="genesetcluster",
            group_id="5dd1f88e17b2bf10f6f25ced",
            group_dict=json.dumps({"control": ["C2", "C3"], "Acute": ["A1", "A2"], "Stable": ["S1", "S3"]}).replace(
                '"', '\\"'),
            # exp_id="5db7c6af17b2bf4a2fefebb2",
            # exp_id="5dccfdf417b2bf67f7fb5426",
            level="T",
            # level="T",
            geneset_id="5dbb9ef317b2bf6d13a3a814",
            # geneset_id="5dbb9ef317b2bf6d13a3a822",
            # geneset_type="G",
            # quant_method='RSEM',
            use_group='yes',
            n_clusters='10',
            sct="hierarchy",
            scm="complete",
            scd="correlation",
            gct="hierarchy",
            gcm="average",
            gcd="euclidean",
            category='mRNA'

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
