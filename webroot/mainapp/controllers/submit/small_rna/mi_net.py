# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.libs.signature import check_sig
import unittest
import os
import re
from mainapp.controllers.project.small_rna_controller import SmallRnaController
from biocluster import file
from biocluster.file import getsize, exists


class MiNetAction(SmallRnaController):
    def __init__(self):
        super(MiNetAction, self).__init__(instant=False)

    def check_exp_file(self, file_path, samples):
        '''
        target_dir = 'data' if self.input_data.client == "client01" else 'tsanger-data'
        base_path = "/mnt/ilustre/{}/".format(target_dir)
        '''

        # file_path = base_path + file_path
        header = ""
        if os.path.exists(file_path):
            with open(file_path, 'r') as exp_f:
                header = exp_f.readline().strip()
        elif re.match(r'^\w+://\S+/.+$',  file_path) or re.match(r'/mnt/ilustre', file_path):
            header = file.get_top_lines(file_path, lines=1)[0]
            pass
        else:
            raise Exception("文件传递格式错误 {}".format(file_path))

        samples_exp = header.split("\t")[1:]

        inter_samples = set(samples_exp).intersection(set(samples))

        if len(inter_samples) < 3:
            pass
            # raise Exception("一致样本数量不够相关性分析{} {}".format(samples, samples_exp ))
        else:
            pass
        self.file_path = file_path
        return list(inter_samples)


    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "target_id", "geneset_id", "mirnaset_id", "submit_location", 'task_type']
        corr_args = ['exp_path', 'corr_method', 'corr_threshold', 'corr_pvalue_type', 'corr_pvalue_threshold']
        # check arg
        tf_args = ["tf_db", "tf_pvalue_type", "tf_pvalue_threshold"]
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s" % arg}
                return json.dumps(info)

        exp_info = self.small_rna.get_exp_tpm_params_info(data.task_id)
        samples = self.small_rna.get_samples_by_task_id(data.task_id)

        try:
            target_info = self.small_rna.get_target_params_info(data.target_id)
        except:
            info = {'success': False, 'info': "Target predict result need choose"}
            return json.dumps(info)

        target_dir = target_info["result_dir"]
        if not target_dir.endswith("/"):
            target_dir += "/"

        use_samples = []
        if hasattr(data, "exp_path"):
            print "samples is"
            print samples
            use_samples = self.check_exp_file(data.exp_path, samples)
            if len(use_samples) < 3:
                info = {"success": False, "info": "基因表达文件中匹配的样本数目不够[样本名称要与项目一致，且一致数目不少于3个]"}
                return json.dumps(info)
            else:
                pass

        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        task_info_dict = self.small_rna.get_task_info(data.task_id)
        if "genome_id" in task_info_dict:
            genome_info = self.small_rna.get_genome_info_by_id(task_info_dict['genome_id'])
        else:
            genome_info = self.small_rna.get_genome_info(task_info_dict['organism_name'], task_info_dict['assembly'], task_info_dict['annot_version'])

        '''
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        exp_level = exp_info['exp_level']
        exp_type, quant_method = exp_info['exp_type'], exp_info['method']
        '''
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            geneset_id=data.geneset_id,
            mirnaset_id=data.mirnaset_id,
            exp_id = str(data.exp_id),
            target_id=data.target_id
        )
        if hasattr(data, "exp_path"):
            params.update({
                "exp_path": data.exp_path,
                "corr_method": data.corr_method,
                "corr_threshold": data.corr_threshold,
                "corr_pvalue_type": data.corr_pvalue_type,
                "corr_pvalue_threshold": data.corr_pvalue_threshold
            })

        if hasattr(data, "tf_db"):
            params.update({
                "tf_db": data.tf_db,
                "tf_pvalue_type": data.tf_pvalue_type,
                "tf_pvalue_threshold": data.tf_pvalue_threshold,
            })


        if hasattr(data, "file_id"):
            params.update({
                "file_id": data.file_id
            })

        if hasattr(data, "tf_sub_db"):
            params.update({
                "tf_sub_db": data.tf_sub_db,
            })
        if hasattr(data, "tss_up"):
            params.update({
                "tss_up": data.tss_up,
            })
        if hasattr(data, "tss_down"):
            params.update({
                "tss_down": data.tss_down,
            })

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Minet" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            version = "v1.1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='small_RNA关联分析',
            params=params,
            status="start"
        )
        main_id = self.small_rna.insert_main_table('sg_minet', main_info)

        # prepare option for workflow
        options = {
            "geneset_id": data.geneset_id,
            "mirnaset_id": data.mirnaset_id,
            "main_id": str(main_id),
            "gene_link": data.task_id,
            "gene_detail": data.task_id,
            "genome": genome_info['dna_fa'],
            "species": genome_info['organism_name'],
            "known_pre": task_info_dict["known_pre"],
            "novol_pre": task_info_dict["novol_pre"],
            "update_info": json.dumps({str(main_id): "sg_minet"})  # to update sg_status
        }

        if exists(target_dir + 'known_target.xls') and exists(target_dir + 'novol_target.xls'):
            options.update({
                "known_target": target_dir + 'known_target.xls',
                "novol_target": target_dir + 'novol_target.xls'
            })
        elif exists(target_dir + 'target_predict_detail.xls'):
            options.update({
                "all_target": target_dir + 'target_predict_detail.xls'
            })
        else:
            info = {'success': False, 'info': "Target predict result can not find"}
            return json.dumps(info)


        if hasattr(data, "exp_path"):
            options.update({
                "exp_target": data.exp_path,
                "corr_way": data.corr_method,
                "corr_cutoff": data.corr_threshold,
                "use_samples": "|".join(use_samples),
                # "exp_matrix": exp_info['main_id'],
                "exp_matrix": str(data.exp_id)
            })
            if data.corr_pvalue_type == "pvalue":
                options.update({
                    "pvalue_cutoff": data.corr_pvalue_threshold,
                })
            elif data.corr_pvalue_type == "padjust":
                options.update({
                    "qvalue_cutoff": data.corr_pvalue_threshold,
                })

        if hasattr(data, "tf_db"):
            options.update({
                "db": data.tf_db,
                "mood_pvalue": data.tf_pvalue_threshold,
            })

        if hasattr(data, "tf_sub_db"):
            options.update({
                "tf_sub_db": data.tf_sub_db,
            })
        if hasattr(data, "tss_up"):
            options.update({
                "tss_up": data.tss_up,
            })
        if hasattr(data, "tss_down"):
            options.update({
                "tss_down": data.tss_down,
            })



        to_files = ["small_rna.get_gene_link(gene_link)",
                    "small_rna.get_gene_detail(gene_detail)"
        ]
        if data.geneset_id != "":
            to_files.append("small_rna.export_geneset_list(geneset_id)")
        if data.mirnaset_id != "":
            to_files.append("small_rna.export_geneset_list(mirnaset_id)")
        if hasattr(data, "exp_path"):
            to_files.append("small_rna.export_exp_matrix2(exp_matrix)")

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'small_rna.report.mi_net'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(MiNetAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # task_info['group_dict'] = group_dict
        # task_info['group_dict'] = group_dict
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    '''
    basic_args = ["task_id", "geneset_id", "mirnaset_id", "submit_location", 'task_type']
    corr_args = ['exp_path', 'corr_method', 'corr_threshold', 'corr_pvalue_type', 'corr_pvalue_threshold']
    # check arg
    tf_args = ["tf_db", "tf_pvalue_type", "tf_pvalue_threshold"]
    '''

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/small_rna/mi_net "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="small_rna",
            task_type="2",
            submit_location="mi_net",
            geneset_id="5c05d745a4e1af590965a295",
            mirnaset_id="5c05f2e2a4e1af0c67a5df4b",
            target_id="5c0633e1a4e1af71351834c8",
            exp_path="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/data5/gene.tpm.matrix",
            corr_method='pearson',
            corr_threshold="0.5",
            corr_pvalue_type='pvalue',
            corr_pvalue_threshold="0.5",
            tf_db="vertebrates",
            tf_pvalue_type="pvalue",
            tf_pvalue_threshold="0.5"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)
