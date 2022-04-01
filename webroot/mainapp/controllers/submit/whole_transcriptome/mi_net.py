# -*- coding: utf-8 -*-

import datetime
import json
import os
import re
import unittest
from collections import OrderedDict

import web
from biocluster import file
from bson.objectid import ObjectId

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class MiNetAction(WholeTranscriptomeController):
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
        elif re.match(r'^\w+://\S+/.+$', file_path) or re.match(r'/mnt/ilustre', file_path):
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
        basic_args = ["task_id", "target_id", "geneset_id", "mirnaset_id", "submit_location", 'task_type', "is_exp"]
        corr_args = ['corr_method', 'corr_threshold', 'corr_pvalue_type', 'corr_pvalue_threshold']
        # check arg
        tf_args = ["tf_db", "tf_pvalue_type", "tf_pvalue_threshold"]
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s" % arg}
                return json.dumps(info)

        exp_info = self.whole_transcriptome.get_main_info_by_record("exp", task_id=data.task_id, level="T")
        samples = self.whole_transcriptome.get_samples_by_task_id(data.task_id)
        samples = list(set(samples))

        try:
            target_info = self.whole_transcriptome.get_main_info_by_record("target_mi",
                                                                           main_id=ObjectId(data.target_id))
        except:
            info = {'success': False, 'info': "Target predict result need choose"}
            return json.dumps(info)

        target_dir = target_info["result_dir"]
        if not target_dir.endswith("/"):
            target_dir += "/"

        use_samples = samples
        '''
        if hasattr(data, "exp_path"):
            print "samples is"
            print samples
            use_samples = self.check_exp_file(data.exp_path, samples)
            if len(use_samples) < 3:
                info = {"success": False, "info": "基因表达文件中匹配的样本数目不够[样本名称要与项目一致，且一致数目不少于3个]"}
                return json.dumps(info)
            else:
                pass
        '''

        print "samples is {}".format(samples)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        task_info_dict = self.whole_transcriptome.get_task_info(data.task_id)
        genome_info = self.whole_transcriptome.get_genome_info_by_genome_id(task_info_dict['genome_id'])

        group_dict = OrderedDict({"all": use_samples})

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
            target_id=data.target_id
        )

        params.update({
            "is_exp": data.is_exp,
            "corr_method": data.corr_method,
            "corr_threshold": data.corr_threshold,
            "corr_pvalue_type": data.corr_pvalue_type,
            "corr_pvalue_threshold": data.corr_pvalue_threshold,

        })

        if hasattr(data, "tf_db"):
            params.update({
                "tf_db": data.tf_db,
                "tf_pvalue_threshold": data.tf_pvalue_threshold,
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
            version="v1.1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='whole_transcriptome关联分析',
            params=params,
            status="start"
        )
        main_id = self.whole_transcriptome.insert_main_table('minet', main_info)

        # prepare option for workflow
        options = {
            "geneset_id": data.geneset_id,
            "mirnaset_id": data.mirnaset_id,
            "task_id": data.task_id,
            "main_id": str(main_id),
            "gene_link": data.task_id,
            "gene_detail": str(exp_info['main_id']) + ";" + data.geneset_id,
            "group_id": "all",
            "group_dict": json.dumps(group_dict),
            "genome": genome_info['dna_fa'],
            "species": genome_info['organism_name'],
            "target": data.target_id + ";" + data.mirnaset_id,
            "known_pre": task_info_dict["known_pre"],
            "novol_pre": task_info_dict["novol_pre"],
            "update_info": json.dumps({str(main_id): "minet"})  # to update sg_status
        }

        options.update({
            "exp_matrix": str(exp_info['main_id']) + ";" + data.mirnaset_id + ";miRNA;T",
            "corr_way": data.corr_method,
            "corr_cutoff": 0,
            "use_samples": "|".join(use_samples),
            "exp_target": str(exp_info['main_id']) + ";" + data.geneset_id + ";mRNA;T",
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

        to_files = [
            "whole_transcriptome.advance.export_geneset_exp_matrix(exp_matrix)",
            "whole_transcriptome.advance.export_geneset_exp_matrix(exp_target)",
            "whole_transcriptome.advance.export_gene_list(geneset_id)",
            "whole_transcriptome.advance.export_gene_list(mirnaset_id)",
            "whole_transcriptome.advance.get_gene_link(gene_link)",
            "whole_transcriptome.advance.export_geneset_genes_name(gene_detail)",
            "whole_transcriptome.advance.export_geneset_mitarget(target)",
        ]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'whole_transcriptome.report.mi_net'
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
        cmd += "s/whole_transcriptome/mi_net "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="whole_transcriptome",
            task_type="2",
            submit_location="mi_net",
            is_exp="yes",
            exp_id="5db7c6ca17b2bf4a2ff053a3",
            target_id="5dac51d617b2bf31566eefaa",
            geneset_id="5db032b317b2bf3cdfb10638",
            mirnaset_id="5dbb9ef317b2bf6d13a3a814",
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
