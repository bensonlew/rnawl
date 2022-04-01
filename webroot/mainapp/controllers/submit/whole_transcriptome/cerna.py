# -*- coding: utf-8 -*-

import datetime
import json
import os
import re
import unittest
from collections import OrderedDict

import web
from biocluster import file

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class CernaAction(WholeTranscriptomeController):
    def __init__(self):
        super(CernaAction, self).__init__(instant=False)

    def check_exp_file(self, file_path, samples):
        '''
        检查上传文件样本名是否一致
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
        return list(inter_samples)

    @check_sig
    def POST(self):
        data = web.input()
        # print "data is {}".format(data)
        basic_args = ["task_id", "miset_id", "target_id", "submit_location", 'task_type',
                      'micorr_cutoff', 'micorr_pvalue_type', 'micorr_pvalue_cutoff', 'micorr_pvalue_method',
                      'cecorr_cutoff', 'cecorr_pvalue_type', 'cecorr_pvalue_cutoff', 'cecorr_pvalue_method',
                      'group_dict', 'group_id', 'retain_type',
                      'mi_min_num'
                      ]

        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s" % arg}
                return json.dumps(info)

        exp_info = self.whole_transcriptome.get_main_info_by_record("exp", task_id=data.task_id, level="T", way="tpm")
        '''
        try:
            geneset_info = self.whole_transcriptome.get_main_info_by_record("geneset", main_id=ObjectId(data.miset_id), task_id=data.task_id)
        except:
            info = {'success': False, 'info': "请选择基因集再分析"}
            return json.dumps(info)
        '''
        # samples = self.whole_transcriptome.get_samples_by_task_id(data.task_id)
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        samples = list()
        for each in group_dict:
            samples += group_dict[each]

        '''
        if geneset_info["type"][-1] != exp_info["exp_level"]:
            info = {"success": False, "info": "基因集属性与表达量不一致"}
            return json.dumps(info)
        '''

        use_samples = []
        if hasattr(data, "smallrna_exp"):
            print "*** samples is"
            print samples
            use_samples = self.check_exp_file(data.smallrna_exp, samples)
            print use_samples
            if len(use_samples) < 3:
                info = {"success": False, "info": "small表达文件样本数量不够只有{}与所选样本不一致".format(samples)}
                return json.dumps(info)
            else:
                pass

        task_info_dict = self.whole_transcriptome.get_task_info(data.task_id)
        # print "task_info_dict is {}".format(task_info_dict)
        # genome_info_dict= self.whole_transcriptome.get_genome_info_by_genome_id(task_info_dict["options"]["genome_id"])
        project_sn = task_info_dict["project_sn"]
        task_id = data.task_id

        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
        )

        for par in ["miset_id", "target_id", "targetset_id",
                    'micorr_cutoff', 'micorr_pvalue_type', 'micorr_pvalue_cutoff', 'micorr_pvalue_method',
                    'cecorr_cutoff', 'cecorr_pvalue_type', 'cecorr_pvalue_cutoff', 'cecorr_pvalue_method',
                    'group_dict', 'group_id', 'mi_min_num', 'retain_type']:
            if par == 'group_dict':
                params.update({par: group_dict})
            else:
                params.update({par: getattr(data, par)})

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "ceRNA" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='ceRNA关联分析',
            params=params,
            status="start",
            version='v1'
        )
        main_id = self.whole_transcriptome.insert_main_table('cerna', main_info)

        '''
        if genome_info_dict["taxon"] in ["Animal"]:
            taxon = "Animal"
        else:
            taxon = "Plant"
        '''

        options = {
            "task_id": data.task_id,
            "main_id": str(main_id),
            "level": "T",
            "mi_fa": data.task_id + ";" + "miRNA;all;all",
            "group_dict": json.dumps(group_dict),
            "annotation": {'task_id': data.task_id, 'type': 'origin'},
            "cerna_exp": str(exp_info['main_id']) + ";all" + ";lncRNA,mRNA,circRNA" + ";T",
            "mirna_exp": str(exp_info['main_id']) + ";" + data.miset_id + ";miRNA" + ";T",
            "target": data.target_id + ";" + data.miset_id + ";" + data.targetset_id,
            "update_info": json.dumps({str(main_id): "cerna"})  # to update sg_status
        }

        for par in ["miset_id",
                    'micorr_cutoff', 'micorr_pvalue_type', 'micorr_pvalue_cutoff', 'micorr_pvalue_method',
                    'cecorr_cutoff', 'cecorr_pvalue_type', 'cecorr_pvalue_cutoff', 'cecorr_pvalue_method',
                    'group_dict', 'group_id', 'mi_min_num', 'retain_type']:
            options.update({par: getattr(data, par)})

        to_files = [
            "whole_transcriptome.advance.get_seq(mi_fa)",
            "whole_transcriptome.advance.export_geneset_exp_matrix(cerna_exp)",
            "whole_transcriptome.advance.export_geneset_exp_matrix(mirna_exp)",
            "whole_transcriptome.advance.export_genes_detail(annotation)",
            "whole_transcriptome.advance.export_geneset_mitarget(target)"
        ]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'whole_transcriptome.report.ce_rna'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(CernaAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            }
        }
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
        cmd += "s/whole_transcriptome/cerna "
        cmd += "-b http://bcl.tsg.com "

        args = dict(
            task_id="tsg_36088",
            task_type="2",
            submit_location="cerna",
            # exp_id="5db7c6ca17b2bf4a2ff053a3",
            target_id="5dd518c217b2bf732808d2e3",
            miset_id="all",
            targetset_id="all",
            group_id="all",
            group_dict=json.dumps({"A": ["A1", "A2"], "C": ["C2", "C3"], "S": ["S1", "S3"]}).replace('"', '\\"'),
            cecorr_cutoff="0.5",
            cecorr_pvalue_type='pvalue',
            cecorr_pvalue_cutoff="0.5",
            cecorr_pvalue_method="bh",
            micorr_cutoff="0.5",
            micorr_pvalue_type='pvalue',
            micorr_pvalue_cutoff="0.5",
            micorr_pvalue_method="bh",
            mi_min_num="1"

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)
