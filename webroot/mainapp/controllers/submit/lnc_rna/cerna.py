# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.libs.signature import check_sig
import unittest
import os
import re
from biocluster.config import Config
from bson.objectid import ObjectId
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from biocluster import file


class CernaAction(LncRnaController):
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
        return list(inter_samples)


    @check_sig
    def POST(self):
        data = web.input()
        # print "data is {}".format(data)
        basic_args = ["task_id", "exp_level", "exp_id", "geneset_id", "lncset_id", "submit_location", 'task_type', 'corr_method', 'corr_cutoff', 'corr_pvalue_type', 'corr_pvalue_cutoff', 'group_dict', 'group_id']
        corr_args = ['smallrna_exp',  'micorr_method', 'micorr_cutoff', 'micorr_pvalue_type', 'micorr_pvalue_cutoff']
        target_args = ['smallrna', 'mirbase_category', 'mirbase_specie', 'miranda_score', 'miranda_energy', 'miranda_strict', 'ps_robot_score']
        time_now = datetime.datetime.now()
        name = time_now.strftime("%Y%m%d_%H%M%S")

        if hasattr(data, 'mirbase_specie'):
            small_rna_dir = Config().WORK_DIR + "/tmp/" + data.task_id + name + "_mirna.fa"
            with open(Config().WORK_DIR + "/tmp/" + data.task_id + name + "_mirna.fa", 'w') as f_o:
                for spe in data.mirbase_specie.split(","):
                    small_sub_dir  = Config().SOFTWARE_DIR + "/database/mirbase/species/{}/{}.mature.fa".format(spe, spe)
                    with open(small_sub_dir, 'r') as f_in:
                        f_o.write(f_in.read())
        elif hasattr(data, 'smallrna'):
            small_rna_dir = data.smallrna
        # print "small is {}".format(small_rna_dir)

        # check arg

        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s" % arg}
                return json.dumps(info)


        exp_info = self.lnc_rna.get_exp_params_info(data.exp_id, data.task_id)
        try:
            geneset_info = self.lnc_rna.get_main_info_by_record("sg_geneset", main_id=ObjectId(data.geneset_id), task_id=data.task_id)
            lncset_info = self.lnc_rna.get_main_info_by_record("sg_geneset", main_id=ObjectId(data.lncset_id), task_id=data.task_id)
        except:
            info = {'success': False, 'info': "请选择基因集再分析"}
            return json.dumps(info)
        # samples = self.lnc_rna.get_samples_by_task_id(data.task_id)
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        samples = list()
        for each in group_dict:
            samples += group_dict[each]

        if geneset_info["type"][-1] != exp_info["exp_level"]:
            info = {"success": False, "info": "基因集属性与表达量不一致"}
            return json.dumps(info)
        if lncset_info["type"][-1] != exp_info["exp_level"]:
            info = {"success": False, "info": "lncRNA集属性与表达量不一致"}
            return json.dumps(info)


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

        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        task_info_dict = self.lnc_rna.get_task_info(data.task_id)
        genome_info_dict= self.lnc_rna.get_genome_info_by_genome_id(task_info_dict["genome_id"])
        annot_info = self.lnc_rna.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), status="end",  type="origin")



        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            group_dict=group_dict,
            group_id=data.group_id,
            geneset_id=data.geneset_id,
            lncset_id=data.lncset_id,
            exp_level=data.exp_level,
            exp_id=data.exp_id,
            corr_method=data.corr_method,
            corr_cutoff=data.corr_cutoff,
            corr_pvalue_type=data.corr_pvalue_type,
            corr_pvalue_cutoff=data.corr_pvalue_cutoff,
        )
        if data.corr_pvalue_type == "padjust":
            params.update({
                "corr_pvalue_method": data.corr_pvalue_method,
            })
        if hasattr(data, "smallrna_exp"):
            params.update({
                "smallrna_exp": data.smallrna_exp,
                "micorr_cutoff" : data.micorr_cutoff,
                "micorr_pvalue_cutoff" : data.micorr_pvalue_cutoff,
                "micorr_pvalue_type" : data.micorr_pvalue_type
            })
            if data.micorr_pvalue_type == "padjust":
                params.update({
                    "micorr_pvalue_method": data.micorr_pvalue_method
                })

        if hasattr(data, "mi_file_id"):
            params.update({
                "mi_file_id": data.mi_file_id
            })
        if hasattr(data, "miexp_file_id"):
            params.update({
                "miexp_file_id": data.miexp_file_id
            })

        for arg in target_args:
            if hasattr(data, arg):
                params.update({
                    arg: data[arg]
                })

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Cerna" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='ceRNA关联分析',
            params=params,
            status="start"
        )
        main_id = self.lnc_rna.insert_main_table('sg_cerna', main_info)

        if genome_info_dict["taxon"] in ["Animal"]:
            taxon = "Animal"
        else:
            taxon = "Plant"

        # prepare option for workflow

        options = {
            "task_id": data.task_id,
            "type": data.exp_level,
            "geneset_id": data.geneset_id,
            "lncset_id": data.lncset_id,
            "main_id": str(main_id),
            "smallrna": small_rna_dir,
            "mrna": task_info_dict['mrna_fa'],
            "mrna_gtf": task_info_dict['mrna_gtf'],
            "known_lnc_fa": task_info_dict['known_lnc_fa'],
            "novol_lnc_fa": task_info_dict['novel_lnc_fa'],
            "known_lnc_gtf": task_info_dict['known_lnc_gtf'],
            "novol_lnc_gtf": task_info_dict['novel_lnc_gtf'],
            "genome": Config().SOFTWARE_DIR + "/database/Genome_DB_finish/" + genome_info_dict['dna_fa'],
            "group_dict": json.dumps(group_dict),
            "annotation": data.geneset_id + ";" + data.lncset_id,
            "cerna_exp": data.exp_id + ";" + data.geneset_id + ";" + data.lncset_id,
            "corr_method": data.corr_method,
            "corr_cutoff": data.corr_cutoff,
            "group_id": data.group_id,
            "corr_pvalue_type": data.corr_pvalue_type,
            "corr_pvalue_cutoff": data.corr_pvalue_cutoff,
            "taxon": taxon,
            "update_info": json.dumps({str(main_id): "sg_cerna"})  # to update sg_status
        }

        if hasattr(data, "smallrna_exp"):
            options.update({
                "smallrna_exp": data.smallrna_exp,
                "micorr_cutoff" : data.micorr_cutoff,
                "micorr_method" : data.corr_method,
                "micorr_pvalue_cutoff" : data.micorr_pvalue_cutoff,
                "micorr_pvalue_type" : data.micorr_pvalue_type,
                "use_samples": "|".join(use_samples)

            })
        else:
            options.update({
                "use_samples": "|".join(samples)
            })

        if data.corr_pvalue_type == "padjust":
            options.update({
                "corr_padjust_method": data.corr_pvalue_method,
            })

        if hasattr(data, "micorr_pvalue_type") and data.micorr_pvalue_type == "padjust":
            options.update({
                "micorr_padjust_method": data.micorr_pvalue_method
            })
        for arg in target_args:
            if arg in ['smallrna', 'mirbase_category', 'mirbase_specie']:
                pass
            else:
                if hasattr(data, arg):
                    options.update({
                        arg: data[arg]
                    })
        options.update({
            "smallrna": small_rna_dir
        })


        to_files = ["lnc_rna.export_geneset_list(geneset_id)",
                    "lnc_rna.export_geneset_list(lncset_id)",
                    "lnc_rna.export_geneset_lncset_exp_matrix(cerna_exp)",
                    "lnc_rna.export_genename_des(annotation)"
        ]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'lnc_rna.report.ce_rna'
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
        cmd += "s/lnc_rna/cerna "
        cmd += "-b http://192.168.12.101:9090 "

        args = dict(
            task_id="lnc_rna",
            task_type="2",
            submit_location="cerna",
            exp_id="5ca4386c17b2bf1791ec0e25",
            geneset_id="5c99dd3417b2bf5886f20be9",
            lncset_id="5c99da9917b2bf284acd32a5",
            exp_level="G",
            group_id='5c90920e17b2bf5bda46c152',
            group_dict=r'{"Con":["Con1","Con2","Con3","Con4","Con5"],"Vit": ["Vit1","Vit2","Vit3","Vit4","Vit5"]}'.replace(
                '"', '\\"'),
            smallrna_exp="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_lnc_rna/cerna/mirna.exp.xls",
            smallrna="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_lnc_rna/cerna/test_mature.fa",
            corr_method='pearson',
            corr_cutoff="0.5",
            corr_pvalue_type='pvalue',
            corr_pvalue_cutoff="0.5",
            corr_pvalue_method="bh",
            micorr_method='pearson',
            micorr_cutoff="0.5",
            micorr_pvalue_type='pvalue',
            micorr_pvalue_cutoff="0.5",
            micorr_pvalue_method="bh"

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)
