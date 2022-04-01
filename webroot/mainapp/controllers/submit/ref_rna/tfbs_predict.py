# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mbio.api.to_file.ref_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os
import re
from bson.objectid import ObjectId


class TfbsPredictAction(RefRnaController):
    """ TfbsPredict controller"""
    def __init__(self):
        super(TfbsPredictAction, self).__init__(instant=False)
        # specify where is the single workflow for the current task
        self.task_name = 'ref_rna.report.tfbs_predict'
        # specify all params needed, they will be saved as value of 'params' in the main table
        self.expected_args = ["task_id", "submit_location", 'task_type']
        self.expected_args += ['tf_selected']  # 如"Gene_ID(TF_id1,family),Gene_ID(TF_id2,family)"时，传回"TF_id1,TF_id2"
        self.expected_args += ['motifs_user']  # 前端传回的相对路劲, 没有时传'none'
        self.expected_args += ['motifs_file_id']
        self.expected_args += ['geneset_id']  # get theses genes' promoter,可以为'all'
        self.expected_args += ['seqdb']  # 客户上传fasta, 无则为'none'
        self.expected_args += ['seqdb_file_id']
        self.expected_args += ['thresh']  # FIMO的pvalue 阈值
        # self.expected_args += ['qv_thresh']  # 0 表示用pvalue， 1表示用qvalue过滤
        self.expected_args += ['exp_id']  # 可以为'none'
        self.expected_args += ['group_dict']
        self.expected_args += ['group_id']
        self.expected_args += ['tf_predict_id']
        self.expected_args += ['corr_cutoff']  # 相关系数阈值
        self.expected_args += ['corr_pvalue']  # pvalue cutoff value
        # self.expected_args += ['corr_use_padjust'] # 0 表示用pvalue， 1表示用qvalue过滤
        self.expected_args += ['backward'] # promoter 提取 tss上游多少个碱基
        self.expected_args += ['forward'] # promoter 提取，tss 下游多少个碱基
        self.expected_args += ['organism'] # 根据物种信息，到本地数据库获取gtf， genome
        # receive the params from web
        self.input_data = web.input()
        # in which the main table will be created
        self.collection = "sg_tfbs_predict"

    def check_params(self):
        print self.input_data
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        # do some other check and modify arg if necessary
        if self.input_data.motifs_user.lower() != 'none':
            if re.match(r'^\w+://\S+/.+$', self.input_data.motifs_user):
                self.motifs_user = self.input_data.motifs_user
            elif self.input_data.motifs_user.startswith("rere"):
                target_dir = 'data' if self.input_data.client == "client01" else 'tsanger-data'
                base_path = "/mnt/ilustre/{}/".format(target_dir)
                self.motifs_user = base_path + self.input_data.motifs_user
        else:
            self.motifs_user = self.input_data.motifs_user
        if self.input_data.seqdb.lower() != 'none':
            if re.match(r'^\w+://\S+/.+$', self.input_data.seqdb):
                pass
            elif self.input_data.motifs_user.startswith("rere"):
                target_dir = 'data' if self.input_data.client == "client01" else 'tsanger-data'
                base_path = "/mnt/ilustre/{}/".format(target_dir)
                self.seqdb = base_path + self.input_data.seqdb
        else:
            self.seqdb = 'none'

        if self.input_data.motifs_user.lower() == 'none' and self.input_data.tf_selected.lower() == 'none':
            info = {'success': False, 'info': "请选择已知转录因子或按照指定要求上传Motif文件"}
            return json.dumps(info)

        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        for each in self.expected_args:
            if each == "task_type":
                # maybe we do not need to do so
                params_dict[each] = int(input_data_dict[each])
            elif each == "group_dict":
                if input_data_dict[each].lower() == 'none' or (not input_data_dict[each]):
                    params_dict[each] = input_data_dict[each]
                else:
                    # deal with dumped dict param from web
                    params_dict[each] = json.loads(input_data_dict[each], object_pairs_hook=OrderedDict)
            else:
                params_dict[each] = input_data_dict[each]
        return json.dumps(params_dict, sort_keys=True, separators=(',',':'))

    def create_main_table(self, packed_params, collection_name):
        result_info = self.ref_rna.get_task_info(self.input_data.task_id)
        self.project_sn = result_info["project_sn"]
        name = "TfbsPredict_"
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            # type=self.input_data.?,  # You may need to save some other info for convenience
            project_sn=self.project_sn,
            task_id=self.input_data.task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='TfbsPredict main table',
            params=packed_params,
            status="start",
            tf_predict_id=ObjectId(self.input_data.tf_predict_id),
        )
        main_id = self.ref_rna.insert_main_table(collection_name, main_info)
        if self.input_data.geneset_id.lower() not in ["all", "refall", "none"]:
            self.ref_rna.insert_geneset_info(self.input_data.geneset_id, collection_name, str(main_id))
        return main_id, name

    def prepare_workflow_options(self, main_id, collection_name, main_table_name):
        input_data_dict = dict(self.input_data)
        options = dict()
        for each in self.expected_args:
            options[each] = input_data_dict[each]
        # get gene detail
        sg_task_info = self.ref_rna.get_task_info(self.input_data.task_id)
        if 'is_demo' in sg_task_info and int(sg_task_info['is_demo']) != 0:
            task_id = sg_task_info['demo_id']
            result = self.ref_rna.get_table_info_by_task_id(task_id, 'sg_express_class_code')
        else:
            result = self.ref_rna.get_table_info_by_task_id(self.input_data.task_id, 'sg_express_class_code')
        if not result:
            raise Exception("query sg_express_class_code failed")
        class_code_id = str(result['_id'])
        # update options
        options.update({
            "main_id": str(main_id),
            "update_info": json.dumps({str(main_id): collection_name}), # to update sg_status
            "exp_level": "gene",
            "geneset_list": self.input_data.geneset_id,
            "class_code_id": class_code_id,
        })
        # get motif species
        table_info = self.ref_rna.get_main_info(self.input_data.tf_predict_id, 'sg_tf_predict')
        tmp_dict = json.loads(table_info["params"])
        options['motif_species'] = tmp_dict['organism']
        options['motif_species_type'] = tmp_dict['s']
        if self.input_data.motifs_user.lower() != 'none':
            options['motifs_user'] = self.motifs_user
        # if self.input_data.seqdb.lower() != 'none':
        options['seqdb'] = self.seqdb
        # prepare to file
        to_files = [
            "ref_rna.export_geneid2tfid_file(tf_predict_id)",
            "ref_rna.get_gene_detail(class_code_id)",
        ]
        if 'all' not in self.input_data.geneset_id.lower():
            to_files.append("ref_rna.export_gene_list(geneset_list)")
        if self.input_data.exp_id.lower() != "none":
            options['exp_matrix'] = self.input_data.exp_id+";all"
            to_files.append("ref_rna.export_geneset_exp_matrix(exp_matrix)",)
        else:
            options['exp_matrix'] = "none"
        # 把参数交给workflow运行相应的tools， 其中to_file用于准备tool的输入文件
        self.set_sheet_data(
            name=self.task_name,
            options=options,
            main_table_name=main_table_name,  # 设置交互分析结果目录名
            module_type="workflow",
            to_file=to_files,
            project_sn=self.project_sn,
            task_id=self.input_data.task_id
        )

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()
        main_id, main_table_name = self.create_main_table(packed_params, self.collection)
        self.prepare_workflow_options(main_id, self.collection, main_table_name)
        task_info = super(TfbsPredictAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the web_api func. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/ref_rna/tfbs_predict "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="tsg_29268",
            task_type="2",  # maybe you need to change it
            submit_location="tfbspredict",
            motifs_user="none",
            tf_selected="gene15353|Cla023007,xx|e_gw1.9.523.1",
            geneset_id="all",
            organism="Myzus_persicae",
            seqdb="None",
            exp_id="5addb8cca4e1af76149d0e52",
            group_id="all",
            group_dict=json.dumps({"A":["5addb8c9a4e1af76149cf926","5addb8c9a4e1af76149cf927","5addb8c9a4e1af76149cf928"],"B":["5addb8c9a4e1af76149cf929","5addb8c9a4e1af76149cf92a","5addb8c9a4e1af76149cf92b"],"C":["5addb8c9a4e1af76149cf92c","5addb8c9a4e1af76149cf92d","5addb8c9a4e1af76149cf92e"],"D":["5addb8c9a4e1af76149cf92f","5addb8c9a4e1af76149cf930","5addb8c9a4e1af76149cf931"],"E":["5addb8c9a4e1af76149cf932","5addb8c9a4e1af76149cf933","5addb8c9a4e1af76149cf934"],"F":["5addb8c9a4e1af76149cf935","5addb8c9a4e1af76149cf936","5addb8c9a4e1af76149cf937"]}).replace('"', '\\"'),
            tf_predict_id="5afd0313a4e1af09e2104936",
            thresh="0.0001",
            qv_thresh="0",
            corr_cutoff="0.5",
            corr_pvalue="0.05",
            corr_use_padjust="0",
            backward="450",
            forward="50",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
