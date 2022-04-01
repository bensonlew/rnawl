# -*- coding: utf-8 -*-

import datetime
import json
import os
import re
import unittest
from collections import OrderedDict
from mbio.api.to_file.medical_transcriptome import *
import web
from bson.objectid import ObjectId

from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig


class TfbsPredictAction(MedicalTranscriptomeController):
    """ TfbsPredict controller"""

    def __init__(self):
        super(TfbsPredictAction, self).__init__(instant=False)
        # specify where is the single workflow for the current task
        self.task_name = 'medical_transcriptome.report.tfbs_predict'
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
        self.expected_args += ['is_corr']
        self.expected_args += ['group_dict']
        self.expected_args += ['group_id']
        self.expected_args += ['tf_predict_id']
        self.expected_args += ['corr_cutoff']  # 相关系数阈值
        self.expected_args += ['corr_pvalue']  # pvalue cutoff value
        # self.expected_args += ['corr_use_padjust'] # 0 表示用pvalue， 1表示用qvalue过滤
        self.expected_args += ['backward']  # promoter 提取 tss上游多少个碱基
        self.expected_args += ['forward']  # promoter 提取，tss 下游多少个碱基
        self.expected_args += ['organism']  # 根据物种信息，到本地数据库获取gtf， genome
        # receive the params from web
        self.input_data = web.input()
        # in which the main table will be created
        self.collection = "sg_tfbs_predict"

    def check_params(self):
        for arg in self.expected_args:
            if not hasattr(self.input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2902801', 'variables': variables}
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
        else:
            self.motifs_user = self.input_data.motifs_user
        if hasattr(self.input_data, "seqdb") and self.input_data.seqdb.lower() != 'none' and self.input_data.seqdb.lower() != '':
            if re.match(r'^\w+://\S+/.+$', self.input_data.seqdb):
                pass
            elif self.input_data.motifs_user.startswith("rere"):
                target_dir = 'data' if self.input_data.client == "client01" else 'tsanger-data'
                base_path = "/mnt/ilustre/{}/".format(target_dir)
                self.seqdb = base_path + self.input_data.seqdb
        else:
            self.seqdb = 'none'

        if self.input_data.motifs_user.lower() == 'none' and self.input_data.tf_selected.lower() == 'none':
            info = {'success': False, 'info': "请选择已知转录因子或按照指定要求上传Motif文件", 'code': 'C2902802', 'variables': ''}
            return json.dumps(info)

        table_info = self.medical_transcriptome.get_main_info(self.input_data.tf_predict_id, 'sg_tf_predict',
                                                            self.input_data.task_id)
        tmp_dict = json.loads(table_info["params"])
        if self.input_data.motifs_user.lower() == 'none' and tmp_dict['s'] == "animal":
            info = {'success': False, 'info': "暂不支持动物靶基因预测", "code": "C2902804"}
            return json.dumps(info)

        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        print "input is {}".format(input_data_dict)
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
        return json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

    def create_main_table(self, packed_params, collection_name):
        result_info = self.medical_transcriptome.get_task_info(self.input_data.task_id)
        self.project_sn = result_info["project_sn"]
        try:
            genome_id = result_info["genome_id"]
        except:
            genome_id = "None"
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
            version='v1'
        )
        main_id = self.medical_transcriptome.insert_main_table(collection_name, main_info)
        if self.input_data.geneset_id.lower() not in ["all", "refall", "none"]:
            self.medical_transcriptome.insert_geneset_info(self.input_data.geneset_id, collection_name, str(main_id))
        return main_id, name, genome_id

    def prepare_workflow_options(self, main_id, collection_name, main_table_name, genome_id):
        input_data_dict = dict(self.input_data)
        options = dict()
        for each in self.expected_args:
            if each not in ['motifs_user', 'seqdb']:
                options[each] = input_data_dict[each]
        # get gene detail
        # result = self.whole_transcriptome.get_table_info_by_task_id(self.input_data.task_id, 'annotation_query')
        # if not result:
        #     raise Exception("query annotation_query failed")
        # class_code_id = str(result['main_id'])
        # update options

        # task_info = self.whole_transcriptome.get_task_info(self.input_data.task_id)
        result_info = self.medical_transcriptome.get_task_info(self.input_data.task_id)
        # long_dir = result_info['sub_output']['long']
        print "result_info is {}".format(result_info)
        new_task_id = self.medical_transcriptome.get_new_id(self.input_data.task_id)
        main_table_data = {'run_id': new_task_id}
        options.update({
            "main_id": str(main_id),
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): collection_name}),  # to update status
            "exp_level": "gene",
            "geneset_list": self.input_data.geneset_id,
            "seq2name": self.input_data.task_id + '|' + 'G',
            "genome_id": genome_id,
            # "gtf":  result_info['output'] + '/other/annotation/all.gtf'
        })
        # get motif species
        table_info = self.medical_transcriptome.get_main_info(self.input_data.tf_predict_id, 'sg_tf_predict',
                                                            self.input_data.task_id)
        self.tf_select_shows = table_info["tf_select_shows"]
        tmp_dict = json.loads(table_info["params"])
        options['motif_species'] = tmp_dict['organism']
        options['motif_species_type'] = tmp_dict['s']
        if self.input_data.motifs_user.lower() != 'none':
            options['motifs_user'] = self.motifs_user
        if self.input_data.seqdb.lower() != 'none' and self.input_data.seqdb.lower() != '':
            options['seqdb'] = self.input_data.seqdb
        else:
            pass
        # prepare to file
        to_files = [
            "medical_transcriptome.export_geneid2tfid_file(tf_predict_id)",
            "medical_transcriptome.get_gene_detail_new(seq2name)",
        ]
        if 'all' not in self.input_data.geneset_id.lower():
            to_files.append("medical_transcriptome.export_gene_list(geneset_list)")
        else:
            pass
        if "all" in input_data_dict['tf_selected'].lower():
            tf_selected_all = [tf.split("(")[0] + "|" + tf.split("(")[1].split(",")[0] for tf in self.tf_select_shows]
            options['tf_selected'] = ",".join(tf_selected_all)
        else:
            options['tf_selected'] = input_data_dict['tf_selected']
        if self.input_data.exp_id.lower() != "none":
            connect_exp = self.db['sg_exp']
            record_exp = connect_exp.find_one({'task_id': self.input_data.task_id, 'main_id': ObjectId(self.input_data.exp_id)})
            is_rmbe = str(record_exp['is_rmbe']).lower()
            if is_rmbe == 'false':
                exp_id = self.input_data.exp_id
            if is_rmbe == 'true':
                exp_id == str(record_exp['batch_main_id'])
            options['exp_matrix'] = exp_id + ";all;G;" + is_rmbe
            to_files.append("medical_transcriptome.export_geneset_exp_matrix_new(exp_matrix)", )
        else:
            options['exp_matrix'] = "none"
        # 把参数交给workflow运行相应的tools， 其中to_file用于准备tool的输入文件
        # options['gtf'] = "none"

        print "option is {}".format(options)
        self.set_sheet_data(
            name=self.task_name,
            options=options,
            main_table_name=main_table_name,  # 设置交互分析结果目录名
            module_type="workflow",
            to_file=to_files,
            project_sn=self.project_sn,
            new_task_id=new_task_id,
            task_id=self.input_data.task_id
        )

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()
        main_id, main_table_name, genome_id = self.create_main_table(packed_params, self.collection)
        self.prepare_workflow_options(main_id, self.collection, main_table_name, genome_id)
        task_info = super(TfbsPredictAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        if 'group_id' in self.input_data and str(self.input_data.group_id).lower() != 'all' and str(
                self.input_data.group_id).lower() != 'none':
            _ = self.medical_transcriptome.update_group_is_use(self.input_data.task_id, self.input_data.group_id)
        if 'control_id' in self.input_data:
            _ = self.medical_transcriptome.update_group_compare_is_use(self.input_data.task_id,
                                                                     self.input_data.control_id)
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
        cmd += "s/medical_transcriptome/tfbs_predict "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",  # maybe you need to change it
            submit_location="tfbspredict",
            # motifs_user="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/whole_rna/test.meme",
            motifs_user="None",
            tf_selected="ENSG00000173153|MA0592.1,ENSG00000008197|MA0003.1",
            geneset_id="all",
            motifs_file_id="",
            organism="Homo_sapiens",
            is_corr="no",
            seqdb="None",
            seqdb_file_id="",
            exp_id="5f50cacf17b2bf5a6c8bfd88",
            group_id="all",
            group_dict=json.dumps({"H1": ["H1581_1", "H1581_2", 'H1581_3'], "H2": ["H1581_4", "H1581_5", 'H1581_6'], "H3": ["H1581_7", "H1581_8", 'H1581_9'], 'S1':["SNU16_1", "SNU16_2", "SNU16_3"], 'S2':["SNU16_4", "SNU16_5", "SNU16_6"], 'S3':["SNU16_7", "SNU16_8", "SNU16_9"]}).replace('"', '\\"'),
            tf_predict_id="5f4e145217b2bf2095bdd6a4",
            thresh="0.0001",
            corr_cutoff="0.5",
            corr_pvalue="0.05",
            backward="450",
            forward="50",

        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
