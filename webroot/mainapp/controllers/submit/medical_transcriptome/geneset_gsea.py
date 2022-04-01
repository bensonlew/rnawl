#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/28 13:09
@file    : gsea.py
"""
from bson.objectid import ObjectId
from biocluster.config import Config
import web
import json
import datetime
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
import os
import re
import unittest
from collections import OrderedDict
from mainapp.models.mongo.medical_transcriptome import *
from mbio.api.to_file.medical_transcriptome import *

class GenesetGseaAction(MedicalTranscriptomeController):
    def __init__(self):
        super(GenesetGseaAction, self).__init__(instant=False)

    def __check_args(self, name, data_obj):
        if not hasattr(data_obj, name):
            error = {'success': False, 'info': '%s参数缺失!' % name, 'code': 'C2904005', 'variables': [name]}
            return json.dumps(error)
        return True
    def check_user_geneset(self, s3file):
        # target_dir = 'data' if self.input_data.client == "client01" else 'tsanger-data'
        # base_path = "/mnt/ilustre/{}/".format(target_dir)

        geneset_path = s3file
        # if os.path.exists(geneset_path):
        #     pass
        # elif re.match(r'^\w+://\S+/.+$', geneset_path) or re.match(r'/mnt/ilustre', geneset_path):
        #     inter_dir = self.create_tmp_dir(self.input_data.task_id, "gsea/")
        #     geneset_path = self.download_from_s3(geneset_path, inter_dir=inter_dir)
        # elif geneset_path.startswith("rerewrweset"):
        #     inter_dir = self.create_tmp_dir(self.input_data.task_id, "wgcna_pipeline/")
        #     geneset_path = self.download_from_s3("/mnt/ilustre/data/" + geneset_path, inter_dir=inter_dir)
        # elif re.match(r'rerewrweset', geneset_path):
        #     inter_dir = self.create_tmp_dir(self.input_data.task_id, "wgcna_pipeline/")
        #     geneset_path = self.download_from_s3(geneset_path, inter_dir=inter_dir)
        # else:
        #     raise Exception("文件传递格式错误 {}".format(geneset_path))

        inter_dir = self.create_tmp_dir(self.input_data.task_id, 'gsea')
        geneset_path = self.download_from_s3(geneset_path, inter_dir=inter_dir)
        print geneset_path
        geneset_name = os.path.basename(geneset_path)
        with open(geneset_path, 'r') as b:
            c = re.findall("[^\w\\t\ \n\-\r\#]", b.read())
            if c:
                d = list(set(c))
                # return json.dumps({'success': False, 'info': "表格中存在以下特殊字符：{}".format("".join(d))})
                raise Exception("表格中存在以下特殊字符{}".format("".join(d)))

        with open(geneset_path, 'r') as f_in, open(inter_dir + '/{}_check.gmt'.format(geneset_name), 'w') as f_o:
            for line in f_in:
                if line.startswith("#"):
                    pass
                else:
                    cols = line.strip().split("\t")
                    cols[2] = re.sub(r"[,; \n]", "\t", cols[2])
                    f_o.write("\t".join(cols)+"\n")
        print  inter_dir + '/{}_check.gmt'.format(geneset_name)
        return inter_dir + '/{}_check.gmt'.format(geneset_name)

    @check_sig
    def POST(self):
        """
        preranked : true or false
        preranked_file: path
        dict(name="rnk", type="infile", format="ref_rna_v2.ref_common"),

        exp_id
        group_id
        group_dict[逗号分隔]
        geneset_id
        dict(name="matrix", type="infile", format="ref_rna_v2.ref_common"),
        # 表型，分组信息
        dict(name="g1", type="string"),
        dict(name="g2", type="string"),
        dict(name="group", type="infile", format="ref_rna_v2.ref_common"),


        # msigdb
        # 基因集来源
        geneset_source
        c1, c2[CC/...], c3[all/list-逗号分隔] 或者 upload_geneset_file
        # 基因集
        dict(name="gmx", type="infile", format="ref_rna_v2.ref_common"),

        sort_method
        dict(name="metric", type="string", default='Signal2Noise'),

        permutation_type [phenotype/gene_set]
        dict(name="permute", type="string", default='phenotype'),

        min_num
        dict(name="set_max", type="int", default=500),
        max_num
        dict(name="set_min", type="int", default=1)
        :return:
        """
        project_type = 'medical_transcriptome'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        data = web.input()
        self.input_data = data
        print "data is {}".format(data)
        params_json = {}
        w_options = {}
        to_file = []
        required = ['preranked', 'submit_location', 'task_type', 'task_id', 'exp_id', 'group_id', 'group_dict', 'geneset_source', 'task_id', 'min_num', 'max_num', 'level', 'geneset_source']
        optional = ['preranked_file', 'preranked_file_id',  'sort_method', 'permutation_type', 'geneset_id',
                    'c1', 'c2', 'genesets',
                    'anno_type', 'go_type', 'go_genesets',
                    'kegg_type', 'kegg_type2', 'kegg_genesets',
                    'upload_genesets_file', 'upload_genesets_file_id',
                    'plot_top_x', 'control_group', 'sub_type', 'anno_value'
        ]
        for name in required:
            info = self.__check_args(name, data)
            if info is not True:
                return info
            else:
                pass

            params_json[name] = getattr(data, name)
            if name in ('task_id', 'min_num', 'max_num', 'level'):
                w_options[name] = getattr(data, name)
        for name in optional:
            if hasattr(data, name):
                params_json[name] = getattr(data, name)

        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # 按对照组排序
        if hasattr(data, 'control_group'):
            if getattr(data, 'control_group') not in group_dict:
                error = {'success': False, 'info': '对照组必须在所选分组内'}
                return json.dumps(error)
            sort_dict = {}
            for g in group_dict.keys():
                if g == data.control_group:
                    sort_dict[g] = 1
                else:
                    sort_dict[g] = 0
            group_dict = OrderedDict(sorted(group_dict.items(), key=lambda x:sort_dict[x[0]]))
        print "gsea preranked"
        # 预排序判断

        if hasattr(data, "preranked") and data.preranked == 'true':
            preranked_file = data.preranked_file
            w_options['rnk'] = preranked_file
            w_options['preranked'] = True
        else:
            for name in ('sort_method', 'permutation_type'):
                info = self.__check_args(name, data)
                if info is not True:
                    return info
                w_options[name] = getattr(data, name)

            if len(group_dict.keys()) != 2:
                error = {'success': False, 'info': '所选分组数量必须为两个', "code" : "C2904006"}
                print "return error {}".format(json.dumps(error))
                return json.dumps(error)
            else:
                if min([len(x) for x in group_dict.values()]) < 2:
                    error = {'success': False, 'info': '每个分组样本数量必须大于二个', "code" : "C2904007"}
                    return json.dumps(error)

        connect = db['sg_exp']
        record = connect.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id)})
        is_rmbe = str(record['is_rmbe']).lower()
        if is_rmbe == 'false':
            exp_id = data.exp_id
        if is_rmbe == 'true':
            exp_id = str(record['batch_main_id'])
        w_options['matrix'] = exp_id + ';' + data.level  + ';' + is_rmbe
        if len(group_dict.keys()) == 2:
            w_options['g1'], w_options['g2'] = group_dict.keys()
        else:
            pass

        samples = []
        for sample in group_dict.values():
            samples.extend(sample)
        print "samples is {}".format(samples)
        w_options['group'] = json.dumps(group_dict)
        w_options['geneset_id'] = data.geneset_id
        to_file.append('medical_transcriptome_new.geneset.export_exp_matrix_medical(matrix)')

        print "gsea geneset"

        # 基因集判断
        w_options['geneset_source'] = data.geneset_source

        if data.geneset_source == 'user_defined':
            info = self.__check_args('upload_genesets_file', data)
            if info is not True:
                return info
            w_options['gmx'] = self.check_user_geneset(data.upload_genesets_file)
        elif data.geneset_source == 'msigdb':
            annot_info = self.medical_transcriptome.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id),  type="origin")
            species = annot_info["species_name"].replace("_", " ")
            if species in ["Mus musculus", "Homo sapiens", "Rattus norvegicus"]:
                pass
            else:
                error = {'success': False, 'info': '该项目物种不可以使用msigdb数据库进行分析', "code" : "C2904008"}
                return json.dumps(error)
            for name in ('c1', 'c2', 'genesets'):
                info = self.__check_args(name, data)
                if info is not True:
                    return info
            # for name in ('sub_type', 'anno_value', 'genesets'):
            #     info = self.__check_args(name, data)
            #     if info is not True:
            #         return info
            c1 = data.c1
            c2 = data.c2
            w_options['gmx'] = ";".join([species, c1, c2, data.genesets])
            to_file.append('medical_transcriptome_new.geneset.export_msigdb_genesets(gmx)')
        elif data.geneset_source == 'annotation':
            info = self.__check_args("anno_type", data)
            if info is not True:
                return info
            params_json["anno_type"] = getattr(data, "anno_type")
            if data.anno_type == "GO":
                # info = self.__check_args("go_type", data)
                info = self.__check_args('sub_type', data)
                if info is not True:
                    return info
                go_type = data.sub_type
                w_options['go_list'] = data.anno_type + ";" + go_type + ";"
                w_options['go_type'] = go_type
                # if go_type == "defined" and hasattr(data, "go_genesets"):
                if go_type == "defined" and hasattr(data, "anno_value"):
                    w_options['go_genesets'] = data.anno_value
                to_file.append('medical_transcriptome_new.geneset.export_anno_genesets(go_list)')
            if data.anno_type == "KEGG":
                # info = self.__check_args("kegg_type", data)
                info = self.__check_args('sub_type', data)
                if info is not True:
                    return info
                kegg_type = data.sub_type
                if kegg_type == "defined":
                    kegg_type2 = "defined"
                else:
                    # info = self.__check_args("kegg_type2", data)
                    info = self.__check_args("anno_value", data)
                    if info is not True:
                        return info
                    kegg_type2 = data.anno_type

                w_options['gmx'] = data.anno_type + ";" + kegg_type + ";" + kegg_type2 + ";"
                if kegg_type == 'defined':
                    if data.anno_type == "":
                        w_options['gmx'] += "all"
                    else:
                        w_options['gmx'] += data.anno_type
                # if hasattr(data, "kegg_genesets"):
                #     if data.kegg_genesets == "":
                #         w_options['gmx'] += "all"
                #     else:
                #         w_options['gmx'] += data.kegg_genesets
                else:
                    w_options['gmx'] += "all"
                to_file.append('medical_transcriptome_new.geneset.export_anno_genesets(gmx)')
            # if data.anno_type.lower() == "do":
            #     w_options['gmx'] = data.anno_type + ";"
            #     if hasattr(data, "do_genesets"):
            #         if data.do_genesets == "":
            #             w_options['gmx'] += "all"
            #         else:
            #             w_options['gmx'] += data.do_genesets
            #     to_file.append('medical_transcriptome_new.geneset.export_anno_genesets(gmx)')
            if data.anno_type.lower() == "reactome":
                w_options['gmx'] = data.anno_type + ";" + exp_id + ';'
                if hasattr(data, "anno_value"):
                    if data.anno_value == "":
                        w_options['gmx'] += "all"
                    else:
                        w_options['gmx'] += data.anno_value
                to_file.append('medical_transcriptome_new.geneset.export_anno_genesets(gmx)')
            if data.anno_type.lower() == "disgenet":
                w_options['gmx'] = data.anno_type + ";" + exp_id + ';'
                if hasattr(data, "anno_value"):
                    if data.anno_value == "":
                        w_options['gmx'] += "all"
                    else:
                        w_options['gmx'] += data.anno_value
                to_file.append('medical_transcriptome_new.geneset.export_anno_genesets(gmx)')
            if data.anno_type.lower() == "do":
                w_options['do_list'] = data.anno_type + ";" + data.anno_value
                if hasattr(data, "anno_value"):
                    w_options['do_genesets'] = data.anno_value
                to_file.append('medical_transcriptome_new.geneset.export_anno_genesets(do_list)')
        w_options['genes_detail'] = {'task_id': data.task_id, 'level': data.level}
        to_file.append('medical_transcriptome_new.geneset.export_genes_detail_medical(genes_detail)')

        if hasattr(data, 'plot_top_x'):
            w_options['plot_top_x'] = data.plot_top_x

        print "updata table"
        # 主表数据生成
        params_json.update({"submit_location": data.submit_location,
                            "task_type": int(data.task_type),
                            "group_dict": group_dict})
        task_info = self.medical_transcriptome.get_task_info(data.task_id)
        main_table_name = 'Gsea_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('version', 'v1'),
            ('group_dict', group_dict),
            ('samples', samples)
        ]
        collection_name = 'sg_geneset_gsea'
        main_table_id = self.medical_transcriptome.insert_main_table(collection_name, mongo_data)
        new_task_id = self.medical_transcriptome.get_new_id(task_info['task_id'])
        main_table_data = {'run_id': new_task_id}
        update_info = {str(main_table_id): collection_name}

        options = {
            'main_table_data': main_table_data,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
        }
        w_options.update(options)
        print "gsea set data sheet"

        task_name = 'medical_transcriptome.report.geneset_gsea'
        self.set_sheet_data(name=task_name, options=w_options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'],
                            new_task_id=new_task_id
                            )

        task_info = super(GenesetGseaAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}

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
            cmd += "s/medical_transcriptome/geneset_gsea "
            cmd += "-b http://bcl.tsg.com "
            args = dict(
                submit_location="gsea",
                level="G",
                geneset_id="all",
                preranked="false",
                # geneset_id="5cf4c84791fc21a3418b4567",
                exp_id="5f865eb817b2bf5a36286dc8",
                geneset_type="G",
                group_id="5f46228c17b2bf20e4e269e1",
                group_dict=json.dumps({"H1":["H1581_1","H1581_2","H1581_3"], "S1":["SNU16_1","SNU16_2","SNU16_3"]}).replace(
                    '"', '\\"'),
                geneset_source="annotation",
                anno_type='Reactome',
                reactome_genesets="R-HSA-4717374,R-HSA-3781865,R-HSA-5668914",
                # geneset_kegg="5b0bfa56a4e1af2d43a6c12a",
                task_id="kkh85oi4fthuvdvjkkjmccv848",
                task_type="2",
                sort_method="Signal2noise",
                permutation_type="phenotype",
                min_num="15",
                max_num="500",

            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


    unittest.main()
