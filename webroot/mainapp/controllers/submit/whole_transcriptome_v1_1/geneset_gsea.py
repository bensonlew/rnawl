#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/28 13:09
@file    : gsea.py
"""

import datetime
import json
import os
import re
import unittest
from collections import OrderedDict
from bson.objectid import ObjectId
import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class GenesetGseaAction(WholeTranscriptomeController):
    def __init__(self):
        super(GenesetGseaAction, self).__init__(instant=False)

    def __check_args(self, name, data_obj):
        if not hasattr(data_obj, name):
            error = {'success': False, 'info': '%s参数缺失!' % name, 'code': 'C2904005', 'variables': [name]}
            return json.dumps(error)
        return True

    def check_user_geneset(self, s3file):
        target_dir = 'data' if self.input_data.client == "client01" else 'tsanger-data'
        base_path = "/mnt/ilustre/{}/".format(target_dir)

        geneset_path = s3file
        if os.path.exists(geneset_path):
            pass
        elif re.match(r'^\w+://\S+/.+$', geneset_path) or re.match(r'/mnt/ilustre', geneset_path):
            inter_dir = self.create_tmp_dir(self.input_data.task_id, "gsea/")
            geneset_path = self.download_from_s3(geneset_path, inter_dir=inter_dir)
        elif geneset_path.startswith("rerewrweset"):
            inter_dir = self.create_tmp_dir(self.input_data.task_id, "wgcna_pipeline/")
            geneset_path = self.download_from_s3("/mnt/ilustre/data/" + geneset_path, inter_dir=inter_dir)
        elif re.match(r'rerewrweset', geneset_path):
            inter_dir = self.create_tmp_dir(self.input_data.task_id, "wgcna_pipeline/")
            geneset_path = self.download_from_s3(geneset_path, inter_dir=inter_dir)
        else:
            raise Exception("文件传递格式错误 {}".format(geneset_path))

        with open(geneset_path, 'r') as f_in, open(geneset_path + 'check.gmt', 'w') as f_o:
            for line in f_in:
                if line.startswith("#"):
                    pass
                else:
                    cols = line.strip().split("\t")
                    cols[2] = re.sub(r"[,; ]", "\t", cols[2])
                    f_o.write("\t".join(cols) + "\n")
            return geneset_path + 'check.gmt'

    @check_sig
    def POST(self):
        """
        preranked : true or false
        preranked_file: path
        dict(name="rnk", type="infile", format="whole_transcriptome.ref_common"),

        exp_id
        group_id
        group_dict[逗号分隔]
        geneset_id
        dict(name="matrix", type="infile", format="whole_transcriptome.ref_common"),
        # 表型，分组信息
        dict(name="g1", type="string"),
        dict(name="g2", type="string"),
        dict(name="group", type="infile", format="whole_transcriptome.ref_common"),


        # msigdb
        # 基因集来源
        geneset_source
        c1, c2[CC/...], c3[all/list-逗号分隔] 或者 upload_geneset_file
        # 基因集
        dict(name="gmx", type="infile", format="whole_transcriptome.ref_common"),

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
        data = web.input()
        self.input_data = data
        print "data is {}".format(data)
        params_json = {}
        w_options = {}
        to_file = []
        required = ['preranked', 'submit_location', 'task_type', 'task_id', 'group_id', 'group_dict', 'geneset_source',
                    'task_id', 'min_num', 'max_num', 'level', 'geneset_source', 'exp_id']
        optional = ['preranked_file', 'preranked_file_id', 'sort_method', 'permutation_type', 'geneset_id',
                    'c1', 'c2', 'genesets',
                    'anno_type', 'go_type', 'go_genesets',
                    'kegg_type', 'kegg_type2', 'kegg_genesets',
                    'upload_genesets_file', 'upload_genesets_file_id',
                    'plot_top_x'
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
        # self.exp_info = self.whole_transcriptome.get_main_info_by_record("exp", task_id=self.input_data.task_id,
        #                                                                  level=self.input_data.level, way="tpm")
        self.exp_dict = self.whole_transcriptome.db['exp'].find_one({'task_id': data.task_id, "level": data.level, 'main_id': ObjectId(data.exp_id)})
        is_rmbe = str(self.exp_dict['is_rmbe']).lower()
        if 'is_rmbe' not in self.exp_dict or is_rmbe == 'false':
            exp_id = str(self.exp_dict['main_id'])
        if is_rmbe == 'true':
            exp_id = str(self.exp_dict['batch_main_id'])

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
                error = {'success': False, 'info': '所选分组数量必须为两个', "code": "C2904006"}
                print "return error {}".format(json.dumps(error))
                return json.dumps(error)
            else:
                if min([len(x) for x in group_dict.values()]) < 2:
                    error = {'success': False, 'info': '每个分组样本数量必须大于二个', "code": "C2904007"}
                    return json.dumps(error)
        w_options['matrix'] = ";".join([exp_id, data.geneset_id, 'mRNA', data.level, is_rmbe])
        if len(group_dict.keys()) == 2:
            w_options['g1'], w_options['g2'] = group_dict.keys()
        else:
            pass

        samples = list()
        for sample in group_dict.values():
            samples.extend(sample)
        print "samples is {}".format(samples)
        w_options['group'] = data.group_dict
        w_options['group_id'] = data.group_id
        w_options['group_dict'] = data.group_dict
        w_options['geneset_id'] = data.geneset_id
        # to_file.append('whole_transcriptome.advance.export_exp_matrix(matrix)')
        to_file.append('whole_transcriptome_v1_1.advance.export_geneset_exp_matrix(matrix)')

        print "gsea geneset"

        # 基因集判断
        w_options['geneset_source'] = data.geneset_source

        if data.geneset_source == 'user_defined':
            info = self.__check_args('upload_genesets_file', data)
            if info is not True:
                return info
            w_options['gmx'] = self.check_user_geneset(data.upload_genesets_file)
        elif data.geneset_source == 'msigdb':
            task_info = self.whole_transcriptome.get_task_info(data.task_id)
            # annot_info = self.whole_transcriptome.get_main_info_by_record("annotation_stat", task_id=str(data.task_id),  type="origin")
            genome_id = task_info["genome_id"]
            genome_dict = self.whole_transcriptome.get_genome_info_by_genome_id(genome_id)
            species = genome_dict["organism_name"].replace("_", " ")
            if species in ["Mus musculus", "Homo sapiens", "Rattus norvegicus"]:
                pass
            else:
                error = {'success': False, 'info': '该项目物种不可以使用msigdb数据库进行分析', "code": "C2904008"}
                return json.dumps(error)
            for name in ('c1', 'c2', 'genesets'):
                info = self.__check_args(name, data)
                if info is not True:
                    return info
            w_options['gmx'] = ";".join([species, data.c1, data.c2, data.genesets])
            to_file.append('whole_transcriptome.advance.export_msigdb_genesets(gmx)')
        elif data.geneset_source == 'annotation':
            info = self.__check_args("anno_type", data)
            if info is not True:
                return info
            params_json["anno_type"] = getattr(data, "anno_type")
            if data.anno_type == "GO":
                info = self.__check_args("go_type", data)
                if info is not True:
                    return info
                w_options['go_list'] = data.anno_type + ";" + data.go_type + ";"
                w_options['go_type'] = data.go_type
                if data.go_type == "defined" and hasattr(data, "go_genesets"):
                    w_options['go_genesets'] = data.go_genesets
                to_file.append('whole_transcriptome.advance.export_anno_genesets(go_list)')
            if data.anno_type == "KEGG":
                info = self.__check_args("kegg_type", data)
                if info is not True:
                    return info
                if data.kegg_type == "defined":
                    kegg_type2 = "defined"
                else:
                    info = self.__check_args("kegg_type2", data)
                    if info is not True:
                        return info
                    kegg_type2 = data.kegg_type2

                w_options['gmx'] = data.anno_type + ";" + data.kegg_type + ";" + kegg_type2 + ";"
                if hasattr(data, "kegg_genesets"):
                    if data.kegg_genesets == "":
                        w_options['gmx'] += "all"
                    else:
                        w_options['gmx'] += data.kegg_genesets
                else:
                    w_options['gmx'] += "all"
                to_file.append('whole_transcriptome.advance.export_anno_genesets(gmx)')
        w_options['genes_detail'] = {'task_id': data.task_id, 'type': 'origin'}
        to_file.append('whole_transcriptome.advance.export_genes_detail(genes_detail)')

        if hasattr(data, 'plot_top_x'):
            w_options['plot_top_x'] = data.plot_top_x

        print "updata table"
        # 主表数据生成
        params_json.update({"submit_location": data.submit_location,
                            "task_type": int(data.task_type),
                            "group_dict": group_dict})
        task_info = self.whole_transcriptome.get_task_info(data.task_id)
        main_table_name = 'Gsea_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('group_dict', group_dict),
            ('samples', samples),
            ('version', 'v1.1')
        ]
        collection_name = 'geneset_gsea'
        main_table_id = self.whole_transcriptome.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
        }
        w_options.update(options)
        print "gsea set data sheet"

        task_name = 'whole_transcriptome.report.geneset_gsea'
        self.set_sheet_data(name=task_name, options=w_options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

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
            cmd += "s/whole_transcriptome/geneset_gsea "
            cmd += "-b http://bcl.tsg.com "
            args = dict(
                submit_location="gsea",
                db_type="whole_transcriptome",
                level="G",
                geneset_id="all",
                preranked="false",
                geneset_type="G",
                group_id="all",
                group_dict=json.dumps({'NFD': ['NFD1', 'NFD2', 'NFD3', 'NFD4'], 'HFD': ['HFD1', 'HFD2', 'HFD3', 'HFD4']}).replace(
                '"', '\\"'),
                geneset_source="msigdb",
                c1="C5: GO gene sets",
                c2="BP: GO biological process",
                genesets="GO_METANEPHRIC_MESENCHYME_DEVELOPMENT,GO_LYMPHOCYTE_CHEMOTAXIS,GO_REGULATION_OF_PHOSPHOLIPASE_A2_ACTIVITY,GO_REGULATION_OF_HEMATOPOIETIC_PROGENITOR_CELL_DIFFERENTIATION,GO_CYCLIC_NUCLEOTIDE_BIOSYNTHETIC_PROCESS,GO_REGULATION_OF_TOR_SIGNALING,GO_RESPONSE_TO_PH",
                # geneset_kegg="5b0bfa56a4e1af2d43a6c12a",
                task_id="tsg_38079",
                task_type="2",
                sort_method="log2_Ratio_of_Classes",
                permutation_type="phenotype",
                min_num="1",
                max_num="500",
                exp_id='5fa3992e17b2bf0df27f4aa0'

            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


    unittest.main()
