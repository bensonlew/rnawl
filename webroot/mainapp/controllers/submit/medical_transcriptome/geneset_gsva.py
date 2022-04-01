# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mbio.api.to_file.medical_transcriptome import *
from mainapp.libs.signature import check_sig
import os
import unittest
from bson.objectid import ObjectId
import re
from biocluster.config import Config

class GenesetGsvaAction(MedicalTranscriptomeController):
    def __init__(self):
        super(GenesetGsvaAction, self).__init__(instant=False)

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

        with open(geneset_path, 'r') as b:
            c = re.findall("[^\w\\t\\n\-\r]", b.read())
            if c:
                d = list(set(c))
                return json.dumps({'success': False, 'info': "表格中存在以下特殊字符：{}".format("".join(d))})

        with open(geneset_path, 'r') as f_in, open(geneset_path + 'check.gmt', 'w') as f_o:
            for line in f_in:
                if line.startswith("#"):
                    pass
                else:
                    cols = line.strip().split("\t")
                    cols[2] = re.sub(r"[,; ]", "\t", cols[2])
                    f_o.write("\t".join(cols))
            return geneset_path + 'check.gmt'
    @check_sig
    def POST(self):
        project_type = 'medical_transcriptome'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        data = web.input()
        print data
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['level', 'exp_id', 'geneset_id', 'group_id', 'group_dict', 'control_id']
        basic_args += ['geneset_source']
        basic_args += ['es_method', 'min_num', 'max_num']
        # basic_args += ['stat_type', 'stat_cutoff', 'fc']
        options = ['c1', 'c2', 'genesets',
                   'anno_type', 'go_type', 'go_genesets',
                   'kegg_type', 'kegg_type2', 'kegg_genesets',
                   'upload_genesets_file', 'upload_genesets_file_id',
                   ]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:   %s" % arg, 'code': 'C2900301', 'variables': variables}
                return json.dumps(info)
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)

        connect = db['sg_exp']
        record = connect.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id)})
        is_rmbe = str(record['is_rmbe']).lower()
        if is_rmbe == 'false':
            exp_id = data.exp_id
        if is_rmbe == 'true':
            exp_id = str(record['batch_main_id'])
        connect = self.db['sg_specimen_group']
        record_ = connect.find_one({'task_id': data.task_id})
        samples = record_['specimen_names']
        sample_list = sum(samples, [])
        sample_list_str = ';'.join(sample_list)
        w_options = dict(
            sample_list_str=sample_list_str,
            task_id=data.task_id,
            min_num=data.min_num,
            max_num=data.max_num,
            level=data.level,
            matrix=exp_id + ';' + data.level + ';' + is_rmbe,
            group=json.dumps(group_dict),
            cmp=data.control_id,
            id2name = exp_id,
            geneset_id=data.geneset_id,
            es_method=data.es_method,
            # stat_type=data.stat_type,
            # stat_cutoff=data.stat_cutoff,
            # fc=data.fc
        )
        to_file = []
        to_file.append('medical_transcriptome_new.geneset.export_exp_matrix_medical(matrix)')
        to_file.append("medical_transcriptome.export_compare(cmp)")
        to_file.append('medical_transcriptome.export_group_gsva(group)')
        to_file.append('medical_transcriptome_new.geneset.export_seq2name(id2name)')

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
        # 主表数据生成
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            geneset_id=data.geneset_id,
            group_id=data.group_id,
            group_dict=group_dict,
            geneset_source=data.geneset_source,
            es_method=data.es_method,
            min_num=data.min_num,
            max_num=data.max_num,
            level=data.level,
            control_id=data.control_id
            # stat_type=data.stat_type,
            # stat_cutoff=data.stat_cutoff
        )
        # msigdb
        if hasattr(data, 'c1'):
            params.update(dict(c1=data.c1))
        if hasattr(data, 'c2'):
            params.update(dict(c2=data.c2))
        if hasattr(data, 'genesets'):
            params.update(dict(genesets=data.genesets))
        # go&kegg
        if hasattr(data, 'anno_type'):
            params.update(dict(anno_type=data.anno_type))
        if hasattr(data, 'sub_type'):
            params.update(dict(sub_type=data.sub_type))
        if hasattr(data, 'anno_value'):
            # params.update(dict(kegg_type=data.anno_value)
            params.update(dict(anno_value=data.anno_value))
        # custom
        if hasattr(data, 'upload_genesets_file'):
            params.update(dict(upload_genesets_file=data.upload_genesets_file))
        if hasattr(data, 'upload_genesets_file_id'):
            params.update(dict(upload_genesets_file_id=data.upload_genesets_file_id))
        params_json = json.dumps(params, sort_keys=True, separators=(',', ':'))
        task_info = self.medical_transcriptome.get_task_info(data.task_id)
        main_table_name = 'Gsva_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = dict(
            project_sn=task_info['project_sn'],
            task_id=data.task_id,
            status='start',
            name=main_table_name,
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=params_json,
            version='v1',
            group_dict=group_dict,
            level=data.level
        )

        collection_name = 'sg_geneset_gsva'
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

        # prepare to file


        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.geneset_gsva'
        self.set_sheet_data(name=task_name,
                            options=w_options,
                            main_table_name=main_table_name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_file,
                            project_sn=task_info['project_sn'],
                            new_task_id=new_task_id,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(GenesetGsvaAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }
        }
        # if 'group_id' in data and str(data.group_id).lower() != 'all':
        #     _ = self.ref_rna_v2.update_group_is_use(data.task_id, data.group_id)
        # if 'control_id' in data:
        #     _ = self.ref_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
        # task_info['group_dict'] = group_dict
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
        cmd += "s/medical_transcriptome/geneset_gsva "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="gsva",
            exp_id="5f50cacf17b2bf5a6c8bfd88",
            level='G',
            geneset_id='all',
            group_id="5f46228c17b2bf20e4e269e1",
            group_dict=json.dumps({"H1":["H1581_1","H1581_2","H1581_3"], "S1":["SNU16_1","SNU16_2","SNU16_3"], 'H2':['H1581_4','H1581_5', 'H1581_6'], 'H3':['H1581_7', 'H1581_8', 'H1581_9'], 'S2':['SNU16_4', 'SNU16_5', 'SNU16_6'], 'S3':['SNU16_7', 'SNU16_8', 'SNU16_9']}).replace(
                    '"', '\\"'),
            geneset_source='msigdb',
            c1='C5: GO gene sets',
            c2='BP: GO biological process|',
            genesets="GO_METANEPHRIC_MESENCHYME_DEVELOPMENT,GO_LYMPHOCYTE_CHEMOTAXIS,GO_REGULATION_OF_PHOSPHOLIPASE_A2_ACTIVITY,GO_REGULATION_OF_HEMATOPOIETIC_PROGENITOR_CELL_DIFFERENTIATION,GO_CYCLIC_NUCLEOTIDE_BIOSYNTHETIC_PROCESS,GO_REGULATION_OF_TOR_SIGNALING,GO_RESPONSE_TO_PH",
            es_method='max',
            min_num='15',
            max_num='500',
            # stat_type='padjust',
            # stat_cutoff='0.05',
            # fc='2',
            control_id="5f46228c17b2bf20e4e269e0",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
