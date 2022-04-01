# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2018-11-19"
# __author2__ = "xuxi 2021-04-14"

import xml.etree.ElementTree as ET
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.protein_transcript_dia_v3_controller import ProteinTranscriptDiaV3Controller
from mainapp.libs.signature import check_sig
import unittest
import os
from bson.objectid import ObjectId
from biocluster.config_trans import ConfigTrans
import json
from mbio.api.to_file.dia import *


class ProteinTranscriptAction(ProteinTranscriptDiaV3Controller):
    def __init__(self):
        super(ProteinTranscriptAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "\ncontroller input data: "+str(data)
        basic_args = ["task_id", "protein_params", 'rna_params']
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
            if arg.lower() == "null":
                info = {'success': False, 'info': "{} : is null or NULL".format(arg)}
                return json.dumps(info)
        exp_info = self.dia.get_exp_params_info_new2(task_id=data.task_id, name="ProteinTable_Origin")
        project_sn = exp_info["project_sn"]
        # data.protein_params = data.protein_params.replace('"',"'")
        # data.rna_params = data.rna_params.replace('"',"'")
        try:
            self.protein_fa = ET.fromstring(data.protein_params).find('stage').find('parameters').find(
                'protein_fasta').text
        except:
            try:
                self.protein_fa = json.loads(data.protein_params)['protein_fasta']
            except:
                info = {'success': False, 'info': "没能找到蛋白序列文件"}
                return json.dumps(info)
            # info = {'success': False, 'info': "没能找到蛋白序列文件"}
            # return json.dumps(info)
        # with open('/mnt/ilustre/users/sanger-dev/fengxml','w') as xml_w:
        #     xml_w.write(data.rna_params)
        #     xml_w.write('\n\n')
        #     xml_w.write(data.protein_params)
        try:
            self.rna_type = ET.fromstring(data.rna_params).find('stage').find('name').text
        except:
            try:
                self.rna_type = json.loads(data.rna_params)['name']
            except:
                info = {'success': False, 'info': "没能解析出转录组流程的路径"}
                return json.dumps(info)
            # info = {'success': False, 'info': "没能解析出转录组流程的路径"}
            # return json.dumps(info)

        try:
            self.rna_task_id = ET.fromstring(data.rna_params).find('task_id').text
        except:
            try:
                self.rna_task_id = json.loads(data.rna_params)['task_id']
            except:
                info = {'success': False, 'info': "没能解析出转录组流程的task_id"}
                return json.dumps(info)
            # info = {'success': False, 'info': "没能解析出转录组流程的task_id"}
            # return json.dumps(info)

        # 根据前端传入的xml解析出rna项目调用的workflow脚本，进而确定出rna来源于那种类型的项目
        if self.rna_type == 'ref_rna.refrna':
            self.rna_type = 'ref_rna'
        if self.rna_type == 'ref_rna_v2.refrna':
            self.rna_type = 'ref_rna_v2'
        if self.rna_type == 'prok_rna.prokrna':
            self.rna_type = 'prok_rna'
        if self.rna_type == 'denovo_rna_v2.denovorna':
            self.rna_type = 'denovo_rna_v2'
        if self.rna_type == 'whole_transcriptome.whole_transcriptome':
            self.rna_type = 'whole_transcriptome'
        if self.rna_type == 'medical_transcriptome.medical_transcriptome':
            self.rna_type = 'medical_transcriptome'

        # 根据rna项目的类型，获得其包含基因与蛋白信息的文件路径
        if self.rna_type == 'denovo_rna_v2':
            task_info, db_version = get_task_info_new('denovo_rna_v2', self.rna_task_id)
            try:
                self.genome_info = self.gene = task_info['assemble_t2g']
            except:
                info = {'success': False, 'info': "没能找到此无参转录组组装信息文件"}
                return json.dumps(info)
            try:
                self.pep = task_info['bedpath'].replace('.bed', '.pep')
            except:
                info = {'success': False, 'info': "没能找到此无参转录组的预测蛋白文件"}
                return json.dumps(info)
            try:
                # self.gene = task_info['unigene_fa'].replace('.fasta', '.unigene.fasta')
                self.gene = task_info['unigene_fa']
            except:
                info = {'success': False, 'info': "没能找到此无参转录组的组装后得到的unigene的fasta文件"}
                return json.dumps(info)
        if self.rna_type == 'ref_rna':
            try:
                self.species = ET.fromstring(data.rna_params).find('stage').find('parameters').find('ref_genome').text
            except:
                try:
                    self.species = json.loads(data.rna_params)['ref_genome']
                except:
                    info = {'success': False, 'info': "没能找到此有参转录组的物种信息"}
                    return json.dumps(info)
                # info = {'success': False, 'info': "没能找到此有参转录组的物种信息"}
                # return json.dumps(info)
            with open(os.path.join(ConfigTrans().SOFTWARE_DIR, 'database/Genome_DB_finish/annot_species.json')) as f:
                js = json.load(f)
                if self.species in js and 'bio_mart_annot' in js[self.species]:
                    self.genome_info = self.gene = os.path.join(ConfigTrans().SOFTWARE_DIR,
                                                                'database/Genome_DB_finish/',
                                                                js[self.species]['bio_mart_annot'])
                    if not os.path.exists(self.genome_info):
                        info = {'success': False, 'info': "没能找到此有参转录组的物种biomart文件"}
                        return json.dumps(info)
                else:
                    info = {'success': False, 'info': "没能找到此有参转录组的物种信息"}
                    return json.dumps(info)
                if self.species in js and 'pep' in js[self.species]:
                    self.pep = os.path.join(ConfigTrans().SOFTWARE_DIR, 'database/Genome_DB_finish/',
                                            js[self.species]['pep'])
                    if not os.path.exists(self.pep):
                        info = {'success': False, 'info': "没能找到此有参转录组的物种pep文件"}
                        return json.dumps(info)
                else:
                    info = {'success': False, 'info': "没能找到此有参转录组的物种信息"}
                    return json.dumps(info)

        if self.rna_type == 'ref_rna_v2':
            task_info, db_version = get_task_info_new('ref_rna_v2', self.rna_task_id)
            try:
                # biomart, type = self.ref_rna_v2.get_des_type(self.rna_task_id)
                biomart, type = get_des_type(task_info, self.rna_task_id, db_version)
                self.genome_info = self.gene = os.path.join(ConfigTrans().SOFTWARE_DIR, 'database/Genome_DB_finish/',
                                                            biomart)
                if not os.path.exists(self.genome_info):
                    info = {'success': False, 'info': "没能找到此有参转录组的物种biomart文件"}
                    return json.dumps(info)
            except:
                info = {'success': False, 'info': "没能找到此有参转录组的物种信息"}
                return json.dumps(info)
            try:
                pep = get_pep(task_info, self.rna_task_id, db_version)
                self.pep = os.path.join(ConfigTrans().SOFTWARE_DIR, 'database/Genome_DB_finish/',
                                        pep)
                if not os.path.exists(self.genome_info):
                    info = {'success': False, 'info': "没能找到此有参转录组的物种pep文件"}
                    return json.dumps(info)
            except:
                info = {'success': False, 'info': "没能找到此有参转录组的物种信息"}
                return json.dumps(info)

        if self.rna_type == 'whole_transcriptome':
            task_info, db_version = get_task_info_new_whole_transcript('whole_transcriptome', self.rna_task_id)
            try:
                # biomart, type = self.ref_rna_v2.get_des_type(self.rna_task_id)
                biomart, type = get_des_type(task_info, self.rna_task_id, db_version)
                self.genome_info = self.gene = os.path.join(ConfigTrans().SOFTWARE_DIR, 'database/Genome_DB_finish/',
                                                            biomart)
                if not os.path.exists(self.genome_info):
                    info = {'success': False, 'info': "没能找到此全转录组的物种biomart文件"}
                    return json.dumps(info)
            except:
                info = {'success': False, 'info': "没能找到此全转录组的物种信息"}
                return json.dumps(info)
            try:
                pep = get_pep(task_info, self.rna_task_id, db_version)
                self.pep = os.path.join(ConfigTrans().SOFTWARE_DIR, 'database/Genome_DB_finish/', pep)
                if not os.path.exists(self.genome_info):
                    info = {'success': False, 'info': "没能找到此全转录组的物种pep文件"}
                    return json.dumps(info)
            except:
                info = {'success': False, 'info': "没能找到此全转录组的物种信息"}
                return json.dumps(info)

        if self.rna_type == 'medical_transcriptome':
            task_info, db_version = get_task_info_new('medical_transcriptome', self.rna_task_id)
            try:
                # biomart, type = self.ref_rna_v2.get_des_type(self.rna_task_id)
                biomart, type = get_des_type(task_info, self.rna_task_id, db_version)
                self.genome_info = self.gene = os.path.join(ConfigTrans().SOFTWARE_DIR, 'database/Genome_DB_finish/',
                                                        biomart)
                if not os.path.exists(self.genome_info):
                    info = {'success': False, 'info': "没能找到此医学转录组的物种biomart文件"}
                    return json.dumps(info)
            except:
                info = {'success': False, 'info': "没能找到此医学转录组的物种信息"}
                return json.dumps(info)
            try:
                pep = get_pep(task_info, self.rna_task_id, db_version)
                self.pep = os.path.join(ConfigTrans().SOFTWARE_DIR, 'database/Genome_DB_finish/',
                                                        pep)
                if not os.path.exists(self.genome_info):
                    info = {'success': False, 'info': "没能找到此医学转录组的物种pep文件"}
                    return json.dumps(info)
            except:
                info = {'success': False, 'info': "没能找到此医学转录组的物种信息"}
                return json.dumps(info)

        if self.rna_type == 'prok_rna':
            # task_info = self.prok_rna.get_task_info(self.rna_task_id)
            print(self.rna_task_id)
            task_info, db_version = get_task_info_new('prok_rna', self.rna_task_id)

            if 'rock_index' in task_info:
                self.genome_info = self.gene = os.path.join(task_info['rock_index'], 'ptt.bed')
                self.pep = os.path.join(task_info['rock_index'], 'cds.faa')
            else:
                info = {'success': False, 'info': "没能找到此原核转录组的基因组信息"}
                return json.dumps(info)

        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            rna_type=self.rna_type,
        )
        # 卡方检验的页面没有单双尾检验
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Protein_transcript" + '_' + self.rna_type + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            rna_task_id=self.rna_task_id,
            name=name,
            rna_type=self.rna_type,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='protein_transcript_dia_v3_analyse',
            status="start",
            protein_fa=self.protein_fa,
            genome_info=self.genome_info,
            rna_pep=self.pep,
            params=params
        )
        # 　self.labelfree来自mainapp/controllers/project
        # /labelfree_and_tmt_controller.py的from
        # mainapp.models.mongo.labelfree_and_tmt import ItraqTmt然后init里面写了self.labelfree = ItraqTmt(bind_object=bind_object)
        # 这里就是接口直接调用函数，而to_file里面的函数是接口造字符串，传递给workflow才使用函数调用
        main_id = self.dia.insert_main_table('sg_p2g_relationship', main_info)

        # prepare option for workflow
        options = {
            "genome_info": self.genome_info,
            "main_id": str(main_id),
            "pep": self.pep,
            "transcript_list": self.gene,
            "protein_faa": self.protein_fa,
            "rna_type": self.rna_type,
            "update_info": json.dumps({str(main_id): "sg_p2g_relationship"}),  # to update sg_status
        }

        to_file_list = []
        if self.rna_type == 'whole_transcriptome':
            options["rna_matrix"] = str(self.rna_task_id) + '__format__' + 'task_id'
            to_file_list.append("protein_transcript_dia_v3.export_exp_matrix_whole_transcriptome(rna_matrix)")
        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'protein_transcript_dia_v3.protein_transcript'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_file_list,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ProteinTranscriptAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            }
        }
        # task_info['group_dict'] = group_dict
        return json.dumps(task_info)


def get_task_info_new(project_type, task_id):
    try:
        ConfigTrans().DBVersion = 0
        db = ConfigTrans().get_mongo_client(mtype=project_type, db_version=0, task_id=task_id)[
            ConfigTrans().get_mongo_dbname(project_type, db_version=0, task_id=task_id)]
        task_info = db['sg_task'].find_one({"task_id": task_id})
        print("task_info: {}".format(task_info))
        task_info['_id']
        return task_info, 0
    except:
        ConfigTrans().DBVersion = 1
        db = ConfigTrans().get_mongo_client(mtype=project_type, db_version=1, task_id=task_id)[
            ConfigTrans().get_mongo_dbname(project_type, db_version=1, task_id=task_id)]
        task_info = db['sg_task'].find_one({"task_id": task_id})
        print("task_info: {}".format(task_info))
        task_info['_id']
        return task_info, 1

def get_task_info_new_whole_transcript(project_type, task_id):
    try:
        ConfigTrans().DBVersion = 0
        db = ConfigTrans().get_mongo_client(mtype=project_type, db_version=0, task_id=task_id)[
            ConfigTrans().get_mongo_dbname(project_type, db_version=0, task_id=task_id)]
        task_info = db['task'].find_one({"task_id": task_id})
        task_info['_id']
        task_info_long_task_id = task_info['options']['long_task_id']
        task_info = db['task'].find_one({"task_id": task_info_long_task_id})['options']

        return task_info, 0
    except:
        ConfigTrans().DBVersion = 1
        db = ConfigTrans().get_mongo_client(mtype=project_type, db_version=1, task_id=task_id)[
            ConfigTrans().get_mongo_dbname(project_type, db_version=1, task_id=task_id)]
        task_info = db['task'].find_one({"task_id": task_id})
        task_info['_id']
        task_info_long_task_id = task_info['options']['long_task_id']
        task_info = db['task'].find_one({"task_id": task_info_long_task_id})['options']
        return task_info, 1

def get_des_type(result_task, task_id, db_version, table_name="sg_annotation_stat"):
    if "genome_id" in result_task.keys():
        genome_id = result_task["genome_id"]
        db = ConfigTrans().get_mongo_client(mtype='ref_rna_v2', db_version=db_version, dydb_forbid=True)[
            ConfigTrans().get_mongo_dbname('ref_rna_v2', db_version=db_version, dydb_forbid=True)]
        col = db["sg_genome_db"]
        genome_info = col.find_one({"genome_id": genome_id})
        db_path = ConfigTrans().SOFTWARE_DIR + "/database/Genome_DB_finish"
        des = os.path.join(db_path, genome_info["bio_mart_annot"])
        des_type = genome_info["biomart_gene_annotype"]
        return des, des_type
    else:
        ConfigTrans().DBVersion = db_version
        db = ConfigTrans().get_mongo_client(mtype='ref_rna_v2', db_version=db_version, task_id=task_id)[
            ConfigTrans().get_mongo_dbname('ref_rna_v2', db_version=db_version, task_id=task_id)]
        collection = db[table_name]
        result = collection.find_one({"task_id": task_id, "type": "origin"})
        species_name = result['species_name']
        json_path = ConfigTrans().SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        json_dict = get_json(json_path)
        des = os.path.join(os.path.split(json_path)[0], json_dict[species_name]["bio_mart_annot"])
        des_type = json_dict[species_name]["biomart_gene_annotype"]
        return des, des_type


def get_pep(result_task, task_id, db_version, table_name="sg_annotation_stat"):
    if "genome_id" in result_task.keys():
        genome_id = result_task["genome_id"]
        db = ConfigTrans().get_mongo_client(mtype='ref_rna_v2', db_version=db_version, dydb_forbid=True)[
            ConfigTrans().get_mongo_dbname('ref_rna_v2', db_version=db_version, dydb_forbid=True)]
        col = db["sg_genome_db"]
        genome_info = col.find_one({"genome_id": genome_id})
        db_path = ConfigTrans().SOFTWARE_DIR + "/database/Genome_DB_finish"
        pep = os.path.join(db_path, genome_info["pep"])
        return pep
    else:
        ConfigTrans().DBVersion = db_version
        db = ConfigTrans().get_mongo_client(mtype='ref_rna_v2', db_version=db_version, task_id=task_id)[
            ConfigTrans().get_mongo_dbname('ref_rna_v2', db_version=db_version, task_id=task_id)]
        collection = db[table_name]
        result = collection.find_one({"task_id": task_id, "type": "origin"})
        species_name = result['species_name']
        json_path = ConfigTrans().SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        json_dict = get_json(json_path)
        pep = os.path.join(os.path.split(json_path)[0], json_dict[species_name]["pep"])
        return pep


def get_json(json_path):
    f = open(json_path, "r")
    json_dict = json.loads(f.read())
    return json_dict


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += '-dbversion {} '.format(1)
        cmd += "s/protein_transcript/protein_transcript "
        cmd += "-b http://wpm2.sanger.com "
        args = dict(
            # task_id="tsg_33097",
            task_id="jssn_svfguokp3g35cbuue5v31s",
            task_type="2",
            rna_type='medical_transcriptome',
            submit_location="protein_transcript",
            protein_params=json.dumps({"protein_fasta":"s3:\\/\\/commonbucket\\/data\\/rerewrweset\\/files\\/va7r02f0841hb52h3epqirta6h\\/upload\\/3ae9c02309d3fe52660a7c7d47e05385.fasta"}).replace('"', '\\"'),
            rna_params=json.dumps({"name":"medical_transcriptome.medical_transcriptome","task_id":"fald_ng6r6in6abtut2tl7ast3e","ref_genome":""}).replace('"', '\\"'),
            contract_id='',
            signature='20eeb1f7c1aa69c2f58c7f36ea2521707ca65cef'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

"""<workflow><member_id>m_188</member_id><member_type>1</member_type><cmd_id>140</cmd_id><project_sn>188_5bea79e58826b</project_sn><task_id>tsg_32890</task_id><bucket>rerewrweset</bucket><version>2</version><mapping_file type="sanger" format="mapping_file" alias="mapping_file.txt">s3://commonbucket/files/m_188/188_5bea79e58826b/mapping_file_1542258555.txt</mapping_file><stage><id>cmd_140_1542258555</id><name>ref_rna.refrna</name><type>workflow</type><parameters><fq_type>PE</fq_type><strand_specific>False</strand_specific><is_duplicate>True</is_duplicate><fastq_dir type=\'sanger\' format=\'fastq_dir\' alias=\'rawfastq\' dir=\'true\'></fastq_dir><taxonmy>Animal</taxonmy><ref_genome>Homo_sapiens</ref_genome><group_table type=\'sanger\' format=\'group_table\' alias=\'group.txt\'>s3://commonbucket/files/m_188/188_5bea79e58826b/trans/10b0a63fff691f7c05a74f2e27aa07f8.txt</group_table><control_file type=\'sanger\' format=\'control_file\' alias=\'control.txt\'>s3://commonbucket/files/m_188/188_5bea79e58826b/trans/6c090453c1946032e55633f17e4336df.txt</control_file><exp_way>fpkm</exp_way><nr_evalue>1e-3</nr_evalue><string_evalue>1e-3</string_evalue><kegg_evalue>1e-3</kegg_evalue><swissprot_evalue>1e-3</swissprot_evalue><kegg_database>All</kegg_database><assemble_method>stringtie</assemble_method><seq_method>Hisat</seq_method><diff_method>DESeq2</diff_method><diff_fdr_ci>0.05</diff_fdr_ci><fc>2</fc><snp_analyze>True</snp_analyze><as_analyze>True</as_analyze></parameters></stage></workflow>
"""  # 有参v1传过来的参数样式

"""<workflow><member_id>m_188</member_id><member_type>1</member_type><cmd_id>172</cmd_id><project_sn>188_5bc7e83458246</project_sn><task_id>tsg_32492</task_id><bucket>rerewrweset</bucket><version>2</version><stage><id>cmd_172_1539829436</id><name>labelfree_and_tmt.labelfreetmt</name><type>workflow</type><parameters><labelfree_or_tmt>labelfree</labelfree_or_tmt><protein type='sanger' format='labelfree_and_tmt.common' alias='protein.txt'>s3://commonbucket/files/m_188/188_5bc7e83458246/labelfree_file_test/c7bb9121480dd04870f92283288a5428.txt</protein><psm type='sanger' format='labelfree_and_tmt.common' alias='psm.txt'>s3://commonbucket/files/m_188/188_5bc7e83458246/labelfree_file_test/db2cdfe209923650d8745b714c966749.txt</psm><peptide type='sanger' format='labelfree_and_tmt.common' alias='peptide.txt'>s3://commonbucket/files/m_188/188_5bc7e83458246/labelfree_file_test/8dc59255b3ddd8d96a2263938bbe6d31.txt</peptide><protein_information type='sanger' format='labelfree_and_tmt.common' alias='Protein_information.xls'>s3://commonbucket/files/m_188/188_5bc7e83458246/labelfree_file_test/23e105d3fc9582f8337392e701cf3662.xls</protein_information><ratio_exp type='sanger' format='labelfree_and_tmt.ratio_exp' alias='exp.txt'>s3://commonbucket/files/m_188/188_5bc7e83458246/labelfree_file_test/99065cf9f137f0dff598eb1059fac9d8.txt</ratio_exp><protein_fasta type='sanger' format='labelfree_and_tmt.common' alias='DB.fasta'>s3://commonbucket/files/m_188/188_5bc7e83458246/labelfree_file_test/6ccd75b599e6e51ee3f44afd69859e48.fasta</protein_fasta><protein_group type='sanger' format='labelfree_and_tmt.group_table' alias='group.txt'>s3://commonbucket/files/m_188/188_5bc7e83458246/labelfree_file_test/10b0a63fff691f7c05a74f2e27aa07f8.txt</protein_group><protein_control type='sanger' format='labelfree_and_tmt.compare_table' alias='control.txt'>s3://commonbucket/files/m_188/188_5bc7e83458246/labelfree_file_test/6c090453c1946032e55633f17e4336df.txt</protein_control><go_evalue>1e-5</go_evalue><go_identity>0.98</go_identity><kegg_class>All</kegg_class><kegg_org></kegg_org><kegg_evalue>1e-5</kegg_evalue><kegg_identity>0.98</kegg_identity><sub_loc>Animal</sub_loc><pfam_evalue>1e-5</pfam_evalue><method_type>t.test</method_type><correct_method></correct_method><pvalue>0.05</pvalue><fc_up>1.2</fc_up><fc_down>0.83</fc_down><sam_analysis>True</sam_analysis><ppi_category>All</ppi_category><ppi_species></ppi_species><data_source>Uniprot</data_source><ref_set>True</ref_set></parameters></stage></workflow>
"""
