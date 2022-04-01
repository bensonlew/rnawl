# -*- coding: utf-8 -*-
# __author__ = 'hongdongxuan'
# last_update liubinxu

import web
import json
import datetime
import unittest
import os
import re
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig


class GenesetPpiAction(LncRnaController):
    """
    蛋白质互作网络的接口
    """

    species_list = [30611, 9598, 61853, 9593, 9606, 9544, 9483, 30608, 9601, 9478, 10141, 10020, 10090, 9986, 10116,
                    43179, 37347, 9685, 9913, 9739, 9669, 9796, 132908, 59463, 9646, 9823, 9785, 9813, 9371, 9361,
                    28377, 9031, 13735, 9103, 59729, 8049, 31033, 8090, 8083, 69293, 99883, 8128, 7955, 13616, 9258,
                    9305, 9315, 7897, 7757, 7719, 51511, 6239, 7227, 4932, 15368, 4513, 4641, 4533, 4538, 4555,
                    4558, 4577, 59689, 3702, 3711, 3847, 3694, 4081, 4113, 29760, 88036, 3218, 3055, 45157]

    def __init__(self):
        super(GenesetPpiAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'taxon', 'species', 'submit_location', 'combine_score', 'geneset_type', 'task_type', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                variables = []
                variables.append(argu)
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        if str(data.species) == '':
            variables = []
            variables.append(data.species)
            info = {"success": False, "info": "%s请选择参考物种" % data.species}
            return json.dumps(info)
        '''
        if int(data.species) not in self.species_list:
            info = {"success": False, "info": "不能进行蛋白质互作分析，因为数据库中不存在该物种的蛋白质互作组数据！"}
            return json.dumps(info)
        '''
        task_name = 'lnc_rna.report.geneset_ppi'
        task_type = 'workflow'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "taxon": data.taxon,
            "geneset_id": data.geneset_id,
            "geneset_type": data.geneset_type,
            "species": data.species,
            "combine_score": data.combine_score,
            "task_id": data.task_id,
        }
        geneset_info = self.lnc_rna.get_main_info(data.geneset_id, 'sg_geneset', data.task_id)
        if geneset_info['type'] != data.geneset_type:
            variables = []
            variables.append(geneset_info['type'])
            variables.append(data.geneset_type)
            info = {"success": False, "info": "基因集属性为%s 但筛选基因集属性为%s" % (geneset_info['type'], data.geneset_type)}
            return json.dumps(info)
        if not geneset_info:
            info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        typee = geneset_info['type']
        task_info = self.lnc_rna.get_task_info(geneset_info['task_id'])

        if task_info.has_key('assemble_fa'):
            seq = task_info['assemble_fa']
        elif task_info.has_key('params'):
            work_flow_params = dict(task_info['params'])
            fa_dir = work_flow_params['gene_fasta']
            fa_dir = fa_dir.split("||")[1]
            if re.match(r'tsanger:',fa_dir):
                seq = fa_dir.replace('tsanger:','/mnt/ilustre/tsanger-data/')
            else:
                seq = fa_dir.replace('sanger:','/mnt/ilustre/data/')
        else:
            info = {"success": False, "info": "找不到序列文件"}
            return json.dumps(info)

        annot_stat = self.lnc_rna.get_annotation_stat_info(geneset_info['task_id'])
        if not annot_stat:
            info = {"success": False, "info": "annot_stat不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        else:
            if annot_stat['has_new']:
                t2g_dir = annot_stat['result_dir'] + '/allannot_class/all_tran2gene.txt'
            else:
                t2g_dir = annot_stat['result_dir'] + '/refannot_class/all_tran2gene.txt'

        main_table_name = 'genesetPPI_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('geneset_type', typee),
            ('name', main_table_name),
            ('t2g_dir', t2g_dir),
            ('geneset_id', ObjectId(data.geneset_id)),
            ('desc', 'ppi_network分析中...'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.lnc_rna.insert_main_table('sg_geneset_ppi', mongo_data)
        update_info = {str(main_table_id): "sg_geneset_ppi"}

        options = {
            'update_info': json.dumps(update_info),
            "ppi_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "geneset_list": data.geneset_id,
            "species": data.species,
            "t2g_dir": self.use_s3(t2g_dir),
            "seq": self.use_s3(seq),
            "type": typee,
            "combine_score": data.combine_score
        }
        to_file = "lnc_rna.export_gene_list_ppi(geneset_list)"
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type=task_type,
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(GenesetPpiAction, self).POST()

        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        geneset_info = self.lnc_rna.insert_geneset_info(data.geneset_id, "sg_ppinetwork", str(main_table_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.lnc_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.lnc_rna.update_group_compare_is_use(data.task_id, data.control_id)

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
        cmd += "s/lnc_rna/geneset_ppi "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="lnc_rna",
            task_type="2",
            submit_location="genesetppi",
            geneset_id="5c73966117b2bf6ef67623e5",
            taxon = 'Animals',
            species= '9606',
            #species = '3702',
            combine_score = '300',
            geneset_type = 'G'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
