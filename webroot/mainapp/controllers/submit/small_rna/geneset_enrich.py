# -*- coding: utf-8 -*-
# __author__ = 'qindanhua, qinjincheng'

from mainapp.controllers.project.small_rna_controller import SmallRnaController
from mainapp.libs.signature import check_sig
import web
import json
import datetime
from bson.objectid import ObjectId
import unittest
import os
from biocluster.config import Config

class GenesetEnrichAction(SmallRnaController):
    def __init__(self):
        super(GenesetEnrichAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()

        # check args
        default_args = ['task_id', 'submit_location', 'task_type']
        # geneset_id ~ sg_geneset.main_id
        # type in ['origin', 'latest']
        # anno_type in ['go', 'kegg']
        # method in ['BH', 'BY', 'holm', 'bonferroni']
        default_args.extend(['geneset_id', 'type', 'anno_type', 'method'])
        for arg in default_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Lack argument: {}'.format(arg)}
                return json.dumps(info)

        # set variables
        task_id = data.task_id
        submit_location = data.submit_location
        task_type = int(data.task_type)
        geneset_id = data.geneset_id
        annotation = data.type
        anno_type = data.anno_type
        method = data.method

        # package incoming data as params
        params = {
            'task_id': task_id,
            'submit_location': submit_location,
            'task_type': task_type,
            'geneset_id': geneset_id,
            'type': annotation,
            'anno_type': anno_type,
            'method': method,
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        # check task_id and collect related variables
        gs_main_id = geneset_id
        geneset_info = self.small_rna.get_main_info(gs_main_id, 'sg_geneset', task_id)
        if not geneset_info:
            info = {'success': False, 'info': 'Geneset not found, task_id: {}'.format(task_id)}
            return json.dumps(info)
        if geneset_info["gene_length"] == 0:
                info = {'success': False, 'info': '该基因集为空，请选择其它的基因集分析'}
                return json.dumps(info)
        task_info = self.small_rna.get_task_info(geneset_info['task_id'])

        # check anno_type and define related variables
        if anno_type == 'go':
            table_name = 'Go'
            collection_name = 'sg_geneset_go_enrich'
            to_file = ['small_rna.export_gene_list(geneset_list)',
                       'small_rna.export_all_list(all_list)',
                       'small_rna.export_go_list(go_list)']
            option = {'go_list': geneset_id}
        elif anno_type == 'kegg':
            table_name = 'Kegg'
            collection_name = 'sg_geneset_kegg_enrich'
            to_file = ['small_rna.export_gene_list(geneset_list)',
                       'small_rna.export_all_list(all_list)',
                       'small_rna.export_kegg_table(kegg_table)',
                       'small_rna.export_multi_gene_list(geneset_kegg)',
                       'small_rna.export_add_info(add_info)']
            option = {'kegg_table': geneset_id,
                      'geneset_kegg': geneset_id,
                      'add_info': '{}\tG'.format(task_id),
                      'task_id': task_id}
        else:
            info = {'success': False, 'info': 'Unsupported functional classification: {}'.format(anno_type)}
            return json.dumps(info)

        # prepare mongo_data for inserting a new document in specific collection
        project_sn = task_info['project_sn']
        time_now = datetime.datetime.now()
        name = 'Geneset{}Enrich_{}'.format(table_name, time_now.strftime("%Y%m%d_%H%M%S"))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'Geneset {} enrich main table built at controller'.format(anno_type)
        mongo_data = [
            # essential keys
            ('task_id', task_id),
            ('project_sn', project_sn),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('name', name),
            ('desc', desc),
            ('params', params),
            # status of process
            ('status', 'start'),
        ]
        main_id = self.small_rna.insert_main_table(collection_name, mongo_data)

        # prepare sheet data for workflow
        options = {
            'geneset_id': geneset_id,
            'type': annotation,
            'anno_type': anno_type,
            'main_id': str(main_id),
            'geneset_type': 'G',
            'method': method,
            'update_info': json.dumps({str(main_id): collection_name}),
            'geneset_list': geneset_id,
            'all_list': geneset_id,
        }
        options.update(option)
        if "database_version" in task_info:
            kegg_version = task_info["database_version"].get("kegg", "")
            options.update({"kegg_version": kegg_version})
        elif "genome_id" in task_info:
            # 判断基因组注释版本
            database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
            collection = database['sg_genome_db']
            genome_doc = collection.find_one({'genome_id': task_info["genome_id"]})
            if "database_version" in genome_doc and "kegg" in genome_doc["database_version"]:
                options.update({
                    "kegg_version": genome_doc["database_version"]["kegg"]
                })
                go_version = genome_doc["database_version"].get("go", "20200628").split("_")[0]
                options.update({"go_version": go_version})

        # prepare sheet data for workflow
        task_name = 'small_rna.report.geneset_enrich'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,
                            module_type='workflow',
                            to_file=to_file,
                            project_sn=project_sn,
                            task_id=task_id)

        # run workflow and obtain return value
        task_info = super(GenesetEnrichAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': name}}

        # insert a document to sg_geneset_info after successfully running workflow
        self.small_rna.insert_geneset_info(geneset_id, collection_name, str(main_id))
        return json.dumps(task_info)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_go(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += '-fr no '
        cmd += '-c {} '.format('client03')
        cmd += 's/small_rna/geneset_enrich '
        cmd += '-b http://192.168.12.102:9090 '
        args = {
            'task_id': 'small_rna',
            'submit_location': 'genesetgoenrich',
            'task_type': '2',
            'geneset_id': '5bf23ee6a4e1af04b0877863',
            'type': 'origin',
            'anno_type': 'go',
            'method': 'BH',
        }
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)

    def test_kegg(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += '-fr no '
        cmd += '-c {} '.format('client03')
        cmd += 's/small_rna/geneset_enrich '
        cmd += '-b http://192.168.12.102:9090 '
        args = {
            'task_id': 'small_rna',
            'submit_location': 'genesetkeggenrich',
            'task_type': '2',
            'geneset_id': '5bf23ee6a4e1af04b0877863',
            'type': 'origin',
            'anno_type': 'kegg',
            'method': 'BH',
        }
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
