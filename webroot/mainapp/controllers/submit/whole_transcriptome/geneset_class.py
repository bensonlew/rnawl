# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,fuwenyao'

import datetime
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig
from mbio.api.to_file.whole_transcriptome.geneset import *


class GenesetClassAction(WholeTranscriptomeController):
    def __init__(self):
        super(GenesetClassAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'level', 'submit_location', 'anno_type', 'task_type']
        for arg in default_argu:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "lack argument: %s" % arg, 'code': 'C2901004', 'variables': variables}
                return json.dumps(info)

        task_name = 'whole_transcriptome.report.geneset_class'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "anno_type": data.anno_type,
            "level": data.level,
            "task_id": data.task_id
        }
        # 判断传入的基因集id是否存在
        # print(data.geneset_id)
        geneset_info = {}
        for gd in data.geneset_id.split(","):
            geneset_info = self.whole_transcriptome.get_main_info(gd, 'geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset not found", 'code': 'C2901005', 'variables': ''}
                return json.dumps(info)
        task_info = self.whole_transcriptome.get_task_info(geneset_info['task_id'])

        if data.anno_type == "go":
            table_name = "GO"
            collection_name = "geneset_go_class"
            to_file = 'whole_transcriptome.geneset.export_go_class(geneset_go)'
            option = {"geneset_go": data.geneset_id}

        elif data.anno_type == "cog":
            table_name = "COG"
            collection_name = "geneset_cog_class"
            to_file = 'whole_transcriptome.geneset.export_cog_class(geneset_cog)'
            option = {"geneset_cog": data.geneset_id}

        else:
            info = {'success': False, 'info': '不支持的功能分类!', 'code': 'C2901006', 'variables': ''}
            return json.dumps(info)

        main_table_name = table_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        # print(main_table_name)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('version', 'v1')
        ]
        main_table_id = self.whole_transcriptome.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "geneset_type": data.level,
            "anno_type": data.anno_type,
        }
        options.update(option)
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(GenesetClassAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}
        geneset_info = self.whole_transcriptome.insert_geneset_info(data.geneset_id, collection_name,
                                                                    str(main_table_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.whole_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.whole_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/whole_transcriptome/geneset_class "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            submit_location="genesetcog",
            db_type="whole_transcriptome",
            geneset_id="5db032b317b2bf3cdfb10638",
            level="T",
            anno_type="cog",
            # geneset_kegg="5b07cc2517b2bf2d03d010e1",
            task_id="whole_transcriptome",
            task_type="2",
            # type="origin"

        )

        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
    # '''{'UPDATE_STATUS_API': 'ref_rna.tupdate_status', 'output': u'tsanger:rerewrweset/files/m_188/188_59f029eed6609/tsg_29390/interaction_results/GenesetCogClass_20180525_123935715', 'id': 'tsg_29390_23175_719852', 'stage_id': 0, 'interaction': True, 'name': 'ref_rna.report.geneset_class', 'db_type': 'ref_rna', 'client': u'client03', 'project_sn': u'188_59f029eed6609', 'IMPORT_REPORT_DATA': True, 'to_file': 'ref_rna.export_cog_class(geneset_cog)', 'type': 'workflow', 'options': {'anno_type': u'cog', 'task_type': u'', 'main_table_id': '5b079387a4e1af6b3ecf326c', 'submit_location': u'genesetcog', 'update_info': '{"5b079387a4e1af6b3ecf326c": "sg_geneset_cog_class"}', 'geneset_type': u'gene', 'geneset_id': u'5ae93962a4e1af7f1cd07dd8,5ae93964a4e1af7f1cd07ddc', 'geneset_cog': u'5ae93962a4e1af7f1cd07dd8,5ae93964a4e1af7f1cd07ddc'}}'''
