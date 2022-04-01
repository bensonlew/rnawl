# -*- coding: utf-8 -*-
# __author__ = 'xuting'  last modify by qindanhua 20170110
import web
import json
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId
from bson import SON


class VennOldAction(MetaController):

    def __init__(self):
        super(VennOldAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        postArgs = ['group_id', 'level_id', "group_detail", 'venn_id']
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'parameters missing!'}
                return json.dumps(info)

        task_name = 'meta.report.venn'
        task_type = 'workflow'
        main_table_name = 'Venn_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        group_detal_dict = json.loads(data.group_detail)
        specimen_ids = list()
        for v in group_detal_dict.values():
            for tmp in v:
                specimen_ids.append(tmp)
        specimen_ids = ",".join(specimen_ids)
        options = {
            "in_otu_table": data.otu_id,
            "group_detail": data.group_detail,
            "group_table": data.group_id,
            "samples": self.meta.sampleIdToName(specimen_ids),
            "level": data.level_id,
            "otu_id": str(data.otu_id),
            "venn_id": str(data.venn_id),
            "old_venn_id": data.venn_id
        }
        to_file = ["meta.export_otu_table_by_level(in_otu_table)", "meta.export_group_table_by_detail(group_table)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="Venn/" + main_table_name,
                            module_type=task_type, to_file=to_file)
        task_info = super(VennOldAction, self).POST()
        if task_info['success']:
            task_info['content'] = {
                'ids': {
                    'id': str(data.venn_id),
                    'name': main_table_name
                }}
        return json.dumps(task_info)
