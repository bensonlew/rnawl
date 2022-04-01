# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

import web
import json
import types
import datetime
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from biocluster.config import Config
from mainapp.models.mongo.submit.denovo_rna.denovo_go_enrich import DenovoEnrich
from mainapp.models.mongo.meta import Meta
from mainapp.models.workflow import Workflow
from mainapp.controllers.project.denovo_controller import DenovoController


class GoEnrichRegulate(DenovoController):
    """
    go富集分析、go调控分析接口
    """
    def __init__(self):
        super(GoEnrichRegulate, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, 'client') else web.ctx.env.get('HTTP_CLIENT')
        #print data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": "+".json(return_result)}
            return json.dumps(info)
        my_param = {'analysis_type': data.analysis_type, "express_diff_id": data.express_diff_id, "compare": data.compare, "submit_location": data.submit_location}
        if data.analysis_type in ["enrich", "both"]:
            my_param["pval"] = float(data.pval)
            my_param["method"] = data.method
            my_param["regulate"] = data.regulate
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        name = data.compare.split(',')[0]
        compare_name = data.compare.split(',')[1]
        express_diff_info = Meta(db=self.mongodb).get_main_info(data.express_diff_id, "sg_denovo_express_diff")
        if not express_diff_info:
            info = {"success": False, "info": "express_diff_id不存在，请检查参数是否正确！"}
            return json.dumps(info)
        task_id = express_diff_info["task_id"]
        project_sn = express_diff_info["project_sn"]
        if data.analysis_type == "enrich":
            main_table_name = "GoRich_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            go_enrich_id = DenovoEnrich().add_go_enrich(name=main_table_name, params=params, project_sn=project_sn, task_id=task_id, go_graph_dir=None)
            update_info = {str(go_enrich_id): "sg_denovo_go_enrich"}
            options = {
                "analysis_type": data.analysis_type,
                "diff_stat": data.express_diff_id,
                "update_info": json.dumps(update_info),
                "go_enrich_id": str(go_enrich_id),
                "all_list": data.express_diff_id,
                "go_list": data.express_diff_id,
                "name": name,
                "compare_name": compare_name,
                "regulate": data.regulate
            }
            to_file = ["denovo.export_diff_express(diff_stat)", "denovo.export_all_gene_list(all_list)", "denovo.export_gos_list(go_list)"]
            self.set_sheet_data(name="denovo_rna.report.go_enrich_regulate", options=options, main_table_name=main_table_name, module_type="workflow", to_file=to_file, main_id=go_enrich_id, collection_name="sg_denovo_go_enrich")
        if data.analysis_type == "regulate":
            main_table_name = "GoRegulate_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            go_regulate_id = DenovoEnrich().add_go_regulate(name=main_table_name, params=params, project_sn=project_sn, task_id=task_id)
            update_info = {str(go_regulate_id): "sg_denovo_go_regulate"}
            options = {
                "analysis_type": data.analysis_type,
                "diff_stat": data.express_diff_id,
                "update_info": json.dumps(update_info),
                "go_regulate_id": str(go_regulate_id),
                "name": name,
                "compare_name": compare_name,
                "go2level": data.express_diff_id
            }
            to_file = ["denovo.export_diff_express(diff_stat)", "denovo.export_go2level(go2level)"]
            self.set_sheet_data(name="denovo_rna.report.go_enrich_regulate", options=options,
                                main_table_name=main_table_name, module_type="workflow", to_file=to_file,
                                main_id=go_regulate_id, collection_name="sg_denovo_go_regulate")
        if data.analysis_type == "both":
            main_table_name = "GoEnrichRegulate_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            sort_id = DenovoEnrich().add_go_pval_sort(name=main_table_name, params=params, project_sn=project_sn, task_id=task_id)
            update_info = {str(sort_id): "sg_denovo_go_pvalue"}
            options = {
                "analysis_type": data.analysis_type,
                "update_info": json.dumps(update_info),
                "sort_id": str(sort_id),
                "name": name,
                "compare_name": compare_name,
                "regulate": data.regulate,
                "diff_stat": data.express_diff_id,
                "all_list": data.express_diff_id,
                "go_list": data.express_diff_id,
                "go2level": data.express_diff_id
            }
            to_file = ["denovo.export_diff_express(diff_stat)", "denovo.export_all_gene_list(all_list)", "denovo.export_gos_list(go_list)", "denovo.export_go2level(go2level)"]
            self.set_sheet_data(name="denovo_rna.report.go_enrich_regulate", options=options,
                                main_table_name=main_table_name, module_type="workflow", to_file=to_file,
                                main_id=sort_id, collection_name="sg_denovo_go_pvalue")
        task_info = super(GoEnrichRegulate, self).POST()
        return json.dumps(task_info)


    def check_options(self, data):
        """
        检查网页端传来的参数是否正确
        """
        params_name = ["analysis_type", "express_diff_id", "compare", "submit_location"]
        success = []
        for name in params_name:
            if not hasattr(data, name):
                success.append("缺少参数：%" % name)
        express_diff_id = str(data.express_diff_id)
        if not isinstance(express_diff_id, ObjectId) and not isinstance(express_diff_id, types.StringType):
            success.append("传入的express_diff_id:%不是一个ObjectId对象或字符串类型!" % express_diff_id)
        analysis_type = data.analysis_type
        if analysis_type not in ["enrich", "regulate", "both"]:
            success.append("%分析不存在" % analysis_type)
        if analysis_type in ["enrich", "both"]:
            for name in ["pval", "method", "regulate"]:
                if not hasattr(data, name):
                    success.append("缺少参数:%" % name)
        return success
