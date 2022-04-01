# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
import datetime
from biocluster.config import Config
from mainapp.models.mongo.denovo import Denovo
from mainapp.libs.param_pack import GetUploadInfo_denovo
from mbio.api.database.denovo_rna_mapping import *
from mainapp.models.mongo.submit.denovo_rna.denovo_mapping import DenovoMapping
import random


class MapAssessment(object):
    """
    比对质量评估接口
    """
    def __init__(self):
        self.db_name = Config().MONGODB + '_rna'

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if not hasattr(data, "express_id"):
            info = {"success": False, "info": "缺少参数express_id!"}
            return json.dumps(info)
        if not hasattr(data, "analysis_type"):
            info = {"success": False, "info": "缺少参数analysis_type!"}
            return json.dumps(info)
        analysis_type = data.analysis_type
        satur_params = ["low_bound", "up_bound", "step", "quality_satur", "orf_id"]
        coverage_params = ["orf_id", "min_len"]
        if analysis_type not in ["saturation", "coverage", "duplication", "correlation"]:
            info = {"success": False, "info": "不参在%s分析类型" % analysis_type}
            return json.dumps(info)
        if data.analysis_type == "saturation":
            for param in satur_params:
                if not hasattr(data, param):
                    info = {"success": False, "info": "缺少%s参数" % param}
                    return json.dumps(info)
        if data.analysis_type == "coverage":
            for param in coverage_params:
                if not hasattr(data, param):
                    info = {"success": False, "info": "缺少%s参数" % param}
                    return json.dumps(info)
        if data.analysis_type == "duplication":
            if not hasattr(data, "quality_dup"):
                info = {"success": False, "info": "缺少quality_dup参数"}
                return json.dumps(info)

        express_info = Denovo().get_main_info(data.express_id, "sg_denovo_express")
        if express_info:
            task_info = Denovo().get_task_info(express_info["task_id"])
            if task_info:
                member_id = task_info["member_id"]
            else:
                info = {"success": False, "info": "这个express_id对应的表达量矩阵对应的task：{}没有member_id!".format(express_info["task_id"])}
                return json.dumps(info)
            insert_data = self.get_insert_data(analysis_type, client, express_info, data, member_id)
            workflow_module = Workflow()
            workflow_module.add_record(insert_data)
            info = {"success": True, "info": "提交成功!"}
            return json.dumps(info)
        else:
            info = {"success": False, "info": "表达量表id不存在！!"}
            return json.dumps(info)

    def get_params(self, data):
        # print(data.express_id)
        my_param = {'analysis_type': data.analysis_type, "express_id": data.express_id}
        if data.analysis_type == "saturation":
            my_param["low_bound"] = int(data.low_bound)
            my_param["up_bound"] = int(data.up_bound)
            my_param["step"] = int(data.step)
            my_param["quality_satur"] = data.quality_satur
        # elif data.analysis_type == "coverage":
        #     my_param["orf_id"] = data.orf_id
        #     my_param["min_len"] = data.min_len
        elif data.analysis_type == "duplication":
            my_param["quality_dup"] = data.quality_dup
        # params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        return my_param

    def get_insert_data(self, analysis_type, client, express_info,  data, member_id):
        my_params = self.get_params(data)
        # params = json.dumps(my_params, sort_keys=True, separators=(',', ':'))
        params = my_params
        options = {'analysis_type': data.analysis_type, "bam": data.express_id, "express_id": data.express_id}
        to_file = ["denovo.export_bam_path(express_id)", "denovo.export_bed_file(orf_id)"]
        name = analysis_type + "_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        if analysis_type == "saturation":
            rpkm_id = DenovoMapping().add_rpkm_table(express_info["project_sn"], express_info["task_id"], params=params, name=name)
            update_info = {str(rpkm_id): "sg_denovo_rpkm", 'database': self.db_name}
            update_info = json.dumps(update_info)
            options["update_info"] = update_info
            options["insert_id"] = str(rpkm_id)
            options["bed"] = data.orf_id
            options.update(my_params)
            to_file = ["denovo.export_bam_path(bam)", "denovo.export_bed_path(bed)"]
        # elif analysis_type == "coverage":
        #     coverage_id = DenovoMapping().add_coverage_table(params=params, name=name, detail=False)
        #     update_info = {str(coverage_id): "sg_denovo_coverage", 'database': self.db_name}
        #     update_info = json.dumps(update_info)
        #     options["update_info"] = update_info
        #     options["insert_id"] = coverage_id
        #     options.update(my_params)
        elif analysis_type == "duplication":
            duplication_id = DenovoMapping().add_duplication_table(express_info["project_sn"], express_info["task_id"], params=params, name=name)
            update_info = {str(duplication_id): "sg_denovo_duplication", 'database': self.db_name}
            update_info = json.dumps(update_info)
            options["update_info"] = update_info
            options["insert_id"] = str(duplication_id)
            options.update(my_params)
            to_file = ["denovo.export_bam_path(bam)"]
        elif analysis_type == "correlation":
            correlation_id = DenovoMapping().add_correlation_table(express_info["project_sn"], express_info["task_id"], params=params, name=name)
            update_info = {str(correlation_id): "sg_denovo_correlation", 'database': self.db_name}
            update_info = json.dumps(update_info)
            options["update_info"] = update_info
            options["insert_id"] = str(correlation_id)
            options["fpkm"] = data.express_id
            options.update(my_params)
            to_file = ["denovo.export_express_matrix(fpkm)", "denovo.export_bam_path(bam)"]

        # workflow_id = Denovo().get_new_id(express_info["task_id"], data.express_id)
        workflow_id = self.get_new_id(express_info["task_id"], data.express_id)
        (output_dir, update_api) = GetUploadInfo_denovo(client, member_id, express_info['project_sn'], express_info['task_id'], 'gene_structure')
        json_data = {
            "id": workflow_id,
            "stage_id": 0,
            "name": "denovo_rna.report.map_assessment",
            "type": "workflow",
            "client": client,
            "project_sn": express_info["project_sn"],
            "to_file": to_file,
            "USE_DB": True,
            "IMPORT_REPORT_DATA": True,
            "UPDATE_STATUS_API": update_api,
            "IMPORT_REPORT_AFTER_END": True,
            "output": output_dir,
            "options": options
        }
        insert_data = {"client": client,
                       "workflow_id": workflow_id,
                       "json": json.dumps(json_data),
                       "ip": web.ctx.ip
                       }
        #print("lllllllllllll")
        #print(options)
        return insert_data

    def get_new_id(self, task_id, main_id):
        new_id = "%s_%s_%s" % (task_id, main_id[-4:], random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id, main_id)
        return new_id
