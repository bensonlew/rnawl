# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import os
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import param_pack, group_detail_sort, filter_json_sort
import datetime
from bson import SON
from biocluster.config import Config


class SampleCheckAction(MetaasvController):
    """
    metaasv 样本检测
    """
    def __init__(self):
        super(SampleCheckAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        self.check_option(data)
        task_name = 'metaasv.report.sample_extract'
        module_type = 'workflow'

        params_json = {}
        if hasattr(data, "fastq_type"):
            params_json["fastq_type"] = data.fastq_type
        if hasattr(data, "fastq_file"):
            params_json["fastq_file"] = data.fastq_file
        if hasattr(data, "submit_location"):
            params_json["submit_location"] = data.submit_location
        if hasattr(data, "task_type"):
            params_json["task_type"] = str(data.task_type)
        if hasattr(data, "task_id"):
            params_json["task_id"] = data.task_id
        if hasattr(data, "project_sn"):
            params_json["project_sn"] = data.project_sn
        if hasattr(data, "from_task"):
            params_json["from_task"] = data.from_task
        if hasattr(data, "query_id"):
            params_json["query_id"] = data.query_id
        if hasattr(data, "fastq_dir"):
            params_json["fastq_dir"] = json.loads(data.fastq_dir)
        if hasattr(data, "file_name"):
            params_json["file_name"] = data.file_name

        main_table_name = 'Sample_check_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        if hasattr(data, "project_sn"):
            mongo_data.append(('project_sn', data.project_sn))
        if hasattr(data, "task_id"):
            mongo_data.append(('task_id', data.task_id))
        if hasattr(data, "query_id"):
            mongo_data.append(('query_id', data.query_id))
        if hasattr(data, "from_task"):
            mongo_data.append(('from_task', data.from_task))
        else:
            mongo_data.append(('from_task', data.task_id))
        if hasattr(data, "file_name"):
            mongo_data.append(('file_name', data.file_name))

        main_table_id = self.metaasv.insert_main_table('sample_check', mongo_data)
        update_info = {str(main_table_id): 'sample_check'}
        options = {
            'update_info': json.dumps(update_info),
            "main_id": str(main_table_id),
            'main_table_data': SON(mongo_data)
        }
        to_file = []
        if hasattr(data, "fastq_type"):
            options["fastq_type"] = data.fastq_type
        if hasattr(data, "fastq_file"):
            options["fastq_file"] = data.fastq_file
        if hasattr(data, "fastq_dir"):
            # json_path = os.path.join(Config().WORK_DIR, "tmp", datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3] + ".fastq.json")
            # with open(json_path, 'w') as fw:
            #     fw.write(data.fastq_dir)
            # options["fastq_dir"] = json_path
            options["fastq_dir"] = data.fastq_dir
            to_file.append('metaasv.export_data_from_fastq_dir(fastq_dir)')
        if hasattr(data, "from_task"):
            if hasattr(data, "task_id"):
                if data.from_task != data.task_id:
                    options["from_task"] = data.from_task
        if hasattr(data, "query_id"):
            options["query_id"] = data.query_id
        if hasattr(data, "file_name"):
            options["file_name"] = data.file_name
        self.set_sheet_data(name=task_name, options=options, main_table_name="Sample_check/" + main_table_name,module_type=module_type, main_id=str(main_table_id),to_file=to_file,collection_name="sample_check")
        task_info = super(SampleCheckAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        if not task_info["success"]:
            task_info["info"] = "程序运行出错，请检查输入的文件夹和文件是否正确！"
        return json.dumps(task_info)

    def check_option(self, data):
        """
        对交互分析的样本检测进行参数检查
        :return:
        """
        """
        注释一下不同的情况:
        1. 分文件和文件夹两种类型进行样本检测；
        2. 分已经做过样本检测和没有做过
        """

        if hasattr(data, "from_task"):#做过样本检测，直接在原有的基础之上进行拷贝
            if not hasattr(data, "query_id"):
                info = {'success': False, 'info': 'parameters missing:%s'%("query_id")}
                return json.dumps(info)
            if not hasattr(data, "task_id"):
                info = {'success': False, 'info': 'parameters missing:%s'%("task_id")}
                return json.dumps(info)
        else:##没有做过样本检测，直接做
            if hasattr(data, "fastq_type"):
                if not hasattr(data, "query_id"):
                    info = {'success': False, 'info': 'parameters missing:%s'%("query_id")}
                    return json.dumps(info)
                if data.fastq_type == "fastq":
                    if not hasattr(data, "fastq_file"):
                        info = {'success': False, 'info': 'parameters missing:%s'%("fastq_file")}
                        return json.dumps(info)
                elif data.fastq_type == "fastq_dir":
                    if not hasattr(data, "fastq_dir"):
                        info = {'success': False, 'info': 'parameters missing:%s'%("fastq_dir")}
                        return json.dumps(info)
