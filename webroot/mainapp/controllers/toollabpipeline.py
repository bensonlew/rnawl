# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20200402

import os
import re
import web
import json
import datetime
from bson import SON
# from biocluster.config import Config
from mainapp.config.db import Config
from mainapp.libs.toollabsignature import check_sig
from mainapp.libs.toollab_input_check import check_format
from mainapp.models.mongo.toollab import Toollab
from mainapp.controllers.project.toollab_controller import ToollabController
import importlib


class ToollabpipelineAction(ToollabController):
    """
    小工具的通用接口，负责接收前端传入的数据，并将任务分发投递到WPM服务
    """
    def __init__(self):
        super(ToollabpipelineAction, self).__init__()

    @check_sig
    @check_format
    def POST(self):
        """
        main_table_name主表名字 最好为sg_pca，这边会统一截取到pca用于main_table_name中命令
        params:中存入的是小工具分析需要的参数，是字典
        basis：中存入的是不是小工具分析需要的参数，比如task_id等
        :return:
        """
        data = web.input()
        print data
        basis = json.loads(data.basis)
        # 加载小工具参数检查模块line 39-44, 注意这里参数检查的文件名称，必须要与tool or module or workflow的名字一致
        if self.check_file_exists('toollabcheck/{}.py'.format(basis['name'].split('.')[-1])):
            module_name = ".".join(self.__module__.split('.')[:-1]) + '.toollabcheck.{}'\
                .format(basis['name'].split('.')[-1])
            tool_params_check = importlib.import_module(module_name)
            if tool_params_check.params_check(json.loads(data.params)):
                return tool_params_check.params_check(json.loads(data.params))
        else:
            print '没有检测到参数检查文件，不进行参数检查!'
        params_json = {
            "task_type": basis["task_type"],
            "task_id": basis["task_id"],
            "project_sn": basis["project_sn"]
        }
        for key in json.loads(data.params):
            params_json[key] = json.loads(data.params)[key]
        instant = True if basis["task_type"] == 'instant' else False
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        analysis_table_name = '{}_'.format(self.format_table_name(basis["main_table_name"])) + \
                              datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ("name", analysis_table_name),
            ("status", "start"),
            ("task_id", basis["task_id"]),
            ('project_sn', basis["project_sn"]),
            ("params", params),
            ("desc", "小工具主表！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('instant', instant)
        ]
        main_id = Toollab('tool_lab').insert_none_table(collection=basis["main_table_name"])  # 这里只获取_id, 主表是在
        # 请求WPM成功后，再插入
        Toollab('tool_lab').update_db_record(collection=basis["main_table_name"], query_dict={"_id": main_id},
                                             update_dict={"main_id": main_id})
        update_info = {str(main_id): basis["main_table_name"]}
        # add by hd 20200508 增加文件类型的判断，并提供下载文件时修改文件名称
        options = {}
        mapping_file = ""
        cal_params = json.loads(data.params)
        print "cal_params:{}".format(cal_params)
        if "mapping_file" in cal_params.keys():
            mapping_file = cal_params["mapping_file"]["path"]
        for key in cal_params:
            if key == "mapping_file":
                continue
            if isinstance(cal_params[key], dict):
                if "dir" in cal_params[key].keys() and cal_params[key]['dir']:
                    if not mapping_file:
                        raise Exception("文件夹参数未找到mapping file!")
                    else:
                        file_path_type = 'sanger' if data.client == "client01" else "tsanger"
                        options[key] = "filelist[%s](%s):%s" % (key, file_path_type, mapping_file)
                        if u"alias" in cal_params[key].keys() and cal_params[key]['alias']:
                            options[key] += "{%s}" % cal_params[key]['alias']
                        if u"format" in cal_params[key].keys() and cal_params[key]["format"]:
                            options[key] = "%s||%s" % (cal_params[key]["format"], options[key])
                else:
                    if u'path' in cal_params[key].keys() and u'alias' in cal_params[key].keys():
                        options[key] = cal_params[key]['path'] + '{' + cal_params[key]['alias'] + '}'
                    else:
                        options[key] = cal_params[key]
            else:
                options[key] = cal_params[key]
        options.update({
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            'main_table_data': SON(mongo_data)
        })
        self.set_sheet_data(name=basis["name"], task_id=basis["task_id"], main_table_name=analysis_table_name,
                            options=options, params=params, project_sn=basis['project_sn'], db_type='tool_lab',
                            instant=instant, module_type=basis["type"], client=data.client)
        task_info = super(ToollabpipelineAction, self).POST()
        try:
            if not task_info['success']:
                temp_info = json.loads(task_info['info'])
                if u'info' in temp_info.keys():
                    task_info['info'] = temp_info[u'info']
        except:
            pass
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def format_table_name(self, main_table_name):
        if re.match('sg_.*', main_table_name):
            return main_table_name[3:]
        else:
            return main_table_name

    def find_client(self):
        # client = Config().rcf.get("toollab", "client")
        client = Config().get_webauth("tool_lab", "client")
        if client not in ['client01', 'client03']:
            raise Exception('请检查config文件，client必须为client01 or clinet03')
        else:
            return client

    def check_file_exists(self, name):
        check_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
        if os.path.exists(check_file):
            return True
        else:
            return False
