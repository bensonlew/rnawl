# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import web
import json
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from biocluster.wpm.client import *
from mainapp.libs.getip import get_ip
from mainapp.models.admin.userlog import UserlogModel
import datetime
import unittest
import os
from biocluster.wpm.client import worker_client


class StopAction(WholeTranscriptomeController):
    def __init__(self):
        super(StopAction, self).__init__(instant=True)
        self.input_data = web.input()
        print(self.input_data)

    @check_sig
    def POST(self):
        check = self.check_params()
        if check is not True:
            return check
        packed_params = self.pack_params()
        info = self.create_main_table(packed_params)
        print(info)
        if not info[0]:
            return info[1]
        else:
            main_id = info[1]
            run_id = info[2]
            workflow_module = Workflow()
            workflow_data = workflow_module.get_by_workflow_id(run_id)
            client = self.input_data.client if hasattr(
                self.input_data, "client") else web.ctx.env.get('HTTP_CLIENT')
            if len(workflow_data) > 0:
                print("往老框架发送任务{} stop信息".format(run_id))
                record = workflow_data[0]
                if client == record.client:
                    if record.has_run == 1 and (record.is_end or record.is_error):
                        if record.is_end:
                            if record.is_error:
                                self.modify_status('failed')
                                self.update_main('end', '任务已改为失败！', {'_id': ObjectId(main_id)})
                                info = {"success": True, "info": "任务已经失败！"}
                            else:
                                # self.modify_status('canceled')
                                self.update_main('end', '任务已经取消！', {'_id': ObjectId(main_id)})
                                info = {"success": True, "info": "任务已经取消！"}
                        else:
                            self.modify_status('failed')
                            self.update_main('end', '任务已改为失败！', {'_id': ObjectId(main_id)})
                            info = {"success": True, "info": "任务已经失败！"}
                        return json.dumps(info)
                    else:
                        insert_data = {"client": client,
                                       "ip": get_ip()
                                       }
                        if workflow_module.set_stop(run_id, insert_data):
                            session = web.config.get('_session')
                            if session and session.is_login:
                                model = UserlogModel()
                                model.add("停止运行Workflow: {}".format(run_id))
                            # self.modify_status('canceled')
                            self.update_main('end', '任务取消成功！', {'_id': ObjectId(main_id)})
                            info = {"success": True, "info": "任务取消成功！"}
                            return json.dumps(info)
                        else:
                            self.update_main('failed', '内部错误！', {'_id': ObjectId(main_id)})
                            info = {"success": False, "info": "内部错误！"}
                            return json.dumps(info)
                else:
                    self.update_main('failed', '没有权限操作！', {'_id': ObjectId(main_id)})
                    info = {"success": False, "info": "没有权限操作！"}
                    return json.dumps(info)
            else:
                try:
                    print("往新框架发送任务{} stop信息".format(run_id))
                    json_data = {"id": run_id.strip(), "msg": "stop"}
                    response = worker_client().run_cmd(json_data)
                    print(response)
                    if response['success']:
                        if 'end' in response['info']:
                            self.modify_status('end')
                            self.update_main('end', '任务已改为成功！', {'_id': ObjectId(main_id)})
                            info = {"success": False, "info": "任务实际已运行结束，已更新为正确状态！"}
                            return json.dumps(info)
                        elif 'error' in response['info']:
                            self.modify_status('failed')
                            self.update_main('end', '任务已改为失败！', {'_id': ObjectId(main_id)})
                            info = {"success": False, "info": "任务实际已运行失败，已更新为正确状态！"}
                            return json.dumps(info)
                        elif 'stop' in response['info']:
                            self.modify_status('terminated')
                            self.update_main('end', '任务已改为终止！', {'_id': ObjectId(main_id)})
                            info = {"success": True, "info": "任务取消成功！"}
                            return json.dumps(info)
                        else:
                            self.modify_status('terminated')
                            self.update_main('end', '任务取消成功！', {'_id': ObjectId(main_id)})
                            info = {"success": True, "info": "任务取消成功！"}
                            return json.dumps(info)
                    else:
                        if 'not found' in response['info']:
                            print("找不到任务ID，现取消任务！")
                            self.modify_status('terminated')
                            self.update_main('end', '任务取消成功！', {'_id': ObjectId(main_id)})
                        else:
                            return json.dumps(response)
                except:
                    self.update_main('failed', '任务ID不存在！', {'_id': ObjectId(main_id)})
                    info = {"success": False, "info": "stop: 任务ID不存在！"}
                    return json.dumps(info)

    def check_params(self):
        requires = ['main_id', 'submit_location', 'task_type']
        for i in requires:
            if not (hasattr(self.input_data, i)):
                variables = list()
                variables.append(i)
                return json.dumps({"success": False, "info": "缺少%s参数!" % i, 'code': 'C2900101', 'variables': variables})
        return True

    def pack_params(self):
        params_dict = dict()
        for each in self.input_data:
            if each == "task_type":
                params_dict[each] = int(self.input_data[each])
            else:
                params_dict[each] = self.input_data[each]
        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        return params

    def create_main_table(self, packed_params):
        status_info = self.whole_transcriptome.get_main_info_by_record('sg_status', _id=ObjectId(self.input_data.main_id))
        print(status_info)
        if status_info:
            if status_info['status'] == 'end':
                info = {"success": False, "info": "任务已经结束，无法取消！"}
                return [False, json.dumps(info)]
            elif status_info['status'] == 'failed':
                info = {"success": False, "info": "任务已经失败，无法取消！"}
                return [False, json.dumps(info)]
            elif status_info['status'] == 'deleted':
                info = {"success": False, "info": "任务已经删除，无法取消！"}
                return [False, json.dumps(info)]
            elif status_info['status'] == 'terminated':
                info = {"success": False, "info": "任务已经取消，请勿重复操作！"}
                return [False, json.dumps(info)]
            else:
                if 'run_id' in status_info and status_info['run_id']:
                    table_name = status_info['table_name']
                    run_id = status_info['run_id']
                    name = "Cancel_" + table_name + "_"
                    time_now = datetime.datetime.now()
                    name += time_now.strftime("%Y%m%d_%H%M%S")
                    main_info = dict(
                        task_id=status_info['task_id'],
                        name=name,
                        created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                        desc='cancel interaction task',
                        params=packed_params,
                        status="start")
                    main_id = self.whole_transcriptome.insert_main_table('sg_cancel', main_info)
                    return [True, main_id, run_id]
                else:
                    info = {"success": False, "info": "找不到任务运行ID，无法取消！"}
                    return [False, json.dumps(info)]
        else:
            info = {"success": False, "info": "mongo数据库中找不到记录！"}
            return [False, json.dumps(info)]

    def modify_status(self, status):
        # 更改sg_status状态
        query_dict = {'_id': ObjectId(self.input_data.main_id)}
        time_now = datetime.datetime.now()
        modified_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        update_dict = {'status': status, 'modified_ts': modified_ts}
        self.whole_transcriptome.update_db_record1('sg_status', query_dict=query_dict, insert_dict=update_dict)
        # 更改主表状态
        status_info = self.whole_transcriptome.get_main_info_by_record('sg_status', _id=ObjectId(self.input_data.main_id))
        table_name = status_info['type_name']
        table_id = status_info['table_id']
        self.whole_transcriptome.update_db_record1(table_name, query_dict={'_id': table_id}, insert_dict=update_dict)

    def update_main(self, status, desc, query_dict):
        update_dict = {'status': status, 'desc': desc}
        self.whole_transcriptome.update_db_record1('sg_cancel', query_dict=query_dict, insert_dict=update_dict)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/whole_transcriptome/stop "
        cmd += "-b http://bcl.tsg.com "
        cmd += "-dbversion 1 "
        args = dict(
            main_id="60b5f19917b2bf141970b5a2",
            submit_location='task_cancel',
            task_type='1',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
