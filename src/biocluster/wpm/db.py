# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from ..core.function import hostname
from .mysql import Mysql
import json
from ..core.function import CJsonEncoder
import traceback
import sys
import time


class WorkflowModel(object):
    """
    操作数据库workflow表
    """
    def __init__(self, wsheet):
        """

        :param wsheet: sheet对象
        """
        self._db = Mysql()
        self.workflow_id = wsheet.id
        self.sheet = wsheet

    def __del__(self):
        self._db.close()

    def close(self):
        self._db.close()

    def save(self, pid=0, report=None):
        """
        添加workflow记录到表格中
        """
        try:
            is_instant = 1 if self.sheet.instant else 0
            if report is None:
                sql = "INSERT INTO workflow (client, workflow_id, json, server, pid, instant, path, type, batch_id) " \
                      "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)"

                data = (self.sheet.client, self.workflow_id, json.dumps(self.sheet.data, cls=CJsonEncoder),
                        hostname, pid, str(is_instant), self.sheet.name, self.sheet.type, self.sheet.batch_id)
            else:
                sql = "INSERT INTO workflow (client, workflow_id, json, server, pid, instant, path, type, batch_id," \
                      " is_report, report_type) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"

                data = (self.sheet.client, self.workflow_id, json.dumps(self.sheet.data, cls=CJsonEncoder),
                        hostname, pid, str(is_instant), self.sheet.name, self.sheet.type, self.sheet.batch_id,
                        1, report)

            count = self._db.insert_one(sql, data)
            if self.sheet.batch_id:
                sql = "update workflow set batch=1 where workflow_id = %s"
                self._db.update(sql, (self.sheet.batch_id, ))
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()
        else:
            return count

    def update_pid(self, pid, work_dir=""):
        sql = "update workflow set has_run=1, run_time=CURRENT_TIMESTAMP(), pid=%s, work_dir=%s " \
              "where workflow_id = %s"
        # print sql
        return self._db.update(sql, (pid, work_dir, self.workflow_id))

    def update_rerun(self, pid):
        sql = "update workflow set has_run=1, is_end=0,is_error=0, pid=%s, rerun=1, rerun_time=CURRENT_TIMESTAMP(), " \
              "server=%s where workflow_id = %s"
        return self._db.update(sql, (pid, hostname, self.workflow_id))

    def find(self):
        sql = "select json,add_time,is_end,has_run,run_time,work_dir from workflow where workflow_id = %s"
        # print sql
        return self._db.get_one(sql, (self.workflow_id, ))

    def update(self):
        sql = "update workflow set has_run=1, last_update =CURRENT_TIMESTAMP() " \
              "where workflow_id = %s"
        # print sql
        return self._db.update(sql, (self.workflow_id, ))

    def error(self, error_msg):
        sql = "update workflow set has_run=1, is_end=1, is_error=1, error=%s, end_time=CURRENT_TIMESTAMP() " \
              "where workflow_id = %s"
        data = (error_msg, self.workflow_id)
        # print sql
        return self._db.update(sql, data)

    def end(self):
        sql = "update workflow set has_run=1, is_end=1, is_error=0, error=' ', end_time=CURRENT_TIMESTAMP() " \
              "where workflow_id = %s"
        # print sql
        return self._db.update(sql, (self.workflow_id, ))

    def stop(self):
        sql = "update tostop set done=1, stoptime=CURRENT_TIMESTAMP() " \
              "where done=0 and workflow_id = %s"
        # print sql
        return self._db.update(sql, (self.workflow_id, ))

    def pause(self):
        try:
            self._db.cursor.execute("SET AUTOCOMMIT = 0")
            sql1 = "update pause set has_pause=1, pause_time=CURRENT_TIMESTAMP() " \
                   "where workflow_id = %s and has_pause=0"
            self._db.query(sql1, (self.workflow_id,))
            sql2 = "update workflow set paused=1 where workflow_id = %s"
            self._db.query(sql2, (self.workflow_id,))
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()
            self._db.end(option="rollback")
        else:
            self._db.end()
        self._db.cursor.execute("SET AUTOCOMMIT = 1")

    def exit_pause(self):
        try:
            self._db.cursor.execute("SET AUTOCOMMIT = 0")
            sql1 = "update pause set has_continue=1,continue_time=CURRENT_TIMESTAMP() " \
                   "where workflow_id = %s and has_pause=1 and exit_pause=1 and has_continue=0"
            self._db.query(sql1, (self.workflow_id, ))
            sql2 = "update workflow set paused=0 where workflow_id = %s"
            self._db.query(sql2, (self.workflow_id, ))
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()
            self._db.end(option="rollback")
        else:
            self._db.end()
        self._db.cursor.execute("SET AUTOCOMMIT = 1")

    def pause_timeout(self):
        try:
            self._db.cursor.execute("SET AUTOCOMMIT = 0")
            sql1 = "update pause set timeout=1,timeout_time=CURRENT_TIMESTAMP() " \
                   "where workflow_id = %s and has_pause=1 and exit_pause=0 and timeout=0"
            self._db.query(sql1, (self.workflow_id,))
            sql2 = "update workflow set paused=0 where workflow_id = %s"
            self._db.query(sql2, (self.workflow_id,))
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()
            self._db.end(option="rollback")
        else:
            self._db.end()
        self._db.cursor.execute("SET AUTOCOMMIT = 1")


class CheckModel(object):
    def __init__(self):
        self._db = Mysql()

    def __del__(self):
        self._db.close()

    def close(self):
        self._db.close()

    def find_stop(self):
        sql = "select workflow_id from tostop where done=0 and time > DATE_SUB(now(),INTERVAL 1 hour)"
        return self._db.get_all(sql)

    def find_pause(self):
        sql = "select workflow_id from pause where has_pause=0 and add_time > DATE_SUB(now(),INTERVAL 1 hour)"
        # print sql
        return self._db.get_all(sql)

    def find_exit_pause(self):
        sql = "select workflow_id from pause where has_pause=1 and exit_pause=1 and has_continue=0 and timeout=0 and " \
              "exit_pause_time > DATE_SUB(now(),INTERVAL 1 hour)"
        # print sql
        return self._db.get_all(sql)

    def update_running(self, running_list):
        if running_list:
            rl = ["\'%s\'" % n for n in running_list]
            ids = ", ".join(rl)
            sql = "update workflow set is_end=0,is_error=0, last_update = CURRENT_TIMESTAMP() where has_run=1 and " \
                  "is_end=0 and workflow_id in (%s)" % ids
            return self._db.update(sql)

    def update_unknown(self):
        sql = "update workflow set is_end=1,is_error=1, error='Unknown', end_time = CURRENT_TIMESTAMP() " \
              "where is_end=0 and is_error=0 and last_update < DATE_SUB(now(),INTERVAL 300 MINUTE) and " \
              "(rerun=0 or (rerun=1 and rerun_time < DATE_SUB(now(),INTERVAL 300 MINUTE)))"
        return self._db.update(sql)


class ApiLogModel(object):
    def __init__(self, log_object=None):
        self._db = Mysql()
        self.log_object = log_object

    def __del__(self):
        self._db.close()

    def close(self):
        self._db.close()

    def save(self):
        is_success = 1 if self.log_object.web_api_success else 0
        has_update_status = 1 if self.log_object.has_update_status else 0
        update_status_success = 1 if self.log_object.update_status_success else 0
        has_update_webapi = 1 if self.log_object.has_update_webapi else 0

        sql = "INSERT INTO apilog (task_id, api, data, update_status, update_status_success, webapi, " \
              "success, server, response, response_code, update_info) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, " \
              "%s, %s)"
        data = (self.log_object.task_id, self.log_object.api, json.dumps(self.log_object.data, cls=CJsonEncoder),
                has_update_status, update_status_success, has_update_webapi, is_success, hostname,
                self.log_object.response, self.log_object.response_code, json.dumps(self.log_object.update_info))
        return self._db.insert_one(sql, data)

    def find(self, id_list=(), start_time=None, end_time=None, host=None, api=None, webapi_only=False,
             update_status_only=False, failed_only=True):
        sql = "select * from apilog where "
        where_str = []
        data = []
        if id_list:
            where_str.append("task_id in (" + ",".join(["\'%s\'" % x for x in id_list]) + ")")
        if start_time:
            where_str.append("addtime >= %s")
            data.append(start_time)
        if end_time:
            where_str.append("addtime <= %s")
            data.append(end_time)
        if host:
            where_str.append("server = %s")
            data.append(host)
        if api:
            where_str.append("api = %s")
            data.append(api)
        if webapi_only:
            if failed_only:
                where_str.append("webapi = 1 and success = 0")
            else:
                where_str.append("webapi = 1")
        if update_status_only:
            if failed_only:
                where_str.append("update_status = 1 and update_status_success = 0")
            else:
                where_str.append("update_status = 1")
        sql += " and ".join(where_str)
        # print sql
        return self._db.get_all(sql, data)


class ClientKeyModel(object):
    def __init__(self):
        self._db = Mysql()

    def __del__(self):
        self._db.close()

    def close(self):
        self._db.close()

    def find_key(self, client):
        sql = "select key from clientkey where client = %s"
        data = self._db.get_one(sql, (client, ))
        if data:
            return data["key"]
        else:
            return None


class ReportModel(object):

    def __init__(self, wid):
        self._db = Mysql()
        self.workflow_id = wid
        # self.auto_id = self.get_workflow_id()

    # def get_workflow_id(self):
    #     sql = "select id from workflow where workflow_id = %s"
    #     result = self._db.get_one(sql, (self.workflow_id,))
    #     if result and "id" in result.keys():
    #         return result["id"]
    #     else:
    #         return 0
    def __del__(self):
        self._db.close()

    def close(self):
        self._db.close()

    def save_workflow(self, cpu, memory, error_info="", rerun=False):
        if rerun:
            sql = "update workflow set cpu_used = %s,memory_used= %s, error_data=%s " \
                  "where workflow_id = %s"
        else:
            sql = "update workflow set cpu_used=%s,memory_used=%s,error_data=%s where workflow_id = %s"
        return self._db.update(sql, (cpu, memory, self.workflow_id, error_info))

    def save_modules(self, data, workflow_id=None):
        if not workflow_id:
            workflow_id = self.workflow_id
        try:
            self._db.cursor.execute("SET AUTOCOMMIT = 0")
            for d in data:
                sql = "insert into module (parent_run_id, run_id, path, work_dir, start_time, end_time, tool_num) " \
                      "values (%s, %s, %s, %s, %s, %s, %s)"
                data = (workflow_id, d["run_id"], d["path"], d["work_dir"],
                        time.strftime("%Y-%m-%d %X", time.localtime(d["start_time"])),
                        time.strftime("%Y-%m-%d %X", time.localtime(d["end_time"])),
                        d["tool_num"])
                self._db.cursor.execute(sql, data)
                # self._db.cursor.execute("SELECT @@IDENTITY AS id")
                # result = self._db.cursor.fetchall()
                # wid = result[0]['id']
                if len(d["tools"]) > 0:
                    self.save_tools(d["run_id"], d["tools"])
                if "modules" in d.keys() and len(d["modules"]) > 0:
                    self.save_modules(d["modules"], d["run_id"])
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()
            self._db.end(option="rollback")
        else:
            self._db.end()
        self._db.cursor.execute("SET AUTOCOMMIT = 1")

    def save_tools(self, parent_run_id, data, commit=False, under_workflow=0):
        try:
            if not commit:
                self._db.cursor.execute("SET AUTOCOMMIT = 0")
            for d in data:
                sql = "insert into tool (parent_run_id, run_id, path, work_dir, run_host, job_type, job_id, " \
                      "request_cpu, request_memory, start_time, wait_spend_time, queue_spend_time,  run_spend_time, " \
                      "end_time, run_times, success, info, under_workflow, process_id, sub_process_num, max_cpu_use, " \
                      "max_rss, average_cpu_use, average_rss, max_vms, average_vms, params) values " \
                      "(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, " \
                      "%s, %s, %s, %s, %s)"
                data = (parent_run_id, d["run_id"], d["path"], d["work_dir"], d["run_host"], d["job_type"], d["job_id"],
                        d["request_cpu"], d["request_memory"],
                        time.strftime("%Y-%m-%d %X", time.localtime(d["start_time"])),
                        d["wait_spend_time"], d["queue_spend_time"], d["run_spend_time"],
                        time.strftime("%Y-%m-%d %X", time.localtime(d["end_time"])),
                        d["run_times"], d["success"], d["info"], under_workflow, d["process_id"], d["sub_process_num"],
                        d["max_cpu_use"], d["max_rss"], d["average_cpu_use"], d["average_rss"], d["max_vms"],
                        d["average_vms"], json.dumps(d["params"]))
                self._db.cursor.execute(sql, data)
                # self._db.cursor.execute("SELECT @@IDENTITY AS id")
                # result = self._db.cursor.fetchall()
                # tid = result[0]['id']
                if len(d["commands"]) > 0:
                    self.save_commands(d["run_id"], d["commands"])
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()
            if commit:
                self._db.end(option="rollback")
        else:
            if commit:
                self._db.end()
        self._db.cursor.execute("SET AUTOCOMMIT = 1")

    def save_commands(self, parent_run_id, data):
        for d in data:
            sql = "insert into command (parent_run_id, name, cmd, start_time, end_time, run_times, main_pid, " \
                  "sub_process_num, max_cpu_use, max_rss, average_cpu_use, average_rss, return_code, max_vms, " \
                  "average_vms) values " \
                  "(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
            data = (parent_run_id, d["name"], d["cmd"], time.strftime("%Y-%m-%d %X", time.localtime(d["start_time"])),
                    time.strftime("%Y-%m-%d %X", time.localtime(d["end_time"])),
                    d["run_times"], d["main_pid"], d["sub_process_num"],  d["max_cpu_use"],
                    d["max_rss"], d["average_cpu_use"], d["average_rss"], d["return_code"], d["max_vms"],
                    d["average_vms"])
            self._db.cursor.execute(sql, data)
