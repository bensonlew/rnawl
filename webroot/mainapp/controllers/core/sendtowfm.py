# -*- coding: utf-8 -*-
# __author__ = 'HD'
# lasted modified @ 20210517

import sys
import gevent
import grpc
from biocluster.config import Config
import os
import json
import re


class Sendtowfm(object):
    def __init__(self, data):
        """
        需要导入"/mnt/ilustre/users/sanger-dev/bc2/sanger_bioinfo/src/biocluster/proto"
        用于将任务自动切换到新旧计算框架中去运行。
        配置文件，读取controllers.run_on_newbiocluster.txt
        # 有参编译检测V2
        wgs_v2.wgs_v2
        # wgs_v2.report中所有的交互分析投递到新平台
        regexp(^wgs_v2.report\\.\\w+$)
        # 单独一个交互分析
        wgs_v2.report.cnv_compare
        """
        self.json_data = data
        import_path = os.path.join(os.path.dirname(Config().WORK_DIR), "wpm2/sanger_bioinfo/src/biocluster/proto")
        import sys
        sys.path.append(import_path)
        self.wfm_port = "7321"
        self.instant = False
        self.timeout = 600
        self._cmdreq_times = 0

    def run_task(self):
        import webapi_pb2
        import webapi_pb2_grpc
        if not isinstance(self.json_data, dict) or "id" not in self.json_data.keys():
            print "add workflow %s format error!" % self.json_data
            return {"success": False, "info": "json格式错误"}
        client = "client03"
        if "instant" in self.json_data.keys():
            self.instant = self.json_data["instant"]
        if "client" in self.json_data.keys():
            client = self.json_data['client']
        self.json_data["WPM"] = True
        try:
            with grpc.insecure_channel('localhost:%s' % self.wfm_port) as channel:
                stub = webapi_pb2_grpc.WebApiStub(channel)
                response = stub.Submit(webapi_pb2.Task(
                    id=self.json_data["id"],
                    client=client,
                    json=json.dumps(self.json_data),
                    instant=self.instant
                ), timeout=self.timeout)
        except Exception as e:
            print "GRPC请求出错:{}".format(e)
            sys.stdout.flush()
            return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
        else:
            return {"success": response.ok, "info": response.reason}

    def run_cmd(self):
        import webapi_pb2
        import webapi_pb2_grpc
        self._cmdreq_times += 1
        if not isinstance(self.json_data, dict) or "id" not in self.json_data.keys():
            print "json格式不正确-样例：{'id': tsg_0001, 'msg': 'stop'}" % self.json_data
            return {"success": False, "info": "json格式错误"}
        try:
            with grpc.insecure_channel('localhost:%s' % self.wfm_port) as channel:
                stub = webapi_pb2_grpc.WebApiStub(channel)
                response = stub.Control(webapi_pb2.Command(
                    id=self.json_data["id"],
                    msg=self.json_data["msg"],
                ), timeout=self.timeout)
                print response
        except Exception as e:
            if self._cmdreq_times > 2:
                print "GRPC请求出错:{}".format(e)
                sys.stdout.flush()
                return {"success": False, "info": "GRPC请求出错,稍后再试:{}".format(e)}
            else:
                print("GRPC请求错误10秒后重试: %s" % e)
                sys.stdout.flush()
                gevent.sleep(5)
                self.run_cmd()
        else:
            return {"success": response.ok, "info": response.reason}


def check_target_biocuster(analysisname, confpath=None):
    """
    检查配置文件，确定要往那个平台上面投递任务。
    :param analysisname: 分析的工作流名字
    :param confpath: 要运行在新框架的流程
    :return:
    """
    runinnewbiocluster = False
    projects = []
    file_path = confpath
    if confpath is None:
        file_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                 "run_on_newbiocluster.txt")
    print "filepath:{}".format(file_path)
    with open(file_path, "r") as r:
        data = r.readlines()
        for line in data:
            if not re.match("^#.*", line):
                line = line.strip().split('\t')
                projects.append(line[0].strip())
    print "project:{}".format(projects)
    for m in projects:
        res = re.match(r"^regexp\((.+)\)$", m)
        if res:
            exp = res.group(1)
            res1 = re.match(exp, analysisname)
            if res1:
                runinnewbiocluster = True
        else:
            if m == analysisname:
                runinnewbiocluster = True
    return runinnewbiocluster


if __name__ == '__main__':
    res = check_target_biocuster("wgs_v2.wgs_v2")
    print res
    res1 = check_target_biocuster("wgs_v2.report.wgs_v21")
    print res1
    a = Sendtowfm({"id": "tsg_0001", "msg": "stop"})
    c = a.run_cmd()
    print c
