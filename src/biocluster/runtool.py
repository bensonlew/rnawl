#!/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import argparse
import pickle
# from biocluster.agent import PickleConfig
import os
from biocluster.core.function import load_class_by_path, daemonize
import traceback
import re
import time
import platform
from biocluster.core.exceptions import RunningError, CodeError
import json
import sys
import grpc
import socket
from biocluster.proto import tool_guide_pb2_grpc, tool_guide_pb2

os.environ["LANG"] ="en_US.UTF-8"

parser = argparse.ArgumentParser(description="run a tool on remote server")
parser.add_argument("-b", "--daemon", action="store_true", help="run in daemon background mode")
parser.add_argument("-d", "--debug", action="store_true", help="run in debug mode,will not use network!")
parser.add_argument("-p", "--port", type=int, help="the ntm service port")
parser.add_argument("-c", "--cpu", type=int, help="request cpu for tool")
parser.add_argument("-m", "--mem", type=int, help="request memory for tool")
parser.add_argument("tool", type=str, help="the tool name to run")
args = parser.parse_args()

os.environ['current_mode'] = "tool"
if args.port:
    os.environ["NTM_PORT"] = str(args.port)
request_cpu = 0
request_memory = 0
if args.cpu:
    request_cpu = args.cpu
if args.mem:
    request_memory = args.mem

def main():
    hostname = socket.gethostname()
    print (hostname)
    ntm_file = "ntm_no"
    if os.path.exists(ntm_file):
        os.remove(ntm_file)
    name = args.tool
    class_file = name + "_class.pk"
    if args.daemon:
        daemonize(stdout="%s.daemon.o" % name, stderr="%s.daemon.e" % name)
        write_pid()
    # with open(class_file, "r") as f:
    with open(class_file, "rb") as f:
        class_list = pickle.load(f)
    for file_class in class_list['files']:
        load_class_by_path(file_class, "File")
    # print class_list['tool']
    paths = class_list['tool'].split(".")
    paths.pop(0)
    paths.pop(0)
    tool = load_class_by_path(".".join(paths), "Tool")
    config_file = name + ".pk"
    # with open(config_file, "r") as f:
    with open(config_file, "rb") as f:
        config = pickle.load(f)
    # print vars(config)
    if args.debug:
        config.DEBUG = True
    else:
        config.DEBUG = False
    if args.port:
        config.ntm_port = args.port
    # else:
        # config.ntm_port = 7322
    config.current_mode = "tool"
    config.request_cpu = request_cpu
    config.request_memory = request_memory
    ss = tool_send_state(config, "runtool", "")
    if ss:
        with open("ntm_no", "w") as f:
            f.write("%s" % platform.uname()[1])
        print ("发送state %s超过3次仍然失败,ntm服务grpc不正常,退出运行" % "runtool")
        raise ss
    t = None
    try:
        t = tool(config)
        t.run()
    except MemoryError as e:
        exstr = traceback.format_exc()
        sys.stderr.write(exstr)
        sys.stdout.write(exstr)
        sys.stderr.flush()
        t.add_state("memory_limit", "MemoryError:检测到内存使用时出现错误!")
    except Exception as e:
        exstr = traceback.format_exc()
        sys.stderr.write(exstr)
        sys.stdout.write(exstr)
        sys.stderr.flush()
        finded = False
        if t:
            t._end = True
            # t.save_report()
            if "SLURM_JOB_ID" in os.environ.keys():
                time.sleep(1)
                t.logger.debug("开始检测SLURM STDERR输出...")
                slurm_error_path = os.path.join(t.work_dir, "%s_%s.err" % (t.name, os.environ["SLURM_JOB_ID"]))
                error_msg = ""
                exceeded = False
                cancelled = False
                with open(slurm_error_path, "r") as f:
                    f.seek(0, 2)
                    size = os.path.getsize(slurm_error_path)
                    point = 5000 if size > 5000 else size
                    f.seek(-point, 2)
                    lines = f.readlines()
                    for line in lines:
                        if re.match(r"^slurmstepd:", line):
                            error_msg += line
                            if re.search(r"memory limit", line):
                                exceeded = True
                            elif re.search(r"CANCELLED", line):
                                cancelled = True
                if exceeded:
                    t.logger.info("检测到内存使用超过申请数被系统杀死!")
                    t.add_state("memory_limit", error_msg)
                    finded = True
                elif cancelled:
                    t.logger.info("检测到任务被取消!")
                    t.add_state("cancelled", error_msg)
                    finded = True
            if not finded:
                t.set_error(e)
        else:
            print "workflow_id:{}, tool-id:{}".format(config.current_workflow_id, config.current_tool_id)
            ss = tool_send_state(config, "error", e)
            if ss:
                raise ss
    sys.exit(0)


send_state_times = 0


def send_state(state_data, config):
    global send_state_times
    send_state_times += 1

    def get_state_data(state):
        yield tool_guide_pb2.State(
            workflow_id=state.workflow_id,
            tool_id=state.tool_id,
            jobid=state.jobid,
            jobtype=state.jobtype,
            state=state.state,
            process_id=state.process_id,
            host=state.host,
            version=state.version,
            data=state.data)
    try:
        with grpc.insecure_channel('localhost:%s' % config.ntm_port) as channel:
            stub = tool_guide_pb2_grpc.ToolGuideStub(channel)
            response = stub.SendState(get_state_data(state_data))
            if response.ok:
                print ("发送state %s成功" % state_data.state)
            else:
                print ("ntm服务拒绝接受state %s, 原因: %s" % (state_data.state, response.reason))
            return ""
    except Exception as e:
        exstr = traceback.format_exc()
        print(exstr)
        sys.stdout.flush()
        if send_state_times > 5:
            print ("发送state %s超过3次仍然失败,退出运行" % state_data.state)
            # raise e
            return e
        else:
            print ("发送state %s失败,30秒后重新尝试" % state_data.state)
            time.sleep(30)
            return send_state(state_data, config)

def tool_send_state(config, state, error):
    """
    给ntm发送状态，检测ntm服务grpc是否正常
    """
    jobid = 0
    if "SLURM_JOB_ID" in os.environ.keys():
        try:
            jobid = int(os.environ["SLURM_JOB_ID"])
        except:
            pass
    if "PBS_JOBID" in os.environ.keys():
        try:
            jobid = int(os.environ["PBS_JOBID"])
        except:
            pass
    if isinstance(error, CodeError):
        error_msg = error.json()
        error_msg["name"] = getattr(config, "_id")
    else:
        error_msg = {"error_type": "running",
                     "name": getattr(config, "_id"),
                     "code": "R001",
                     "variables": None,
                     "info": str(error)
                     }
    state_data = tool_guide_pb2.State(
        workflow_id=config.current_workflow_id,
        tool_id=config.current_tool_id,
        jobid=jobid,
        jobtype="SLURM",
        state=state,
        process_id=int(os.getpid()),
        host=platform.uname()[1],
        version=0,
        data=json.dumps(error_msg),
    )
    print "workflow_id:{}, tool-id:{}".format(config.current_workflow_id, config.current_tool_id)
    ss = send_state(state_data, config)
    return ss

def write_pid():
    with open("run.pid", "w") as f:
        f.write("%s" % os.getpid())


if __name__ == '__main__':
    main()
