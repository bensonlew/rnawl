#!/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import argparse
import pickle
from biocluster.agent import PickleConfig
import os
from biocluster.core.function import load_class_by_path, daemonize
import traceback
import re
import time
import platform
from biocluster.core.exceptions import RunningError, CodeError
import json
import zerorpc
import sys

parser = argparse.ArgumentParser(description="run a tool on remote server")
parser.add_argument("-b", "--daemon", action="store_true", help="run in daemon background mode")
parser.add_argument("-d", "--debug", action="store_true", help="run in debug mode,will not use network!")
parser.add_argument("tool", type=str,
                    help="the tool name to run")
args = parser.parse_args()


def main():
    name = args.tool
    class_file = name + "_class.pk"
    if args.daemon:
        daemonize(stdout="%s.daemon.o" % name, stderr="%s.daemon.e" % name)
        write_pid()
    with open(class_file, "r") as f:
        class_list = pickle.load(f)
    for file_class in class_list['files']:
        load_class_by_path(file_class, "File")
    # print class_list['tool']
    paths = class_list['tool'].split(".")
    paths.pop(0)
    paths.pop(0)
    tool = load_class_by_path(".".join(paths), "Tool")
    config_file = name + ".pk"
    with open(config_file, "r") as f:
        config = pickle.load(f)
    # print vars(config)
    if args.debug:
        config.DEBUG = True
    else:
        config.DEBUG = False
    t = None
    try:
        t = tool(config)
        t.run()
    except MemoryError, e:
        exstr = traceback.format_exc()
        print >> sys.stderr, exstr
        print >> sys.stderr, exstr
        sys.stderr.flush()
        t.add_state("memory_limit", "检测到内存使用时出现错误!")
    except Exception, e:
        exstr = traceback.format_exc()
        print >>sys.stderr, exstr
        print >> sys.stderr, exstr
        sys.stderr.flush()
        finded = False
        if t:
            t._end = True
            t.save_report()
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
            jobid = None
            if "SLURM_JOB_ID" in os.environ.keys():
                jobid = os.environ["SLURM_JOB_ID"]
            if "PBS_JOBID" in os.environ.keys():
                jobid = os.environ["PBS_JOBID"]
            if isinstance(e, CodeError):
                error_msg = e.json()
                error_msg["name"] = getattr(config, "_id")
            else:
                error_msg = {"error_type": "running",
                             "name": getattr(config, "_id"),
                             "code": "R001",
                             "variables": None,
                             "info": str(e)
                             }
            data = {"id": getattr(config, "_id"),
                    "state": "error",
                    "data": json.dumps(error_msg),
                    "jobid": jobid,
                    "host": platform.uname()[1],
                    "version": 0
                    }
            client = zerorpc.Client()
            client.connect(config.endpoint)
            client.report(data)
            client.close()
    sys.exit(0)


def write_pid():
    with open("run.pid", "w") as f:
        f.write("%s" % os.getpid())

if __name__ == '__main__':
    main()
