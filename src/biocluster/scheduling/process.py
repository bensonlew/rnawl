# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .job import Job
import os
import gevent
import re
from multiprocessing import Manager
import pickle
from biocluster.agent import PickleConfig
import os
from biocluster.core.function import load_class_by_path,  get_classpath_by_object, hostname
# from gipc.gipc import _GProcess as Process
from biocluster.logger import Wlog
import gipc
from biocluster.config import Config
import time
import sys
import traceback
from biocluster.core.exceptions import CodeError, RunningError


class PROCESS(Job):
    """
    用于本地直接运行
    """
    def __init__(self, agent):
        super(PROCESS, self).__init__(agent)
        self.agent = agent
        self.workflow = agent.get_workflow()
        if not hasattr(self.workflow, "process_share_manager"):
            self.workflow.process_share_manager = Manager()
        self.shared_callback_action = self.workflow.process_share_manager.dict()
        agent.shared_callback_action = self.shared_callback_action

        self.process = None

    def submit(self):
        super(PROCESS, self).submit()
        # self.process.start()
        self.process = gipc.start_process(local_process_run, args=(self.agent,
                                                                   self.workflow.rpc_server.process_queue,
                                                                   self.shared_callback_action,), daemon=True)
        self.id = self.process.pid
        if self.id:
            self.agent.fire("runstart", hostname)
        return self.id

    def delete(self):
            if self.process.is_alive():
                self.process.terminate()
                # self.manager.shutdown()
            else:
                self.process.join()
                # self.manager.shutdown()

    def set_end(self):
        super(PROCESS, self).set_end()
        self.delete()


def local_process_run(agent, process_queue, shared_callback_action):
    #     # super(LocalProcess, self).__init__()
    #     self.agent = agent
    #     self._process_pipe_writer = process_pipe_writer
    #     self._shared_callback_action = shared_callback_action
    #
    # def run(self):
        # super(LocalProcess, self).run()
        workflow = agent.get_workflow()
        workflow.rpc_server.close()
        # Watcher().stopall()
        gevent.sleep(0)
        print("agent work_dir {}".format(agent.work_dir))
        print("agent name {}".format(agent.name))
        os.chdir(agent.work_dir)
        log = os.path.join(agent.work_dir, "%s_%s.err" % (agent.name, os.getpid()))
        if os.path.exists(log):
            pass
        else:
            os.system("touch %s" % log)
        print("log: %s" % log)
        so = file(log, 'a+')
        se = file(log, 'a+', 0)
        os.dup2(so.fileno(), sys.stdout.fileno())
        os.dup2(se.fileno(), sys.stderr.fileno())
        file_class_paths = []  #
        for option in agent.get_option_object().values():
            if option.type in {'outfile', 'infile'}:
                if option.format:
                    file_class_paths.append(option.format)
                else:
                    for f in option.format_list:
                        file_class_paths.append(f)
        for file_class in file_class_paths:
            load_class_by_path(file_class, "File")

        tool_path = get_classpath_by_object(agent)
        paths = tool_path.split(".")
        paths.pop(0)
        paths.pop(0)
        tool = load_class_by_path(".".join(paths), "Tool")

        config = PickleConfig()
        config.clone(agent)
        config.DEBUG = False
        config.instant = True
        itool = None
        try:
            itool = tool(config)
            itool.logger = Wlog(itool).get_logger('Tool子进程 %s (parent: %s )' % (os.getpid(), os.getppid()))
            itool.shared_callback_action = shared_callback_action
            itool.process_queue = process_queue
            itool.run()
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stderr.flush()
            if itool:
                if isinstance(e, CodeError):
                    e.bind_object = itool
                itool.save_report()
                itool.set_error(e)
            else:
                e1 = RunningError(str(e))
                e1.bind_object = agent
                data = {"id": agent.id,
                        "state": "error",
                        "data": e1.json(),
                        "version": 0
                        }
                process_queue.put(data)

