# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import psutil
import time
import threading

__author__ = 'liubinxu'


class Trinity2DistributeAgent(Agent):
    """
    Trinity
    """
    def __init__(self, parent):
        super(Trinity2DistributeAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'distribute_cmd', 'format': 'denovo_rna_v2.common'},
            {'default': 2, 'type': 'int', 'name': 'cpu'}, # cpu不要太大否则Java容易报错
            {'default': 0, 'type': 'float', 'name': 'mis_rate'},
            {"name": "trinity_version", "type": "string", "default": "2.8.5"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("distribute_cmd").is_set:
            raise OptionError("必须设置参数distribute_cmd", code = "32006201")
        if self.option("cpu")>5:
            raise OptionError("cpu设置过大，不要大于 5", code = "32006202")

    def set_resource(self):
        self._cpu = self.option('cpu')
        self._memory = "{}G".format(self.option('cpu')  * 20)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        super(Trinity2DistributeAgent, self).end()

def parafly_process_control(tool_pid, para_cmd):
    """
    通过tool的process id 监控子线程中的butterfly
    """
    with open(para_cmd, 'r') as para_f:
        para_num = len(para_f.readlines())
    no_bfly = 0
    time_start = time.time()


    with open("para_control.txt", 'a+') as control_log:
        while 1:
            c_thread_pid = os.getpid()
            time_step = 10
            p = psutil.Process(tool_pid)
            p_children = p.children(recursive=True)
            bfly_process_num = 0
            # 监控子进程任务数, butterfly有无超时
            is_run_para = False
            for p_child in p_children:
                try:
                    cmd_name = str(p_child.cmdline())
                except:
                    control_log.write("进程不存在 {}\n".format(p_child.pid))
                    cmd_name = ""
                try:
                    cmd_time = float(p_child.cpu_times().user)
                except:
                    cmd_time = 0

                if 'ParaFly' in cmd_name:
                    is_run_para = True
                if "Butterfly.jar" in cmd_name and cmd_time > 300:
                    control_log.write("进程{} {}使用超时{}s，自动结束\n".format(str(p_child.pid), cmd_name, cmd_time))
                    p_child.terminate()
                else:
                    control_log.write("进程{} {} 运行时间{}s\n".format(str(p_child.pid), cmd_name, cmd_time))
                    if "Butterfly.jar" in cmd_name:
                        bfly_process_num += 1
            # 计算parafly完成数
            finish_num = 0
            if os.path.exists(para_cmd + ".completed"):
                with open(para_cmd + ".completed", 'r') as para_finish_f:
                    finish_num = len(para_finish_f.readlines())

            if bfly_process_num == 0:
                no_bfly += time_step
            else:
                no_bfly = 0

            time_now = time.time()
            time_str = time.localtime()

            # control_log.write("time is {} thread id is {}\n".format(time_str, c_thread_pid))
            '''
            if bfly_process_num == 0 and float(finish_num)/float(para_num) > 0.8:
                control_log.write("程序完成,监控结束\n")
                break
            elif time_now - time_start > 3000 and no_bfly > 1000:
                control_log.write("程序已运行超过3000秒 且 连续1000秒没有butterfly脚本运行,监控结束\n")
                break
            '''
            if not is_run_para:
                control_log.write("程序完成,监控结束\n")
                break
            elif time_now - time_start > 3000 and no_bfly > 1000:
                control_log.write("程序已运行超过3000秒 且 连续1000秒没有butterfly脚本运行,监控结束\n")
                break
            else:
                time.sleep(time_step)
        control_log.write("循环结束")
    return

class Trinity2DistributeTool(Tool):
    """
    Trinity
    """
    def __init__(self, config):
        super(Trinity2DistributeTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        if self.option("trinity_version") == "2.5.0":
            self.parafly = '/bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/trinity-plugins/ParaFly-0.1.0/bin/ParaFly'
        elif self.option("trinity_version") == "2.8.5":
            self.parafly = '/bioinfo/denovo_rna_v2/trinityrnaseq-2.8.5/trinity-plugins/ParaFly-0.1.0/bin/ParaFly'

        else:
            self.set_error("trinity version error", code="32006205")

        self.perl = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/bin/'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/denovo_rna_v2/bowtie2-2.3.3.1-linux-x86_64/'
        self.java = self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin'
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.set_environ(PATH=self.perl)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.java)

        if self.option("trinity_version") == "2.8.5":
            self.salmon = self.config.SOFTWARE_DIR + '/bioinfo/rna/salmon-0.14.1/bin/'
            self.jellyfish = self.config.SOFTWARE_DIR + '/bioinfo/denovo_rna_v2/jellyfish-2.3.0/bin/'
            self.set_environ(PATH=self.salmon)
            self.set_environ(PATH=self.jellyfish)
            self.gcc = software_dir + '/gcc/5.1.0/bin'
            self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
            self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)



    def run_parafly(self):
        cmd = '{} '.format(self.parafly)
        cmd += '-{} {} '.format("c", self.option("distribute_cmd").prop['path'])
        cmd += '-{} {} '.format("CPU", self.option("cpu"))
        cmd += '-v -shuffle'
        cmd_name = 'trinity2distribute'
        command = self.add_command(cmd_name, cmd, ignore_error=True)

        # 启动监控
        tool_pid = os.getpid()
        command.run()
        try:
            control_t = threading.Thread(target=parafly_process_control, args=(tool_pid, self.option("distribute_cmd").prop['path']))
            control_t.start()
            self.logger.info("启动监控线程")
        except:
            self.set_error("无法启动监控线程", code="32006206")

        self.wait()
        control_t.join()

        self.logger.info("启动监控线程结束")
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif os.path.exists(os.path.join(self.work_dir, 'FailedCommands')):
            fail_file = open(os.path.join(self.work_dir, 'FailedCommands'))
            fail_num = len(fail_file.readlines())
            if fail_num < 50:
                self.logger.info("unfinished command little than 50 will jumped")
            else:
                self.logger.info("unfinished command more than 50 will rerun")
                command.rerun()
                self.wait()
                # self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32006203")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32006204")

    def set_output(self):
        pass

    def run(self):
        super(Trinity2DistributeTool, self).run()
        self.run_parafly()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20180102/Single_de_tr_data2.2'
        data = {
            "id": "Trinity2Distribute" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.trinity2_distribute",
            "instant": False,
            "options": dict(
                distribute_cmd=test_dir + "/" + "trinity_cmd_distribute09",
                cpu="4",
                mis_rate="0",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
