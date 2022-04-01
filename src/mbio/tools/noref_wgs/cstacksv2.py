# -*- coding: utf-8 -*-
# __author__ = 'hongdongxuan'
# last modify 20190104

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import datetime
import random
import json


class Cstacksv2Agent(Agent):
    """
    cstacks  -s ./yanxiaoling_MJ20181017046/04.noref/03.ustacks/CG_07 -n 4 -p 16 
     -o ./0/catalog 2>./04.cstacks-new/0/cstacks.US01_1.log
    """
    def __init__(self, parent):
        super(Cstacksv2Agent, self).__init__(parent)
        options = [
            {"name": "ustacksfiles", "type": "string"},
            {"name": 'outnum', 'type': 'string'}
        ]
        self.add_option(options)
        self.step.add_steps('Cstacks')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Cstacks.start()
        self.step.update()

    def step_end(self):
        self.step.Cstacks.finish()
        self.step.update()

    def check_options(self):
        if not self.option('ustacksfiles'):
            raise OptionError('必须输入:ustacksfiles')
        if not self.option('outnum'):
            raise OptionError('必须输入输出文件夹名字')

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 16
        self._memory = "130G"

    def end(self):
        super(Cstacksv2Agent, self).end()


class Cstacksv2Tool(Tool):
    """
    cstacks  -s /mnt/ilustre/centos7users/dongmei.fu/dna_project/genetic_
    evolution/yanxiaoling_MJ20181017046/04.noref/03.ustacks/CG_07 -n 4 -p 16  -o
     /mnt/ilustre/centos7users/dna/PYRAD/STACKS/04.cstacks-new/0 2>/mnt/ilustre/centos7
     users/dna/PYRAD/STACKS/04.cstacks-new/0/cstacks.CG_07.log
    """
    def __init__(self, config):
        super(Cstacksv2Tool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.cstacks_path = self.config.SOFTWARE_DIR + "/bioinfo/noRefWGS/cstacks"

    def run_cstacks(self):
        """
        {} -s {} -n 4 -p 16 --max_gaps 2 -o {} 2>{} &&
        modified by hd 上述参数中增加了--max_gaps 2
        :return:
        """
        outdir = os.path.join(self.output_dir, self.option('outnum'))
        if os.path.exists(outdir):
            os.system('rm -rf {}'.format(outdir))
            os.mkdir(outdir)
        else:
            os.mkdir(outdir)
        cmd = ""
        i = 0
        for m in json.loads(self.option('ustacksfiles')):
            sample = os.path.basename(m)
            log_path = os.path.join(self.output_dir, ("cstacks." + sample + ".log"))
            catalog_path = os.path.join(outdir, "catalog")
            if i == 0:
                cmd += "{} -s {} -n 4 -p 16 -o {} 2>{} && ".format(self.cstacks_path, m, outdir, log_path)
            else:
                cmd += "{} -s {} -n 4 -p 16 -o {} -c {} 2>{} && "\
                    .format(self.cstacks_path, m, outdir, catalog_path, log_path)
            i += 1
        cmd1 = cmd.strip("&& ")
        self.logger.info("cmd:{}".format(cmd1))
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "cstacksv2_{}.sh".format(now_time)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd1)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始进行run_cstacksv2")
        command1 = self.add_command("run_cstacksv2", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("run_cstacksv2完成！")
        else:
            self.set_error("run_cstacksv2出错！")
        os.system('rm {}'.format(file_path))

    def run(self):
        super(Cstacksv2Tool, self).run()
        self.run_cstacks()
        self.end()
