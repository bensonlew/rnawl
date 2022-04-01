# -*- coding: utf-8 -*-
# __author__ = 'hongdongxuan'
# last modify 20200224

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import datetime
import random


class CstacksMergeAgent(Agent):
    """
    输入文件为/mnt/ilustre/centos7users/dna/PYRAD/STACKS/04.cstacks-new/0/catalog
    """
    def __init__(self, parent):
        super(CstacksMergeAgent, self).__init__(parent)
        options = [
            {"name": "catalog_path", "type": "string"}  # 该参数是cstacks分割后生成的文件所在的文件夹
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
        if not self.option('catalog_path'):
            raise OptionError('必须输入:catalog_path', code="35500203")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 16
        self._memory = "200G"

    def end(self):
        super(CstacksMergeAgent, self).end()


class CstacksMergeTool(Tool):
    def __init__(self, config):
        super(CstacksMergeTool, self).__init__(config)
        self.file_names = []
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.cstacks_path = self.config.SOFTWARE_DIR + "/bioinfo/noRefWGS/cstacks"

    def make_file_list(self):
        """
        重组文件列表
        """
        for m in os.listdir(self.option("catalog_path")):
            if os.path.isdir(os.path.join(self.option("catalog_path"), m)):
                try:
                    if int(m) not in self.file_names:
                        self.file_names.append(int(m))
                except Exception:
                    continue
                else:
                    pass
            else:
                continue
        self.file_names.sort()
        self.logger.info(self.file_names)
        pass

    def run_script(self):
        """
        cstacks  -s /mnt/ilustre/centos7users/dna/PYRAD/STACKS/04.cstacks-new/0/catalog 
        -n 4 -p 16  -o /mnt/ilustre/centos7users/dna/PYRAD/STACKS/04.cstacks-new/ 2>
        /mnt/ilustre/centos7users/dna/PYRAD/STACKS/04.cstacks-new/cstacks.0.log
        :return:
        """
        os.system('rm -rf {}/*'.format(self.output_dir))
        cmd1 = ''
        i = 0
        for n in self.file_names:
            s_params = os.path.join(self.option("catalog_path"), str(n), "catalog")
            l_params = os.path.join(self.output_dir, "cstacks." + str(n) + '.log')
            if i == 0:
                cmd1 += "{} -s {} -n 4 -p 16 -o {} 2>{} && ".format(self.cstacks_path, s_params,
                                                                    self.output_dir, l_params)
            else:
                cmd1 += "{} -s {} -n 4 -p 16 -o {} -c {} 2>{} && ".format(self.cstacks_path, s_params,
                                                                          self.output_dir,
                                                                          os.path.join(self.output_dir, 'catalog'),
                                                                          l_params)
            i += 1
        cmd1 = cmd1.strip("&& ")
        self.logger.info("cmd:{}".format(cmd1))
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "cstacks_{}.sh".format(now_time)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd1)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始进行run_cstacks_merge")
        command1 = self.add_command("run_cstacks_merge", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("run_cstacks_merge完成！")
        else:
            self.set_error("run_cstacks_merge出错！")
        os.system('rm {}'.format(file_path))

    def run(self):
        super(CstacksMergeTool, self).run()
        self.make_file_list()
        self.run_script()
        self.end()
