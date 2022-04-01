# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 2018080902

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import random
import os


class NjTreeAgent(Agent):
    """
    群体结构子模块，计算nj tree
    """
    def __init__(self, parent):
        super(NjTreeAgent, self).__init__(parent)
        options = [
            {"name": "pop_fasta", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)
        self.step.add_steps('script')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.script.start()
        self.step.update()

    def step_end(self):
        self.step.script.finish()
        self.step.update()

    def check_options(self):
        if not self.option("pop_fasta").is_set:
            raise OptionError("缺少%s，请添加%s参数!", variables=("pop_fasta", "fasta sequence"), code="11111")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 17
        self._memory = '50G'

    def end(self):
        super(NjTreeAgent, self).end()


class NjTreeTool(Tool):
    def __init__(self, config):
        super(NjTreeTool, self).__init__(config)
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/dna_evolution/FastTreeMP"

    def vcf2tree(self):
        """
        /mnt/ilustre/users/sanger-dev/app/bioinfo/dna_evolution/FastTreeMP
        /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/jinhua/tree/pop.fasta >
        /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/jinhua/tree/FastTreeMP/pop.nj.tree
        :return:
        """
        os.system("export OMP_NUM_THREADS=16")
        cmd = "{} {} > {}".format(self.script_path, self.option("pop_fasta").prop['path'],
                                  self.output_dir + "/pop.nj.tree")
        self.logger.info(cmd)
        self.sbatch_task(cmd)

    def sbatch_task(self, cmd, cmd_name="script"):
        """
        该模块是为了解决add_command中有> |等定向符号导致的报错，所以这里就写一个外包程序
        :params: cmd 正常的tool中组建的运行shell命令，这里可以包含有定向符，具体执行命令一定要写成绝对路径，或者加载对应的环境变量
        :params: cmd_name 同一个tool中 如果要执行多个cmd要保证cmd_name不唯一
        :return:
        """
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/dna_evolution/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "script_{}.sh".format(now_time)
        self.logger.info("执行脚本路径:{}".format(file_path))
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path))
        shell = "/bioinfo/dna_evolution/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始执行:{}".format(shell))
        command1 = self.add_command(cmd_name, shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("脚本{}执行完成！".format(shell))
        else:
            self.set_error("脚本{}执行出错！".format(shell))
        os.system('rm {}'.format(file_path))

    def run(self):
        super(NjTreeTool, self).run()
        self.vcf2tree()
        self.end()
