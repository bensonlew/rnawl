# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180830

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime
import random


class BayesTreeAgent(Agent):
    """
    群体结构子模块，进化树中的贝叶斯的算法计算树，--该脚本运行会有很长时间，甚至一个月两个月
    """
    def __init__(self, parent):
        super(BayesTreeAgent, self).__init__(parent)
        options = [
            {"name": "pop_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "model_test_out", "type": "infile", "format": "dna_evolution.model_test"}
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
            raise OptionError("缺少%s，请添加%s!", variables=("pop_fasta", "pop_fasta_file"), code="11111")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 17
        self._memory = '100G'

    def end(self):
        super(BayesTreeAgent, self).end()


class BayesTreeTool(Tool):
    def __init__(self, config):
        super(BayesTreeTool, self).__init__(config)
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.script_path = self.config.PACKAGE_DIR + "/dna_evolution/bayes-prepair.pl"
        self.script_path1 = self.config.SOFTWARE_DIR + "/bioinfo/dna_evolution/mpirun"
        self.script_path2 = self.config.SOFTWARE_DIR + "/bioinfo/dna_evolution/MrBayes-3.1.2h/mb"

    def bayesprepair(self):
        """
        perl /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/dna_evolution/bayes-prepair.pl
        -i /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/jinhua/tree/pop.fasta
        -m /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/jinhua/tree/new/pop.model.test.out
         -o /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/jinhua/tree/bayes/bayes
        :return:
        """
        cmd1 = "{} {} -i {} -m {} -o {}" \
            .format(self.perl_path, self.script_path, self.option("pop_fasta").prop['path'],
                    self.option("model_test_out").prop['path'], self.output_dir + "/bayes")
        self.logger.info(cmd1)
        self.logger.info("开始进行bayes-prepair")
        command = self.add_command("bayes_prepair", cmd1).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bayes-prepair完成！")
        else:
            self.set_error("bayes-prepair出错！")

    def bayes_run(self):
        """
        命令2：
        mpirun -np 16 mb /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/jinhua/tree/bayes/bayes.run.nex
         > /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/jinhua/tree/bayes/bayes.log
        :return:
        """
        cmd2 = "cd {} && {} -np 16 {} {} > {}".format(self.output_dir, self.script_path1, self.script_path2,
                                                      self.output_dir + "/bayes.run.nex", self.work_dir + "/bayes.log")
        self.logger.info(cmd2)
        self.sbatch_task(cmd2, "bayes_run")

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
        super(BayesTreeTool, self).run()
        self.bayesprepair()
        self.bayes_run()
        self.end()
