# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180830

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import random
import os
import re


class MlTreeAgent(Agent):
    """
    群体结构子模块，ml tree的计算
    """
    def __init__(self, parent):
        super(MlTreeAgent, self).__init__(parent)
        options = [
            {"name": "model_test_out", "type": "infile", "format": "dna_evolution.model_test"},
            {"name": "bs_trees", "type": "int", "default": 1000},
            {"name": "threads", "type": "string", "default": "16"},
            # {"name": "group", "type": "infile", "format": "dna_evolution.group_table"},
            # {"name": "outgroup", "type": "infile", "format": "dna_evolution.outgroup"}
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
        if not self.option("model_test_out").is_set:
            raise OptionError("缺少%s参数，请添加%s!", variables=("model_test_out", "model_test_out参数"), code="11111")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 17
        self._memory = '50G'

    def end(self):
        super(MlTreeAgent, self).end()


class MlTreeTool(Tool):
    def __init__(self, config):
        super(MlTreeTool, self).__init__(config)
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/dna_evolution/"
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/7.2.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/7.2.0/lib64')
        self.phylip_path = ""

    def vcf2tree(self):
        """
        tail -n6 /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/jinhua/tree/new/pop.model.test.out|
        grep raxml-ng|sed 's/>//g'|tr -d '\n'

        /mnt/ilustre/users/sanger-dev/app/bioinfo/dna_evolution/raxml-ng
        --msa /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/evolution/jinhua/tree/pop.phylip --model TVM+G4
        --threads 16 --bs-trees 1000
        :return:
        """
        cmd1 = "tail -n6 {}|grep raxml-ng|sed 's/>//g'|tr -d '\n' > {}"\
            .format(self.option("model_test_out").prop['path'], self.work_dir + "/cmd1.txt")
        self.logger.info(cmd1)
        self.sbatch_task(cmd1, "raxml_cmd1")
        with open(self.work_dir + "/cmd1.txt", "r") as r:
            data = r.readline()
            tmp = data.strip().split(" ")
            self.phylip_path = tmp[2]
            self.logger.info(self.phylip_path)
            cmd2 = "{}{} --threads {} --bs-trees {} --all --force".format(self.script_path, data.strip(),
                                                                          self.option("threads"),
                                                                          self.option("bs_trees"))
        self.logger.info(cmd2)
        self.sbatch_task(cmd2, "raxml_cmd2")

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

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(os.path.dirname(self.phylip_path))
        for f in results:
            if re.match(".*\.raxml\..*", f):
                os.link(os.path.dirname(self.phylip_path) + '/' + f, self.output_dir + '/' + f)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(MlTreeTool, self).run()
        self.vcf2tree()
        self.set_output()
        self.end()
