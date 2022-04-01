# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 2018080907

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class GctaPcaAgent(Agent):
    """
    群体结构子模块，计算pca -- # 输入的是plink--make-bed的目录
    """
    def __init__(self, parent):
        super(GctaPcaAgent, self).__init__(parent)
        options = [
            {"name": "bfile_dir", "type": "infile", "format": "dna_evolution.gcta_dir", "required": True}
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
        pass

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 3
        self._memory = '50G'

    def end(self):
        super(GctaPcaAgent, self).end()


class GctaPcaTool(Tool):
    def __init__(self, config):
        super(GctaPcaTool, self).__init__(config)
        self.script_path = "/bioinfo/dna_evolution/gcta"

    def gcta_bfile(self):
        """
        gcta --bfile /mnt/ilustre/users/minghao.zhang/newmdt/Project/lidonghua/pop_20180621/step03.pca-generic/pop
        --make-grm --out /mnt/ilustre/users/minghao.zhang/newmdt/Project/lidonghua/pop_20180621/step03.pca-generic/pop
        生成pop.grm.id,pop.grm.N.bin,pop.grm.bin
        :return:
        """
        cmd = "{} --bfile {}/pop --make-grm --out {}".format(self.script_path, self.option("bfile_dir").prop['path'],
                                                            self.output_dir + "/pop")
        self.run_cmd(cmd, "gcta_bfile")

    def gcta_grm(self):
        """
        gcta --grm /mnt/ilustre/users/minghao.zhang/newmdt/Project/lidonghua/pop_20180621/step03.pca-generic/pop
        --pca 20 --out /mnt/ilustre/users/minghao.zhang/newmdt/Project/lidonghua/pop_20180621/step03.pca-generic/pop.pca
        生成pop.pca.eigenval，pop.pca.eigenvec
        :return:
        """
        cmd = "{} --grm {} --pca 20 --out {}".format(self.script_path, self.output_dir + "/pop",
                                                     self.output_dir + "/pop.pca")
        self.run_cmd(cmd, "gcta_grm")

    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("{}运行失败".format(cmd_name))

    def run(self):
        super(GctaPcaTool, self).run()
        self.gcta_bfile()
        self.gcta_grm()
        self.end()
