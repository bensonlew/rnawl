# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180830

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class ModeltestAgent(Agent):
    """
    群体结构子模块，使用modeltest生成后续分析需要的配置文件
    """
    def __init__(self, parent):
        super(ModeltestAgent, self).__init__(parent)
        options = [
            {"name": "phylip", "type": "infile", "format": "dna_evolution.phylip"},
            {"name": "model_test", "type": "outfile", "format": "dna_evolution.model_test"}
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
        if not self.option("phylip").is_set:
            raise OptionError("缺少%s, 请添加参数%s!", variables=("phylip", "phylip_file"), code="11111")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '15G'

    def end(self):
        super(ModeltestAgent, self).end()


class ModeltestTool(Tool):
    def __init__(self, config):
        super(ModeltestTool, self).__init__(config)
        self.script_path = "bioinfo/dna_evolution/modeltest-ng"

    def model_test(self):
        """
        modeltest-ng  -p 4 -d nt -i pop.phylip  -o pop.model.test
        :return:
        """
        cmd = "{} -p 4 -d nt -i {} -o {}"\
            .format(self.script_path, self.option("phylip").prop['path'],
                    self.output_dir + "/pop.model.test")
        self.logger.info(cmd)
        self.logger.info("开始进行script")
        command = self.add_command("modeltest", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("script完成！")
            if os.path.exists(self.output_dir + "/pop.model.test.out"):
                self.option("model_test", self.output_dir + "/pop.model.test.out")
        else:
            self.set_error("script出错！")

    def run(self):
        super(ModeltestTool, self).run()
        self.model_test()
        self.end()
