# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modify 20180515


from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class CircosChrlistAgent(Agent):
    """
    参考基因组:
    ref.chrlist:全部chrID和length；或者排名长度从大到小前20个的scaID和length；
    total.chrlist:全部chrID\tlength
    """
    def __init__(self, parent):
        super(CircosChrlistAgent, self).__init__(parent)
        options = [
            {"name": "dict", "type": "string"}  # one col
        ]
        self.add_option(options)
        # self.step.add_steps('CircosChrlist')
    #     self.on('start', self.step_start)
    #     self.on('end', self.step_end)

    # def step_start(self):
    #     self.step.CircosChrlist.start()
    #     self.step.update()

    # def step_end(self):
    #     self.step.CircosChrlist.finish()
    #     self.step.update()

    def check_options(self):
        if not self.option("dict"):
            raise OptionError("必须输入ref.dict文件", code="34501201")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(CircosChrlistAgent, self).end()


class CircosChrlistTool(Tool):
    def __init__(self, config):
        super(CircosChrlistTool, self).__init__(config)
        self.chr_path = self.config.PACKAGE_DIR + "/wgs/chr.pl"  
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '

    def totallist(self):
        """
        dict生成total.chrlist
        """
        os.system("less {}|grep -v 'HD'|sed 's/:/\t/g'|cut -f 3,5 >{}/total.chrlist".format(self.option("dict"), self.output_dir))

    def circoschrlist(self):
        """
        :return:
        """
        cmd = "{} {} -i {} -o {}" \
            .format(self.perl_path, self.chr_path, self.option("dict"),
                    self.output_dir+"/ref.chrlist")
        self.logger.info(cmd)
        self.logger.info("开始进行CircosChrlist")
        command = self.add_command("circoschrlist", cmd).run()  # circoschrlist必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("CircosChrlist完成！")
        else:
            self.set_error("CircosChrlist出错！", code="34501201")
            self.set_error("CircosChrlist出错！", code="34501204")

    def run(self):
        super(CircosChrlistTool, self).run()
        self.circoschrlist()
        self.totallist()
        self.end()
