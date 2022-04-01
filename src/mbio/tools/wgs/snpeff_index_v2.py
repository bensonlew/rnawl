# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modify 20180515
# modify 20180523

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SnpeffIndexV2Agent(Agent):
    """
    参考基因组snpeff的注释;
    结果产生在ref.fa所在目录了。
    """
    def __init__(self, parent):
        super(SnpeffIndexV2Agent, self).__init__(parent)
        options = [
            {"name": "reffa", "type": "string"},    # 链接到指定目录的ref里变成sequences.fa
            {"name": "genesgtf", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("reffa"):
            raise OptionError("请设置reffa", code="34505901")
        if not self.option("genesgtf"):
            raise OptionError("请设置genesgtf", code="34505902")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 10
        self._memory = '20G'

    def end(self):
        super(SnpeffIndexV2Agent, self).end()


class SnpeffIndexV2Tool(Tool):
    def __init__(self, config):
        super(SnpeffIndexV2Tool, self).__init__(config)
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin/') # 加载java
        self.java_path = "program/sun_jdk1.8.0/bin/java"
        self.snpeff_path = self.config.SOFTWARE_DIR + '/bioinfo/annotation/snpeff/'
        self.reffa_path, self.reffa_name = os.path.split(self.option("reffa"))

    def refmkdir(self):
        self.reffa_path = self.reffa_path.strip()   # 去除首位空格
        path = self.reffa_path + "/ref"
        isExists = os.path.exists(path)
        if not isExists:       # 判断路径是否存在;存在     True;不存在   False
            os.makedirs(path)
            print path + ' 创建成功'
            return True
        else:
            print path + ' 目录已存在'       # 如果目录存在则不创建，并提示目录已存在
            os.system("rm -fr {}".format(path))
            os.makedirs(path)
            print path + ' 目录已存在,重新建立成功'
            return True

    def snpeff_config(self):
        if os.path.exists(self.reffa_path + "/ref"):    # 检查
            os.link(self.option("reffa"), self.reffa_path + "/ref/sequences.fa")
            os.link(self.option("genesgtf"), self.reffa_path + "/ref/genes.gtf")  # 链接前的ref目录已删除
            filename = self.reffa_path + "/snpEff.config"
            with open(filename, 'w') as f:       # 如果filename不存在会自动创建， 'w'表示写数据，
                # 写之前会清空文件中的原有数据！
                f.write("data.dir ={}/\n".format(self.reffa_path))
                f.write("ref.genome : ref\n")

    def SnpeffIndexV2(self):
        """
        java -jar snpEff.jar build -gff3 -v ref -c snpEff.config
        """
        cmd = "{} -jar {} build -gtf22 -v ref -c {}"\
            .format(self.java_path, self.snpeff_path + "/snpEff.jar", self.reffa_path + "/snpEff.config")
        self.logger.info(cmd)
        self.logger.info("开始进行SnpeffIndexV2")
        command = self.add_command("snpeffbin", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("SnpeffIndexV2完成！")
        else:
            self.set_error("SnpeffIndexV2出错！", code="34505901")

    def run(self):
        super(SnpeffIndexV2Tool, self).run()
        self.refmkdir()
        self.snpeff_config()
        self.SnpeffIndexV2()
        self.end()
