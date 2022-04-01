# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180428

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os 


class EmapperAgent(Agent):
    """
    EGGNOG数据库标准注释流程
    """
    def __init__(self,parent):
        super(EmapperAgent,self).__init__(parent)
        options = [
            {"name": "dmnd_db", "type": "string"},  # db file
            {"name": "db_path", "type": "string"},  # emapper的path
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},  # split's fasta
        ]
        self.add_option(options)
        self.step.add_steps('EggNOGmapper')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.EggNOGmapper.start()
        self.step.update()

    def step_end(self):
        self.step.EggNOGmapper.finish()
        self.step.update()

    def check_options(self):
        if not self.option("dmnd_db"):
            raise OptionError("请设置dmnd_db", code="34502401")
        if not self.option("db_path"):
            raise OptionError("db_path EggNog.fa等所在路径", code="34502402")
        if not self.option("fasta"):
            raise OptionError("请输入待比对fasta", code="34502403")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 13
        self._memory = '30G'

    def end(self):
        super(EmapperAgent, self).end()


class EmapperTool(Tool):
    def __init__(self, config):
        super(EmapperTool, self).__init__(config)   # diamond加入环境变量
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/align/diamond-0.8.35/")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/bioinfo/align/diamond-0.8.35/")
        self.emapper_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/eggnog-mapper/emapper.py"
        self.python_path = "/program/Python/bin/python"

    def Emapper(self):
        """
        要重新写下！！！
        :return:
        """
        name = os.path.basename(self.option("fasta").prop["path"])+".blast"    # 需要加进去数据库的type吗?如nr.blast
        cmd = "{} {} --translate --dmnd_db {} --dbtype seqdb --override  -m diamond  --data_dir  {} -i {} -o {} " \
              "--cpu 12".format(self.python_path, self.emapper_path, self.option("dmnd_db"),
                                self.option("db_path"), self.option("fasta").prop["path"], self.output_dir + "/" + name)
        self.logger.info(cmd)
        self.logger.info("开始进行Emapper")
        command = self.add_command("emapper", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Emapper完成！")
        else:
            self.set_error("Emapper出错！", code="34502401")
            self.set_error("Emapper出错！", code="34502404")

    def run(self):
        super(EmapperTool, self).run()
        self.Emapper()
        self.end()
