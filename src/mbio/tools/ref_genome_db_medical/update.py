# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import subprocess
import unittest


class UpdateAgent(Agent):
    """
    用于参考基因组更新上线
    """

    def __init__(self, parent):
        super(UpdateAgent, self).__init__(parent)
        options = [
            {"name": "genome_id", "type": "string", 'default': ""},
            {"name": "result_dir", "type": "string", 'default': ""},
        ]
        self.add_option(options)
        self.step.add_steps("update_genome")
        self.on('start', self.start_update_genome)
        self.on('end', self.end_update_genome)

    def start_update_genome(self):
        self.step.update_genome.start()
        self.step.update()

    def end_update_genome(self):
        self.step.update_genome.finish()
        self.step.update()

    def check_option(self):
        if not self.option("genome_id"):
            raise OptionError("Genome id not exists!")
        if not os.path.exists(self.opiton("result_dir")):
            raise OptionError("Genome file not exists!")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = "2G"


class UpdateTool(Tool):
    def __init__(self, config):
        super(UpdateTool, self).__init__(config)
        self.update_mongo_script = self.config.PACKAGE_DIR + "/ref_genome_db_v2/update_db_to_sanger.py"
        self.update_disk_script = self.config.PACKAGE_DIR + "/ref_genome_db_v2/syn2nb2.sh"

    def update_disk(self):
        cmd = "{} {}".format(self.update_disk_script, self.option("result_dir"))
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("文件上传完成")
        except subprocess.CalledProcessError:
            self.set_error("文件上传失败!!")

    def update_mongo(self):
        cmd = "{} {} --genome_ids {}".format(self.config.SOFTWARE_DIR + '/program/Python/bin/python', self.update_mongo_script, self.option("genome_id"))
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("数据库更新完成")
        except subprocess.CalledProcessError:
            self.set_error("数据库更新失败")

    def run(self):
        super(UpdateTool, self).run()
        self.update_disk()
        self.update_mongo()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Update_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db_v2.update",
            "instant": False,
            "options": dict(
                genome_id="GM0578",
                result_dir="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Sclerotinia_sclerotiorum_test3/GCA_001857865.1_ASM185786v1",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()