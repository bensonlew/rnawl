# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest


class BlastIndexAgent(Agent):
    """
    构建blast索引文件
    """
    def __init__(self, parent):
        super(BlastIndexAgent, self).__init__(parent)
        options = [
            {"name": "reference", "type": "infile", "format": "small_rna.fasta"}, # 数据库参考文件
            {"name": "reference_type", "type": "string", "default": "nucl"}, # 数据库参考文件格式
        ]
        self.add_option(options)
        self.step.add_steps("blast_index")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.blast_index.start()
        self.step.update()

    def stepfinish(self):
        self.step.blast_index.finish()
        self.step.update()

    def check_options(self):
        if not self.option("reference").is_set:
            raise OptionError("必须提供输入文件")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(BlastIndexAgent, self).end()

class BlastIndexTool(Tool):
    def __init__(self, config):
        super(BlastIndexTool, self).__init__(config)
        self.blast_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"
        self.set_environ(PATH=self.blast_path)

    def run(self):
        """
        运行
        :return:
        """
        super(BlastIndexTool, self).run()
        self.blast_index()
        self.end()

    def blast_index(self):
        """
        构建索引文件
        """
        cmd = os.path.join(self.blast_path, "makeblastdb")
        self.db_path = self.option("reference").prop['path']
        cmd += " -dbtype %s -in %s -parse_seqids -out %s " % (self.option('reference_type'),
                                                              self.option("reference").prop['path'],
                                                              self.db_path)
        self.logger.info("开始运行makeblastdb，生成结果库文件放在工作目录的customer_blastdb下")
        makedb_obj = self.add_command("makeblastdb", cmd).run()
        self.wait(makedb_obj)
        if makedb_obj.return_code == 0:
            self.logger.info("makeblastdb运行完成")
        else:
            self.set_error("makeblastdb运行出错!")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "BlastIndex_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.blast_index",
            "instant": False,
            "options": dict(
                reference="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/repeatmasker_test/Trinity.fasta",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()