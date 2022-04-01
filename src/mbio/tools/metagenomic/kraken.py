# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __last_modify__ = '20191227'


from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
from mbio.packages.toolapps.common import link_dir


class KrakenAgent(Agent):
    """
    Kraken基于reads进行丰度统计和功能计算
    """

    def __init__(self, parent):
        super(KrakenAgent, self).__init__(parent)
        options = [
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},
            {"name": "sample", "type": "string"},
            {"name": "readlen", "type": "int", "default": 100},
            {"name": "confidence", "type": "float", "default": 0.5},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("read1").is_set:
            raise OptionError("请输入read1文件！")
        if not self.option("read2").is_set:
            raise OptionError("请输入read2文件！")
        if not self.option("sample"):
            raise OptionError("请输入sample的样品名称！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 8
        self._memory = '50G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
        ])
        super(KrakenAgent, self).end()

class KrakenTool(Tool):
    def __init__(self, config):
        super(KrakenTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/kraken2-2.0.8-beta/bin:" + self.config.SOFTWARE_DIR + "/gcc/7.2.0/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/nhmmer:"
        self.lib = self.config.SOFTWARE_DIR + "/gcc/7.2.0/lib64:"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.kraken2 = "/bioinfo/metaGenomic/kraken2-2.0.8-beta/bin/kraken2"
        self.out = self.work_dir + "/" + self.option("sample") + ".taxon.xls"
        self.db =self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/kraken2-2.0.8-beta/database"
        self.brakren = "/bioinfo/metaGenomic/Bracken-master/"
        self.kraken = "/bioinfo/metaGenomic/kraken2-2.0.8-beta/bin/"

    def run_kraken2(self):
        """
        description
        :return:
        """
        cmd = "{} --threads 8 --db {} --report {} --confidence {} --output {}".format(self.kraken2, self.db, self.out, self.option("confidence"), self.work_dir + "/" + self.option("sample") + ".mapping.stdout")
        cmd +=" --paired {} {}".format(self.option("read1").prop['path'], self.option("read2").prop['path'])
        command = self.add_command("run_kraken2", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            os.system("grep '^C' {} >{}".format(self.work_dir + "/" + self.option("sample") + ".mapping.stdout", self.output_dir + "/" + self.option("sample") + ".mapping.xls"))
            self.logger.info("run_kraken2运行完成！")
        else:
            self.set_error("run_kraken2运行完成运行出错!")

    def run_brakren_build(self):
        """
        构建brakren的索引
        """
        cmd = "{}bracken-build -d {} -t 8 -x {} -l {} ".format(self.brakren, self.db, self.kraken, self.option("readlen"))
        command = self.add_command("run_brakren_build", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_brakren_build运行完成！")
        else:
            self.set_error("run_brakren_build运行完成运行出错!")

    def run_brakren(self):
        """
        运行brakren的，得到物种各水平分类
        """
        if os.path.exists(self.work_dir + "/result"):
            shutil.rmtree(self.work_dir + "/result")
        os.mkdir(self.work_dir + "/result")
        for i in ["D","P","C","O","F","G","S"]:
            cmd = "{}bracken -i {} -d {} -o {} -r {} -l {} -t 0".format(self.brakren, self.out, self.db, self.work_dir + "/result/" + self.option("sample") + "_" + i + ".xls", self.option("readlen"), i)
            psd = "run_brakren_" + i.lower()
            command = self.add_command(psd, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("{}运行完成!".format(psd))
            else:
                self.set_error("{}运行完成运行出错!".format(psd))

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        if len(os.listdir(self.output_dir)) >=1:
            shutil.rmtree(self.output_dir)
        link_dir(self.work_dir + "/result", self.output_dir)
        os.link(self.out, self.output_dir + "/" + self.option("sample") + ".taxon.xls")

    def run(self):
        super(KrakenTool, self).run()
        self.run_kraken2()
        self.run_brakren_build()
        self.run_brakren()
        self.set_output()
        self.end()