# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __last_modify__ = '20200109'


from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import subprocess
from mbio.packages.toolapps.common import link_dir

class CentrifugeAgent(Agent):
    """
    Centrifuge基于reads进行物种分类统计
    """
    def __init__(self, parent):
        super(CentrifugeAgent, self).__init__(parent)
        options = [
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},
            {"name": "sample", "type": "string"},
            {"name": "readlen", "type": "int", "default": 100},
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
        if os.path.getsize(self.option("read1").prop['path']) <=10000000000:
            self._cpu = 8
            self._memory = '30G'
        else:
            num = os.path.getsize(self.option("read1").prop['path']) / 1000000000
            self._cpu = 10
            self._memory = str(num * 4) + 'G'

    def end(self):
        super(CentrifugeAgent, self).end()

class CentrifugeTool(Tool):
    def __init__(self, config):
        super(CentrifugeTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/kraken2-2.0.8-beta/bin:" + self.config.SOFTWARE_DIR + "/gcc/7.2.0/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/nhmmer:"
        self.lib = self.config.SOFTWARE_DIR + "/gcc/7.2.0/lib64:"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.centrifuge = "/bioinfo/metaGenomic/centrifuge-1.0.4-beta/bin/centrifuge"
        self.centrifuge_kreport = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/centrifuge-1.0.4-beta/bin/centrifuge-kreport"
        self.db = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/centrifuge-1.0.4-beta/indices/data/p_compressed"
        self.brakren = "/bioinfo/metaGenomic/Bracken-master/"
        self.kraken = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/kraken2-2.0.8-beta/bin/"
        self.kraken2_db = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/kraken2-2.0.8-beta/database/"

    def run_centrifuge(self):
        """
        description
        :return:
        """
        cmd = "{} -q -p 8 -x {} -1 {} -2 {}".format(self.centrifuge, self.db, self.option("read1").prop['path'], self.option("read2").prop['path'])
        command = self.add_command("run_centrifuge", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_centrifuge运行完成！")
        else:
            self.set_error("run_centrifuge运行完成运行出错!")

    def run_kreport(self):
        """
        解析centrifuge软件的结果，其结果类似kraken2结果
        """
        cmd = "{} -x {} {} > {}".format(self.centrifuge_kreport, self.db, self.work_dir + "/centrifuge_report.tsv", self.work_dir + "/report.tsv")
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("run_kreport done")
        except subprocess.CalledProcessError:
            self.set_error("run_kreport error")

    def run_brakren_build(self):
        """
        构建brakren的索引
        """
        cmd = "{}bracken-build -d {} -t 8 -x {} -l {} ".format(self.brakren, self.kraken2_db, self.kraken, self.option("readlen"))
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
        for i in ["D", "P", "C", "O", "F", "G", "S"]:
            cmd = "{}bracken -i {} -d {} -o {} -r {} -l {} -t 0".format(self.brakren, self.work_dir + "/report.tsv", self.kraken2_db, self.work_dir + "/result/" + self.option("sample") + "_" + i + ".xls",self.option("readlen"), i)
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
            for i in os.listdir(self.output_dir):
                os.remove(self.output_dir + "/" + i)
        link_dir(self.work_dir + "/result", self.output_dir)
        os.link(self.work_dir + "/centrifuge_report.tsv", self.output_dir + "/" + self.option("sample") +".centrifuge_report.tsv")

    def run(self):
        super(CentrifugeTool, self).run()
        self.run_centrifuge()
        self.run_kreport()
        self.run_brakren_build()
        self.run_brakren()
        self.set_output()
        self.end()