#-*- coding: utf-8 -*-

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bacgenome.common import get_num,link_file


class IntegronPredictAgent(Agent):
    """
    细菌基因组，进行integron的预测
    """
    def __init__(self, parent):
        super(IntegronPredictAgent, self).__init__(parent)
        options = [
            {"name": "input_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "file_dir", "type": "string"},  # 该文件夹下有*.gff,*.faa,*.fna(基因的gff文件，基因蛋白文件，基因的核酸文件)
            {"name": "sample", "type": "string"}, ## 传入的样本名称
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_fa").is_set:
            raise OptionError("必须设置参数input_fa文件!")

    def set_resource(self):
        self._cpu = 4
        self._memory = '10G'

    def end(self):
        super(IntegronPredictAgent, self).end()


class IntegronPredictTool(Tool):
    def __init__(self, config):
        super(IntegronPredictTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/hmmer-3.1b2/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/infernal-1.1.3-linux-intel-gcc/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/mobile_genetic_elements/Integron_Finder-master/data:"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib:"
        self.set_environ(PATH=self.path,LD_LIBRARY_PATH=self.lib)
        self.python = "/program/Python35/bin/integron_finder"
        self.names = self.option("sample")
        self.fasta = self.work_dir + "/" + self.names
        self.python2 = "/program/Python/bin/python"
        self.packages = self.config.PACKAGE_DIR + "/bacgenome/tiqu_sequence.py"
        self.stat = self.work_dir + "/Results_Integron_Finder_" + self.names+ "/"+ self.names +".integrons"

    def run_integron(self):
        """
        运行Integron_Finder软件进行预测
        :return:
        """
        if os.path.exists(self.work_dir + "/" + self.names):
            os.remove(self.work_dir + "/" + self.names)
        os.link(self.option("input_fa").prop["path"],self.work_dir + "/" + self.names)
        if os.path.exists(self.work_dir + "/Results_Integron_Finder_" + self.names):
            shutil.rmtree(self.work_dir + "/Results_Integron_Finder_" + self.names)
        if self.option("file_dir"):
            shutil.copytree(self.option("file_dir"), self.work_dir + "/Results_Integron_Finder_" + self.names)
        cmd = "{} --local-max --func-annot --keep-tmp --cpu 4 --outdir {} {}".format(self.python, self.work_dir, self.fasta)
        self.logger.info(cmd)
        command = self.add_command("run_integron", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_integron运行完成！")
        else:
            self.set_error("run_integron运行完成运行出错!")

    def run_stat(self):
        """
        运行脚本得到序列和完成统计
        :return:
        """
        stat_file = self.work_dir + "/Results_Integron_Finder_" + self.names+ "/"+ self.names +".summary"
        out_stat = os.path.join(self.work_dir + "/Results_Integron_Finder_" + self.names,self.names +".stat.xls")
        cmd = "{} {} {} {} {} {} {}".format(self.python2,self.packages, self.output_dir + "/" + self.names + ".integrons", self.fasta, self.names, stat_file, out_stat)
        self.logger.info(cmd)
        command = self.add_command("run_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("tiqu_sequence运行完成！")
        else:
            self.set_error("tiqu_sequence运行完成运行出错!")

    def integron_stat(self):
        scaffload_list =[]
        scaffload_num = 0
        with open(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names + ".integrons","r") as v,open(self.output_dir + "/" + self.names + ".integrons","w") as t:
            data = v.readlines()
            if len(data) > 1:
                t.write(data[0])
                t.write(data[1])
                for i in data[2:]:
                    if i:
                        if i.split("\t")[1] not in scaffload_list:
                            scaffload_list.append(i.split("\t")[1])
                            scaffload_num += 1
                        t.write("integron_" + str(scaffload_num).zfill(2) + "\t" + "\t".join(i.strip().split("\t")[1:]) + "\n" )

    def set_output(self):
        self.logger.info(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names +".integrons")
        num =get_num(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names +".integrons")
        if num >2:
            #if os.path.exists(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names + ".integrons"):
            #   link_file(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names + ".integrons", self.output_dir + "/" + self.names + ".integrons")
            if os.path.exists(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names + ".summary"):
                link_file(self.work_dir + "/Results_Integron_Finder_" + self.names + "/" + self.names + ".summary", self.output_dir + "/" + self.names + ".summary")
            if os.path.exists(self.work_dir + "/"+ self.names + "_sequence.fna"):
                link_file(self.work_dir + "/"+ self.names + "_sequence.fna", self.output_dir + "/"+ self.names + "_sequence.fna")
            if os.path.exists(self.work_dir + "/"+ self.names + ".stat.xls"):
                link_file(self.work_dir + "/"+ self.names + ".stat.xls", self.output_dir + "/"+ self.names + ".stat.xls")
            if os.path.exists(self.work_dir + "/"+ self.names + ".sample.xls"):
                link_file(self.work_dir + "/"+ self.names + ".sample.xls", self.output_dir + "/"+ self.names + ".sample.xls")
            if os.path.exists(self.work_dir + "/"+ self.names + ".integron.fna"):
                link_file(self.work_dir + "/"+ self.names + ".integron.fna", self.output_dir + "/"+ self.names + ".integron.fna")

    def run(self):
        super(IntegronPredictTool, self).run()
        self.run_integron()
        if os.path.exists(self.stat) and get_num(self.stat) > 2:
            self.integron_stat()
            self.run_stat()
        self.set_output()
        self.end()

