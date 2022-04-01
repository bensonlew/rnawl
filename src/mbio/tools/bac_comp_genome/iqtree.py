# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2019.10.14

import re,os,shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.bac_comp_genome.common_function import link_file


class IqtreeAgent(Agent):
    """
    iqtree进行进化树构建的tool
    """

    def __init__(self, parent):
        super(IqtreeAgent, self).__init__(parent)
        options = [
            {"name": "fa", "type": "infile", "format": "sequence.fasta"},  # 所有16s的序列
            {"name": "out_group", "type": "string"}, #外类群在进化文件中名称
            {"name": "bootstrap", "type": "int", "default": "1000"},  # 外类群在进化文件中名称
            {"name": "merit", "type": "string", "default": "BIC"},#模型选择评估标准，BIC、AIC和AICc
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fa").is_set:
            raise OptionError("必须输入fa比对文件！")
        if self.option("bootstrap") < 1000:
            raise OptionError("输入bootstrap必须大于1000！")
        return True

    def set_resource(self):
        self._cpu = 8
        self._memory = '20G'

    def end(self):
        super(IqtreeAgent, self).end()


class IqtreeTool(Tool):
    def __init__(self, config):
        super(IqtreeTool, self).__init__(config)
        self._version = "1.0"
        self.sof_path ="/bioinfo/compare_genome/software/iqtree-1.6.12-Linux/bin/iqtree"
        self.fa = self.option("fa").prop["path"]
        self.bootstrap = self.option("bootstrap")
        self.merit = self.option("merit")

    def run(self):
        """
        运行
        :return:
        """
        super(IqtreeTool, self).run()
        self.run_tree()
        self.set_output()
        self.end()

    def run_tree(self):
        if os.path.exists(self.work_dir + "/temp"):
            shutil.rmtree(self.work_dir + "/temp")
        os.mkdir(self.work_dir + "/temp")
        cmd = '{} -s {} -pre {} -m MFP -nt AUTO -ntmax 8 -bb {} -merit {}'.format(self.sof_path, self.fa, self.work_dir + "/temp/all", self.bootstrap, self.merit)
        if self.option("out_group"):
            cmd += " -o {}".format(self.option("out_group"))
        command = self.add_command("run_tree", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.chuli_file(self.work_dir + "/temp/all.iqtree", self.work_dir + "/all.assess_result.xls")
            self.logger.info("run_tree运行完成！")
        else:
            self.set_error("run_tree运行完成运行出错!")

    def set_output(self):
        link_file(self.work_dir + "/temp/all.treefile", self.output_dir + "/all.tree.nwk")
        link_file(self.work_dir + "/all.assess_result.xls", self.output_dir + "/all.assess_result.xls")

    def get_num(self, file):
        list = []
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if re.search('^>', line):
                    list.append(line)
        return len(list)

    def chuli_file(self, input, output):
        with open(input, 'r') as f, open(output, "w") as g:
            lines = f.readlines()
            n = 0
            j = 1
            pat = re.compile("\s*")
            for line in lines:
                if re.search("^Model\s+LogL", line):
                    n += 1
                elif re.search("^AIC,\sw-AIC", line):
                    n = 0
                    break
                if n == 1 and not re.search("^\n", line):
                    lin = pat.split(line.strip())
                    if j == 1:
                        g.write("\t".join(lin) + "\n")
                    else:
                        g.write(lin[0] + "\t" + lin[1] + "\t" + lin[2] + "\t" + lin[2] + lin[3] + lin[4] + "\t" + lin[
                            5] + "\t" + lin[5] + lin[6] + lin[7] + "\t" + lin[8] + "\t" + lin[8] + lin[9] + lin[
                                    10] + "\n")
                    j += 1
