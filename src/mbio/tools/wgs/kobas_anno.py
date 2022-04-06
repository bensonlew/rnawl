# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20181009

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class KobasAnnoAgent(Agent):
    """
    diamond注释后的文件再用kobas的annotion.py，即src内的script内程序完成注释后续内容
    """
    def __init__(self, parent):
        super(KobasAnnoAgent, self).__init__(parent)
        options = [
            {"name": "kegg_diamond_result", "type": "string"},  # diamond的kegg的blast结果文件
            # {"name": "sqlite3_path", "type": "string"},     # kobas软件安装sqlite3的db路径
            {"name": "intype", "type": "string", "default": "blastout:tab"},
            {"name": "species", "type": "string", "default": "ko"}
        ]
        self.add_option(options)
        self.step.add_steps('KobasAnno')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.KobasAnno.start()
        self.step.update()

    def step_end(self):
        self.step.KobasAnno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("kegg_diamond_result"):
            raise OptionError("请设置kegg_diamond_result")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 8
        self._memory = '8G'

    def end(self):
        super(KobasAnnoAgent, self).end()


class KobasAnnoTool(Tool):
    def __init__(self, config):
        super(KobasAnnoTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/kobas-3.0")
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/kobas-3.0/src")
        if 'HOME' not in os.environ.keys():
            self.logger.info("环境变量中没有HOME将加载环境变量")
            self.set_environ(HOME=os.path.dirname(self.config.SOFTWARE_DIR))  # 这一步一定要加上，不然找不到.kobasrc，
        else:
            self.logger.info("环境变量中有HOME无需加载环境变量")
        # 会报configuration file does not exist
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/R-3.3.1/lib64/R/lib')
        self.python_path = '/miniconda2/bin/python'
        self.annotation_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/kobas-3.0/scripts/annotate.py"
        self.sqlite3_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/kobas-3.0/sqlite3"

    def KobasAnno(self):
        """
        要重新写下！！！
        :return:
        """
        print os.environ
        name = os.path.basename(self.option('kegg_diamond_result'))
        cmd = "{} {} -i {}".format(self.python_path, self.annotation_path, self.option('kegg_diamond_result'))
        cmd += " -t {} -s {}".format(self.option('intype'), self.option('species'))
        cmd += " -o {} -q {}".format(self.output_dir + "/{}.kobas".format(name), self.sqlite3_path)
        cmd += ' -k {} -v {} -y {} -p {} -x {}'.format(self.config.SOFTWARE_DIR + "/bioinfo/WGS/kobas-3.0",
                                                       self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin",
                                                       self.config.SOFTWARE_DIR + "/bioinfo/WGS/kobas-3.0/seq_pep",
                                                       self.config.SOFTWARE_DIR +
                                                       "/bioinfo/align/ncbi-blast-2.3.0+/bin/blastp",
                                                       self.config.SOFTWARE_DIR +
                                                       "/bioinfo/align/ncbi-blast-2.3.0+/bin/blastx")
        self.logger.info(cmd)
        self.logger.info("开始进行KobasAnno")
        # command = self.add_command("kobasanno", cmd, ignore_error=False).run()  # 必须小写，
        command = self.add_command("kobasanno", cmd, ignore_error=True).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("KobasAnno完成！")
        else:
            if not self.check_error:
                self.set_error("KobasAnno出错！")
            else:
                self.logger.info("KobasAnno没有注释到结果，正常退出！")

    def run(self):
        super(KobasAnnoTool, self).run()
        self.KobasAnno()
        self.end()

    def check_error(self, file_path):
        with open(file_path, 'r') as r:
            for line in r:
                if re.match('None of the entries was annotated successfully.*', line):
                    return True
        return False
