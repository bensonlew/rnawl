# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class HumanAnnoAgent(Agent):
    """
    author: shenghe
    last_modify: 2017.05.12
    """

    def __init__(self, parent):
        super(HumanAnnoAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_table_dir"}  # 输入文件
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("blastout").is_set:
            raise OptionError("必须设置输入文件")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['pathway_profile.xls', 'xls', 'kegg pathway丰度文件'],
            ['module_profile.xls', 'xls', 'kegg module丰度文件'],
            ])
        super(HumanAnnoAgent, self).end()


class HumanAnnoTool(Tool):
    def __init__(self, config):
        super(HumanAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.humann_scon = self.config.SOFTWARE_DIR + "/bioinfo/annotation/humann-0.99/SConstruct"


    def run(self):
        """
        运行
        :return:
        """
        super(HumanAnnoTool, self).run()
        table_dir = self.work_dir + '/input'
        self.logger.info(self.option("blastout").format)
        if self.option("blastout").format == 'align.blast.blast_xml_dir':
            if not os.path.exists("input"):
                os.makedirs("input")
            files = self.option("blastout").convert2blast6default(table_dir)
        else:
            files = self.option('blastout').prop['files']
        self.humann_anno(files)

    def humann_anno(self, blast_table):
        if not os.path.exists("input"):
            os.makedirs("input")
        for i in blast_table:
            if not os.path.exists("input/" + os.path.basename(i)):
                os.link(i, "input/" + os.path.basename(i))
        with open(self.humann_scon) as f, open("SConstruct", 'w') as w:
            w.write(f.read())
        cmd = "/program/Python/bin/scons " + "--site-dir=" + os.path.dirname(self.humann_scon) + \
            "/site_scons " + "-C " + self.work_dir + " -j 2"
        command = self.add_command("humann", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("humann注释完成")
            self.format_result()
        else:
            self.set_error("humann注释出错")

    def format_result(self):
