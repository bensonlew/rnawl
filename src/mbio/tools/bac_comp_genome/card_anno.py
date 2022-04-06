# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import subprocess
import os


class CardAnnoAgent(Agent):
    """
    比较基因组比对结果文件注释,先将xml文件转换成table并合并成一张table，进行anno详细注释
    author: gaohao
    last_modify:
    """

    def __init__(self, parent):
        super(CardAnnoAgent, self).__init__(parent)
        options = [
            {"name": "card_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # 比对到card库的xml文件夹
            {"name": "card_anno_result", "type": "outfile", 'format': "sequence.profile_table"},  # 注释详细结果表
            ]
        self.add_option(options)
        self.result_name = ''

    def check_options(self):
        if not self.option("card_xml_dir").is_set:
            raise OptionError("必须设置输入文件")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(CardAnnoAgent, self).end()


class CardAnnoTool(Tool):
    def __init__(self, config):
        super(CardAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/python"
        self.db = self.config.SOFTWARE_DIR + "/database/card/card_function.xls"
        self.python_script = self.config.PACKAGE_DIR + '/bac_comp_genome/anno_card.py'
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.result_name = ''

    def run(self):
        """
        运行
        :return:
        """
        super(CardAnnoTool, self).run()
        self.set_output()
        self.merge_table()
        self.run_card_anno()
        self.end()

    def merge_table(self):
        self.card_number = 0
        xml_file = os.listdir(self.option('card_xml_dir').prop['path'])
        self.result_name = os.path.join(self.output_dir, "card_align_table.xls")
        if os.path.exists(self.result_name):
            os.remove(self.result_name)
        for i in xml_file:
            self.card_number += 1
            file_path = os.path.join(self.option('card_xml_dir').prop['path'], i)
            table = xml2table(file_path,
                              self.work_dir + "/tmp_card_anno/" + "card_" + str(self.card_number) + "_table.xls")
            cmd = '{} {} {}'.format(self.sh_path, table, self.result_name)
            self.logger.info("start cat {}".format(i))
            command_name = "cat" + str(self.card_number)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(i))
            else:
                self.set_error("cat %s error", variables=(i))
        with open(self.result_name, "r") as f:
            head = f.next().strip()
        os.system("sed -i '/^Score\t/d ' " + self.result_name)
        os.system("sed -i '1i" + head + "\' " + self.result_name)

    def run_card_anno(self):
        cmd = '{} {} -f {} -i {} -o {}'.format(self.python_path, self.python_script, self.db, self.result_name , self.output_dir + "/gene_card_anno.xls")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行card_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行card_anno出错', code="31200502")

    def set_output(self):
        if os.path.exists(self.work_dir + '/tmp_card_anno'):
            pass
        else:
            os.mkdir(self.work_dir + '/tmp_card_anno')