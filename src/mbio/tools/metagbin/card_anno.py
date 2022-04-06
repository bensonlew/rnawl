# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import subprocess
import os,time
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout


class CardAnnoAgent(Agent):
    """
    宏基因组比对结果文件注释,先将xml文件转换成table并合并成一张table，进行anno详细注释
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(CardAnnoAgent, self).__init__(parent)
        options = [
            {"name": "card_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # 比对到card库的xml文件夹
            {"name": "card_anno_result", "type": "outfile", 'format': "sequence.profile_table"},  # 注释详细结果表
            {"name": "result_format", "type": "string", 'default': "meta"},  # 输出结果类型，meta或者dna add by zouxuan
            {"name": "category", "type": "string" ,"default":""},  #"default": "aro_category"
            {"name": "database", "type": "string", "default": "card_v3.0.9"}, ## 增加此字段是为了加版本信息
            # 结果类型，aro_category or drug_class 细菌升级使用
            ]
        self.add_option(options)
        self.result_name = ''

    def check_options(self):
        if not self.option("card_xml_dir").is_set:
            raise OptionError("必须设置输入文件", code="31200501")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        if os.path.exists(os.path.join(self.output_dir,"gene_card_anno.xls")):
            self.option('card_anno_result',os.path.join(self.output_dir,"gene_card_anno.xls"))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['gene_card_anno.xls', 'xls', '序列详细物种分类文件']
            ])
        super(CardAnnoAgent, self).end()


class CardAnnoTool(Tool):
    def __init__(self, config):
        super(CardAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/python"
        #self.python_path = "miniconda2/bin/python"
        # self.python_script = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/meta_card_mongo.py'
        self.python_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/mg_anno_card.py' ##fix_by qingchen.zhang @20200812 card
        self.python_dna_script = self.config.PACKAGE_DIR + '/metagbin/dna_card.py' ## 老的
        self.python_dna_script2 = self.config.PACKAGE_DIR + '/annotation/dna_card.py' ## 新的
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.result_name = ''
        self.process_rerun = 0

    def run(self):
        """
        运行
        :return:
        """
        super(CardAnnoTool, self).run()
        self.set_output()
        self.merge_table()
        self.logger.info(self.num)
        if self.num > 1:
            self.run_card_anno()
        self.end()

    def get_num(self):
        with open(self.result_name, "r") as f:
            lines = f.readlines()
            return len(lines)

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
            #if self.card_number > 1:
            #    os.system("sed -i '/^Score\t/d ' " + table)
            cmd = '{} {} {}'.format(self.sh_path, table, self.result_name)
            self.logger.info("start cat {}".format(i))
            command_name = "cat" + str(self.card_number)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(i))
            else:
                self.set_error("cat %s error", variables=(i), code="31200501")
                raise Exception("cat {} error".format(i))
        with open(self.result_name, "r") as f:
            head = f.next().strip()
        os.system("sed -i '/^Score\t/d ' " + self.result_name)
        os.system("sed -i '1i" + head + "\' " + self.result_name)
        self.num = self.get_num()
        self.logger.info(self.num)

    def run_card_anno(self):
        if "result_format" not in self.get_option_object().keys() or self.option("result_format") == "meta":
            script_path = self.python_script
        else:
            if self.option('category') != "drug_class": ## 老的
                script_path = self.python_dna_script
            else: ## 新的
                script_path = self.python_dna_script2
        cmd = '{} {} -i {} -o {}'.format(self.python_path, script_path,self.result_name , self.output_dir + "/gene_card_anno.xls")
        #command = self.add_command("anno", cmd).run()

        if self.option('category') != "":
        #if "category" in self.get_option_object().keys():
            cmd += ' -c {} '.format(self.option("category"))
        self.logger.info(cmd)
        #self.wait(command)
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行card_anno完成')
        except (ServerSelectionTimeoutError, NetworkTimeout): # 捕获因为mongo服务器问题导致的异常后重运行此方法 fix by # qingchen.zhang@20210121
            if self.process_rerun < 5:
                self.process_rerun += 1
                self.logger.info("检测到Time out Error, 第{}次重运行方法".format(self.process_rerun))
                time.sleep(5)
                self.run_card_anno()
            else:
                self.add_state('memory_limit', '检测到Time out Error, 重运行tool')
        except subprocess.CalledProcessError:
            self.set_error('运行card_anno出错', code="31200502")

    def set_output(self):
        if os.path.exists(self.work_dir + '/tmp_card_anno'):
            pass
        else:
            os.mkdir(self.work_dir + '/tmp_card_anno')
            #self.option('card_anno_result',os.path.join(self.output_dir,"gene_card_anno.xls"))
            #if len(os.listdir(self.output_dir)) == 1:
            #    self.logger.info("output right")
