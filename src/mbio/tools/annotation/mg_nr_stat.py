# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last modify by qingchen.zhang @20190410
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.annotation.mg_annotation.mg_taxon import mg_taxon
import os
import threading
import re
import unittest

class MgNrStatAgent(Agent):
    """
    nr_stat.py 等功能 v1.0
    author: zhouxuan
    last_modify: 2017.0913
    last modify by: shaohua.yuan
    """

    def __init__(self, parent):
        super(MgNrStatAgent, self).__init__(parent)
        options = [
            {"name": "taxon_out", "type": "infile", "format": "sequence.profile_table"},# 比对到nr库的结果文件query_taxons.xls
            {"name": "nr_method", "type": "string", "default": "best_hit"}, # nr注释结果选择，eg. best_hit,lca,deunclassified
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self._memory_increase_step = 40  # 每次重运行增加内存20G by guhaidong @ 20190111

    def check_options(self):
        if not self.option("taxon_out").is_set:
            raise OptionError("必须设置输入taxon文件", code="31203101")
        if not self.option("reads_profile_table").is_set:
            raise OptionError("必须设置输入基因丰度文件", code="31203102")
        return True

    def set_resource(self):
        """
        内存改为变动方式 fix byqingchen.zhang@20200403
        将内存20G改为变动内存
        :return:
        """
        self._cpu = 2
        file_size = os.path.getsize(self.option("reads_profile_table").prop['path']) / (1024*1024*1024)   ###(查看文件有多少M)
        memory = int(float(file_size) * 5 + 5)
        if memory < 20 :
            self._memory = "20G"
        elif memory >= 250: ##限制一个正常机器的内存
            self._memory = "250G"
        else:
            self._memory = '{}G'.format(memory)
        # tmp_mem = 10 * (self._rerun_time + 1)  # 每次因拼接失败而重运行的内存增加10G by GHD @ 20180320
        # self._memory = '%sG' % tmp_mem
        # self.logger.info('mg_nr_stat use memory : ' + self._memory)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['query_taxons_detail.xls', 'xls', '序列详细物种分类文件']
        ])
        super(MgNrStatAgent, self).end()


class MgNrStatTool(Tool):
    def __init__(self, config):
        super(MgNrStatTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = "program/Python/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + '/bioinfo/taxon/scripts/mg_nr_profile.py'

    def run(self):
        """
        运行
        :return:
        """
        super(MgNrStatTool, self).run()
        #self.run_nr_stat()
        self.tax_profile()
        self.set_output()
        self.end()

    def run_nr_stat(self):
        #nr = nr_stat()
        nr = mg_taxon()
        self.logger.info("start mg_taxon(detail_to_level)")
        try:
            if "nr_method" in self.get_option_object().keys() and self.option("nr_method") == "best_hit":
                nr.detail_to_level(self.option('taxon_out').prop['path'], self.work_dir)
            elif "nr_method" in self.get_option_object().keys() and self.option("nr_method") == "deunclassified":
                nr.detail_to_level_deunclassified(self.option('taxon_out').prop['path'], self.work_dir)
            elif "nr_method" in self.get_option_object().keys() and self.option("nr_method") == "lca":
                nr.detail_to_level_lca(self.option('taxon_out').prop['path'], self.work_dir)
            else:
                nr.detail_to_level(self.option('taxon_out').prop['path'], self.work_dir) # 兼容正在跑的老项目
        except Exception as e:
            self.set_error("mg_taxon(detail_to_level) failed %s", variables=(e), code="31203101")
            self.set_error("mg_taxon(detail_to_level) failed %s", variables=(e), code="31203103")
        #os.remove(self.work_dir + "/gene_taxons.xls")  # 删除没有用的结果文件
        #os.remove(self.work_dir + "/gene_taxons_detail.xls")
        self.new_query_taxons = os.path.join(self.work_dir, "new_query_taxons.xls")
        self.rm_1(self.work_dir + "/query_taxons.xls", self.new_query_taxons)
        os.remove(self.work_dir + "/query_taxons.xls")
        self.tax_profile()


    def rm_1(self, old_path, new_path):
        with open(old_path, "r") as r, open(new_path, "w") as w:
            for line in r:
                line = line.strip('\n').split("_1\t")
                new_line = ('\t').join(line)
                w.write(new_line + "\n")

    def tax_profile(self):
        self.logger.info("start nr_tax_profile")
        cmd1 = "{} {} -i {} -r {} -o {}".format(self.python_path, self.python_script, self.option("taxon_out").prop['path'],
                                                self.option('reads_profile_table').prop['path'], self.output_dir)
        command1 = self.add_command('tax_profile', cmd1, ignore_error=True).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("tax_profile succeed")
        elif command1.return_code in [-9, 1, -7]:  # change return_code by ghd @ 20190110
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.logger.info("return code: %s" % command1.return_code)
            self.set_error("tax_profile failed", code="31203102")
            self.set_error("tax_profile failed", code="31203104")


    def set_output(self):
        self.logger.info("start set_output")

