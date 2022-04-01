# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re

class VfdbAnnoStatAgent(Agent):
    """
    宏基因vfdb注释结果丰度统计表
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(VfdbAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "vfdb_core_anno", "type": "infile", "format": "sequence.profile_table"},
            # 核心库注释表
            {"name": "vfdb_predict_anno", "type": "infile", "format": "sequence.profile_table"},
            #预测数据库注释表
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"}
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("vfdb_core_anno").is_set:
            raise OptionError("找不到核心注释文件", code="31203801")
        if not self.option("vfdb_predict_anno").is_set:
            raise OptionError("找不到预测注释文件", code="31203802")
        if not self.option('reads_profile_table').is_set:
            raise OptionError("必须设置基因丰度文件", code="31203803")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'  # 内存2G增加到10G  by GHD @20180309 改回 by GHD @20180428
        # tmp_mem = 5 * (self._rerun_time + 1)  # 每次因拼接失败而重运行的内存增加5G by GHD @ 20180320
        # self._memory = '%sG' % tmp_mem
        # self.logger.info('vfdb_anno_stat use memory : ' + self._memory)

    def end(self):
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['query_taxons_detail.xls', 'xls', '序列详细物种分类文件']
            ])
        """
        super(VfdbAnnoStatAgent, self).end()

class VfdbAnnoStatTool(Tool):
    def __init__(self, config):
        super(VfdbAnnoStatTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.script = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/vfdb_anno_abudance.pl'

    def run(self):
        """
        运行
        :return:
        """
        super(VfdbAnnoStatTool, self).run()
        self.run_total_anno()
        self.run_vfdb_stat()
        self.set_output()
        self.end()

    def run_total_anno(self):
        self.logger.info("start merge core and predict anno table.")
        core_anno = self.option('vfdb_core_anno').prop['path']
        predict_anno = self.option('vfdb_predict_anno').prop['path']
        outfile = self.output_dir + "/gene_vfdb_total_anno.xls"
        with open(core_anno, "r") as f1, open(predict_anno, "r") as f2, open(outfile, "w") as outf:
            for line in f1:
                line = line.strip()
                line1 = line.split("\t")
                if line1[0] == "#Query":
                    head = line
                    outf.write(head + "\t" + "Database\n")
                else:
                    outf.write(line + "\tcore\n")
            for line in f2:
                line = line.strip()
                line1 = line.split("\t")
                if line1[0] != "#Query":
                    outf.write(line + "\tpredict\n")
        self.logger.info("finish merge core and predict anno table.")

    def run_vfdb_stat(self):
        self.logger.info("start vfdb_stat")
        cmd = "{} {} -c {} -pre {} -p {} -o {}".format(self.perl_path, self.script, self.option('vfdb_core_anno').prop['path'],self.option('vfdb_predict_anno').prop['path'],self.option('reads_profile_table').prop['path'], self.output_dir)
        command = self.add_command('vfdb_profile', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("vfdb_stat succeed")
        elif command.return_code == -9:  # change code by guhaidong @ 20190110
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.set_error("vfdb_stat failed", code="31203801")
            raise Exception("vfdb_stat failed")

    def set_output(self):
        if len(os.listdir(self.output_dir)) == 8:
            self.logger.info("结果文件正确生成")
        else:
            self.set_error("文件个数不正确，请检查", code="31203802")
            raise Exception("文件个数不正确，请检查")
