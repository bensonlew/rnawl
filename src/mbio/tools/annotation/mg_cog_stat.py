# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class MgCogStatAgent(Agent):
    """
    宏基因cog注释结果统计.py v1.0
    author: zhouxuan
    last_modify: 2017.09.14
    last_modify by :shaohua.yuan
    """

    def __init__(self, parent):
        super(MgCogStatAgent, self).__init__(parent)
        options = [
            {"name": "cog_table_dir", "type": "infile", "format": "annotation.mg_anno_dir"},
            # 比对到eggNOG库的注释结果文件
            {"name": "align_table_dir", "type": "infile", "format": "annotation.mg_anno_dir"},
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self._memory_increase_step = 40  # 每次重运行增加内存20G by qingchen.zhang @ 20190529

    def check_options(self):
        if not self.option("cog_table_dir").is_set:
            raise OptionError("必须设置输入cog注释文件", code="31202701")
        if not self.option('reads_profile_table').is_set:
            raise OptionError("必须设置输入丰度文件", code="31202702")
        return True

    def set_resource(self):
        """
        内存改为变动方式 fix byqingchen.zhang@20200403
        将内存20G改为变动内存
        :return:
        """
        self._cpu = 5
        file_size = os.path.getsize(self.option("reads_profile_table").prop['path']) / (1024*1024*1024)   ###(查看文件有多少G)
        memory = int(float(file_size) * 5 + 20)
        if memory < 20 :
            self._memory = "20G"
        elif memory >= 250: ##限制一个正常机器的内存
            self._memory = "250G"
        else:
            self._memory = '{}G'.format(memory)
        # self._memory = '20G'  # 改回 by GHD @ 20180428
        # tmp_mem = 10 * (self._rerun_time + 1)  # 每次因拼接失败而重运行的内存增加10G by GHD @ 20180320
        # self._memory = '%sG' % tmp_mem
        # self.logger.info('mg_cog_stat use memory : ' + self._memory)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(MgCogStatAgent, self).end()

class MgCogStatTool(Tool):
    def __init__(self, config):
        super(MgCogStatTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        # self.script = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/eggNOG_anno_abundance.pl'
        self.script = self.config.PACKAGE_DIR + '/metagenomic/anno_profile/cog_profile.pl'
        self.result_name = ''
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.cat_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/cat_seq.sh'##add by qingchen.zhang@20190529

    def run(self):
        """
        运行
        :return:
        """
        super(MgCogStatTool, self).run()
        if self.option("align_table_dir").is_set:
            self.merge_align_table()
        self.merge_table()
        self.run_cog_stat()
        self.set_output()
        self.end()

    def merge_align_table(self):
        """
        合并所有cog比对结果
        :return:
        """
        profile_file = os.listdir(self.option('align_table_dir').prop['path'])
        self.align_table = os.path.join(self.output_dir, "cog_align_table.xls")
        if os.path.exists(self.align_table):
            os.remove(self.align_table)
        cmd = '{}'.format(self.cat_path)
        for i in profile_file:
            file_path = os.path.join(self.option('align_table_dir').prop['path'], i)
            with open(file_path, 'r') as f:
                line = f.readline()
                if re.search(r'^Score\t', line):
                    os.system("sed -i '/^Score\t/d ' " + file_path)
            cmd += ' ' + file_path
            #cmd = '{} {} {}'.format(self.sh_path, file_path, self.align_table)
        cmd += ' ' + self.align_table ##add by qingchen.zhang@20190529
        self.logger.info("start cat align")
        command_name = "cat_align"
        command = self.add_command(command_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{} done".format(command_name))
        else:
            self.set_error("cat align error", code="31202701")
            raise Exception("cat align error")
        os.system("sed -i '/^Score\t/d ' "+ self.align_table)
        os.system("sed -i '1iScore\tE-Value\tHSP-Len\tIdentity-%\tSimilarity-%\tQuery-Name\tQ-Len\tQ-Begin\tQ-End\tQ-Frame\tHit-Name\tHit-Len\tHsp-Begin\tHsp-End\tHsp-Frame\tHit-Description'" + " %s" %(self.align_table))

    def merge_table(self):
        """
        合并所有注释结果文件
        :return:
        """
        profile_file = os.listdir(self.option('cog_table_dir').prop['path'])
        self.result_name = os.path.join(self.output_dir, "gene_cog_anno.xls")
        if os.path.exists(self.result_name):
            os.remove(self.result_name)
        cmd = '{}'.format(self.cat_path)
        for i in profile_file:
            file_path = os.path.join(self.option('cog_table_dir').prop['path'], i)
            #if cog_number > 1:
            #with open(file_path, 'r') as f:
                #line = f.readline()
                #if re.search(r'^#', line):
                    #os.system("sed -i '/^#/d ' " + file_path)
            cmd += ' ' + file_path
            #cmd = '{} {} {}'.format(self.sh_path, file_path, self.result_name)
        cmd += ' ' + self.result_name ##add by qingchen.zhang@20190529
        self.logger.info("start cat align")
        command_name = "cat_anno"
        command = self.add_command(command_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cat_anno done")
        else:
            self.set_error("cat anno error", code="31202702")
            raise Exception("cat anno error")
        os.system("sed -i '/^#/d ' "+ self.result_name)
        os.system("sed -i '1i#Query\tNOG\tNOG_description\tFunction\tFun_description\tCategory\tIdentity(%)\tAlign_len'" + " %s" %(self.result_name))


    def run_cog_stat(self):
        """
        进行功能水平的丰度计算
        :return:
        """
        self.logger.info("start cog_stat")
        # cmd = "{} {} -q {} -p {} -o {}".format(self.perl_path, self.script, self.result_name,
        #                                        self.option('reads_profile_table').prop['path'], self.output_dir)
        cmd = "{} {} {} {} {}".format(self.perl_path, self.script, self.result_name,
                                      self.option('reads_profile_table').path, self.output_dir)
        command = self.add_command('tax_profile', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cog_stat succeed")
        elif command.return_code in [137, -9]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("cog_stat failed", code="31202703")
            raise Exception("cog_stat failed")

    def set_output(self):
        self.logger.info("start set_output")
        for f in os.listdir(self.output_dir):  # 删除sed的中间文件
            if f.startswith('sed'):
                fp = os.path.join(self.output_dir, f)
                os.system("rm -f " + fp)
