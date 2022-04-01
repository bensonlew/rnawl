# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import subprocess
import pandas as pd
import os,re


class MgKeggLevelAgent(Agent):
    """
    宏基因kegg注释level水平丰度计算
    author: zhouxuan
    last_modify: 2017.09.25
    last_modify by: shaohua.yuan
    """

    def __init__(self, parent):
        super(MgKeggLevelAgent, self).__init__(parent)
        options = [
            {"name": "kegg_result_dir", "type": "infile", "format": "annotation.mg_anno_dir"},
            {"name": "reads_profile", "type": "infile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self._memory_increase_step = 20  # 每次重运行增加内存20G by qingchen.zhang @ 20190529

    def check_options(self):
        if not self.option("kegg_result_dir").is_set:
            raise OptionError("必须设置输入文件", code="31202901")
        if not self.option("reads_profile").is_set:
            raise OptionError("必须设置基因丰度表", code="31202902")
        return True

    def set_resource(self):
        """
        内存改为变动方式 fix byqingchen.zhang@20200403
        将内存5G转为变动内存
        :return:
        """
        self._cpu = 2
        file_size = os.path.getsize(self.option("reads_profile").prop['path']) / (1024*1024*1024)   ###(查看文件有多少G)
        memory = int(float(file_size) * 5 + 20)
        if memory < 20 :
            self._memory = "20G"
        elif memory >= 250: ##限制一个正常机器的内存
            self._memory = "250G"
        else:
            self._memory = '{}G'.format(memory)
        # tmp_mem = 5 + 10 * self._rerun_time  # 每次投递加10G内存
        # self._memory = "%sG" % tmp_mem

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(MgKeggLevelAgent, self).end()


class MgKeggLevelTool(Tool):
    def __init__(self, config):
        super(MgKeggLevelTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        # self.perl_script = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/mg_kegg_level_abu.pl'
        self.perl_script = self.config.PACKAGE_DIR + '/metagenomic/anno_profile/kegg_level_profile.pl'
        self.python_path = "program/Python/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/mg_kegg_level_mongo.py'
        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.cat_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/cat_seq.sh'
        self.anno_result = ''
        self.level_anno = ''

    def run(self):
        """
        运行
        :return:
        """
        super(MgKeggLevelTool, self).run()
        self.merge_table()
        self.run_kegg_level_anno()
        self.run_kegg_stat()
        self.set_output()
        self.end()

    def merge_table(self):
        """
        合并注释结果文件
        :return:
        """
        profile_file = os.listdir(self.option('kegg_result_dir').prop['path'])
        self.anno_result = os.path.join(self.work_dir, "tmp_kegg_anno.xls")
        if os.path.exists(self.anno_result):
            os.remove(self.anno_result)
        cmd = '{}'.format(self.cat_path)
        for i in profile_file:
            if "kegg_anno_result" in i:
                file_path = os.path.join(self.option('kegg_result_dir').prop['path'], i)
                #cmd = '{} {} {}'.format(self.sh_path, file_path, self.anno_result)
                with open(file_path, 'r') as f:
                    line = f.readline()
                    if re.search(r'^#', line):
                        os.system("sed -i '/^#\t/d ' " + file_path)
                cmd += ' ' + file_path
        cmd += ' ' + self.anno_result ##add by qingchen.zhang@20190529
        self.logger.info("start cat_anno")
        command_name = "cat_anno"
        command = self.add_command(command_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cat {} done".format(command_name))
        else:
            self.set_error("cat_anno error", code="31202901")
            raise Exception("cat_anno error")
        os.system("sed -i '/^#\t/d ' "+ self.anno_result)
        os.system("sed -i '1i#Query\tGene\tKO\tDefinition\tPathway\tEnzyme\tModule\tHyperlink\tIdentity(%)\tAlign_len'" + " %s" %(self.anno_result))

    def run_kegg_level_anno(self):
        self.level_anno =  self.work_dir + "/kegg_level_anno.xls"
        anno_all = self.work_dir + "/gene_kegg_anno_all.xls"
        ko_path = self.config.SOFTWARE_DIR + "/database/KEGG/metag_database/kegg_v94.2/kegg_ko_v94.2.xls"
        info_col = ["ko_name", "level1", "level2", "level3"]
        self.ko_info = pd.read_csv(ko_path, sep='\t', index_col=0)[info_col]
        with open(self.anno_result, 'r') as infile, open(self.level_anno, 'wb') as outfile, open(anno_all, "wb") as outfile2:
            outfile.write('#Query\tlevel1\tlevel2\tlevel3\n')
            head = infile.next().strip()
            outfile2.write(head + "\t" +"KEGG_Name\tLevel1\tLevel2\tLevel3\n")
            for line in infile:
                line = line.strip()
                line1 = line.split("\t")
                if not "#Query" in line1[0]:
                    outfile2.write(line + "\t")
                    query = line1[0]
                    KO = line1[2]
                    if KO in self.ko_info.index:
                        ko_name, le1, le2, le3 = list(self.ko_info.loc[KO, ])
                        le1 = re.sub(r'^ ', '', le1)
                        le2 = re.sub(r'^ ', '', le2)
                        le3 = re.sub(r'^ ', '', le3)
                        outfile2.write("{}\t{}\t{}\t{}\n".format(ko_name, le1, le2, le3))
                        le1 = le1.replace('; ', ';').split(';')
                        le2 = le2.replace('; ', ';').split(';')
                        le3 = le3.replace('; ', ';').split(';')
                        for i in range(len(le1)):
                            outfile.write("{}\t{}\t{}\t{}\n".format(query, le1[i], le2[i], le3[i]))
                    else:
                        outfile2.write("-\t-\t-\t-" + "\n")
                        print "wrong ID", KO
        self.logger.info("done kegg level")
        return
        kegg_anno = self.anno_result
        self.level_anno =  self.work_dir + "/kegg_level_anno.xls"
        cmd = self.python_path + ' {} -i {} -o {} '.format(self.python_script, kegg_anno, self.work_dir)
        self.logger.info(cmd)
        command = self.add_command('kegg_level_anno', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("kegg_level_anno succeed")
        else:
            self.set_error("kegg_level_anno failed", code="31202902")
            raise Exception("kegg_level_anno failed")

    def run_kegg_stat(self):
        # cmd = self.perl_path + ' {} -q {} -p {} -o {} '. \
        #     format(self.perl_script, self.level_anno, self.option('reads_profile').prop['path'], self.output_dir)
        cmd = self.perl_path + ' {} {} {} {}'.format(self.perl_script, self.level_anno,
                                                     self.option("reads_profile").path,
                                                     self.output_dir)
        self.logger.info(cmd)
        command = self.add_command('kegg_stat', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("kegg_stat succeed")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)  # add memory limit error by guhaidong @ 20180327
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("kegg_stat failed", code="31202903")
            raise Exception("kegg_stat failed")

    def set_output(self):
        newfile = self.output_dir + "/gene_kegg_anno.xls"
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(self.work_dir + "/gene_kegg_anno_all.xls", newfile)
        self.logger.info("set_output")
        for f in os.listdir(self.output_dir):  # 删除sed的中间文件
            if f.startswith('sed'):
                fp = os.path.join(self.output_dir, f)
                os.system("rm -f " + fp)
        if len(os.listdir(self.output_dir)) == 4:
            self.logger.info("OUTPUT RIGHT")
        else:
            raise Exception("OUTPUT WRONG")
