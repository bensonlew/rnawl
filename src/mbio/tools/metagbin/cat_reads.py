# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
from biocluster.core.exceptions import OptionError
import subprocess


class CatReadsAgent(Agent):
    """
    catreads:通过cat命令将多个序列合并为一个
    version 1.0
    author: guhaidong
    last_modify: 2017.09.06
    """

    def __init__(self, parent):
        super(CatReadsAgent, self).__init__(parent)
        options = [
            {"name": "map_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 输入合并序列fq路径
            {"name": "fastq1", "type": "outfile", "format": "sequence.fastq"},  # 输入文件,read1
            {"name": "fastq2", "type": "outfile", "format": "sequence.fastq"},  # 输入文件,read2
            {"name": "fastqs", "type": "outfile", "format": "sequence.fastq"},  # 输入文件,reads
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("map_dir").is_set:
            raise OptionError("请传入fq序列路径", code="31400701")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'


class CatReadsTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(CatReadsTool, self).__init__(config)
        # self.sh_path = '/bioinfo/metaGenomic/scripts/'
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'

    def cat_reads(self):
        """
        合并序列
        :return:
        """
        samp_list, file_dict = self.get_info()
        self.logger.info('{}'.format(file_dict))
        have_s_fq = False
        cmd1 = self.sh_path + 'cat_seq.sh'
        cmd2 = self.sh_path + 'cat_seq.sh'
        cmd3 = self.sh_path + 'cat_seq.sh'
        for sample in samp_list:
            if 'l' in file_dict[sample].keys():
                cmd1 += ' ' + self.option('map_dir').prop['path'] + '/' + file_dict[sample]['l']
                self.logger.info('{}'.format(file_dict[sample]['l']))
            if 'r' in file_dict[sample].keys():
                cmd2 += ' ' + self.option('map_dir').prop['path'] + '/' + file_dict[sample]['r']
                self.logger.info('{}'.format(file_dict[sample]['r']))
            if 's' in file_dict[sample].keys():
                cmd3 += ' ' + self.option('map_dir').prop['path'] + '/' + file_dict[sample]['s']
                self.logger.info('{}'.format(file_dict[sample]['s']))
                have_s_fq = True
        cmd1 += ' ' + self.output_dir + '/PE.l.fq'
        cmd2 += ' ' + self.output_dir + '/PE.r.fq'
        self.logger.info('运行cat_reads，将reads进行合并')
        command1 = self.add_command("cat_reads_cmd1", cmd1).run()
        command2 = self.add_command("cat_reads_cmd2", cmd2).run()
        self.wait(command1)
        self.wait(command2)
        if command1.return_code == 0:
            self.logger.info("左端fq合并完成")
        else:
            self.set_error("左端fq合并失败！", code="31400701")
        if command2.return_code == 0:
            self.logger.info("右端fq合并完成")
        else:
            self.set_error("右端fq合并失败！", code="31400702")
        if have_s_fq:
            cmd3 += ' ' + self.output_dir + '/PE.s.fq'
            command3 = self.add_command("cat_reads_cmd3", cmd3).run()
            self.wait(command3)
            if command3.return_code == 0:
                self.logger.info("single端fq合并完成")
            else:
                self.set_error("single端fq合并失败！", code="31400703")

    def get_info(self):
        """
        从list表中获取文件对应PSE信息
        :return:
        samp_list = [sample1, sample2, sample3 ...]
        file_dict[samp]['l|r|s'] = file_name
        """
        file_dict = dict()
        samp_list = list()
        info_list = os.path.join(self.option('map_dir').prop['path'], "list.txt")
        with open(info_list, "r") as files:
            for line in files:
                line = line.strip().split()
                if line[1] not in file_dict.keys():
                    file_dict[line[1]] = {line[2]: line[0]}
                else:
                    file_dict[line[1]][line[2]] = line[0]
                if line[1] not in samp_list:
                    samp_list.append(line[1])
        return samp_list, file_dict

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.logger.info("设置结果目录")
        self.option('fastq1').set_path(self.output_dir + '/PE.l.fq')
        self.option('fastq2').set_path(self.output_dir + '/PE.r.fq')
        self.option('fastqs').set_path(self.output_dir + '/PE.s.fq')
        self.logger.info("设置结果目录成功")

    def run(self):
        super(CatReadsTool, self).run()
        self.cat_reads()
        self.set_output()
        self.end()
