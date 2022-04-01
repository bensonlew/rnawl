#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re


class IslandIslanderAgent(Agent):
    """
    用于扫描图的gbk文件生成
    version 1.0
    author: gaohao
    last_modify: 2018.04.02
    """

    def __init__(self, parent):
        super(IslandIslanderAgent, self).__init__(parent)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fasta_dir"},  #
            {"name": "out", "type": "outfile", "format": "bacgenome.island"},
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('fa_dir').is_set:
            raise OptionError("请设置基因组学列文件夹！", code="31402301")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./fastq_stat.xls", "xls", "fastq信息统计表"]
        ])
        super(IslandIslanderAgent, self).end()


class IslandIslanderTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(IslandIslanderTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.islander = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/Islander_software/Islander.pl"

    def run_islander(self):
        filelist =os.listdir(self.option('fa_dir').prop['path'])
        for file in filelist:
            de =file.split('.')[0]
            if os.path.exists(self.work_dir + '/' + de + '.fna'):
                os.remove(self.work_dir + '/' + de + '.fna')
            os.link(self.option('fa_dir').prop['path']+ '/' + file,self.work_dir + '/' + de + '.fna')
            fa=self.work_dir + '/' + de
            cmd = "{} {} {} --verbose --translate --trna --annotate --reisland --table 11 --nocheck".format(
                self.perl_path, self.islander,fa)
            self.logger.info(cmd)
            path = str(de.lower()) + '_run_islander'
            self.logger.info("开始运行%s "%path)
            command = self.add_command(path, cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行%s完成" %path)
                dd=de + '*'
                os.system('rm %s'%dd)
            else:
                self.set_error("运行%s运行出错!" , variables=(path), code="31402301")



    def set_output(self):
        path =self.work_dir + '/all_final_arrayed_islands.gff'
        self.option('out').set_path(path)

    def run(self):
        """
        运行
        """
        super(IslandIslanderTool, self).run()
        self.run_islander()
        self.set_output()
        self.end()
