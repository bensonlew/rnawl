#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re


class IslandDimobAgent(Agent):
    """
    用于扫描图的gbk文件生成
    version 1.0
    author: gaohao
    last_modify: 2018.04.02
    """

    def __init__(self, parent):
        super(IslandDimobAgent, self).__init__(parent)
        options = [
            {"name": "gbk_dir", "type": "infile", "format": "gene_structure.gbk_dir"},  #
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('gbk_dir').is_set:
            raise OptionError("请设置基因组基因gbk文件夹不存在！", code="31402201")

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
        super(IslandDimobAgent, self).end()


class IslandDimobTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(IslandDimobTool, self).__init__(config)
        self.gbk =self.option('gbk_dir').prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.dimob = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/islandpath_dimob/Dimob.pl"


    def run_dimob(self):
        des = 'location' + '\t' + 'sofware' + '\t' + 'GI_num' + '\t' + 'start' + '\t' + 'end'
        os.system("echo \'%s\' >all.dimob.xls" % des)
        files =os.listdir(self.gbk)
        for file in files:
            spf = file.split('.')
            de = file.split('.gbk')[0]##gaohao 20191010
            #de = spf[0]
            self.output = self.work_dir + '/' + de + '.IG.xls'
            num =self.get_len(self.gbk + '/' + file)
            if num in ['true']:
                newfile ='_'.join(spf[:-1]) + '.' + spf[-1]
                new_path = self.work_dir + '/' + newfile
                if os.path.exists(new_path):
                    os.remove(new_path)
                os.link(self.gbk + '/' + file, new_path)
                os.system("sed 's/*/-/g' -i {}".format(new_path))
                cmd = "{} {} {} {}".format(self.perl_path, self.dimob, new_path, self.output)
                self.logger.info(cmd)
                path = de.lower() + '_run_dimob'
                self.logger.info("开始运行%s " % path)
                command = self.add_command(str(path), cmd)
                command.run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("运行%s完成" % path)
                    if os.path.getsize(self.output) != 0:
                        path = self.work_dir + '/' + de + '.IG.txt'
                        with open(self.output, "r") as f, open(path, 'w') as d:
                            lines = f.readlines()
                            for line in lines:
                                d.write(de + '\t' + 'IslandPath-DIMOB' + '\t' + line)
                        f.close()
                        d.close()
                        os.system("cat %s >>all.dimob.xls" % path)
                else:
                    self.set_error("运行%s运行出错!" , variables=( path), code="31402201")


    def set_output(self):
        path =self.work_dir + '/all.dimob.xls'
        if os.path.exists(path):
            self.option('out').set_path(path)


    def run(self):
        """
        运行
        """
        super(IslandDimobTool, self).run()
        self.run_dimob()
        self.set_output()
        self.end()

    def get_len(self,file):
        num =''
        with open(file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line =line.rstrip('\r\n')
                if re.search('     CDS',line):
                    num = 'true'
                    break
                else:
                    num = 'false'
        return num

