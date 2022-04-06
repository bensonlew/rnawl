# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.05.25
import os
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import shutil

class CompleteAsseStatAgent(Agent):
    """
    细菌基因组完成图组装结果评估
    """
    def __init__(self, parent):
        super(CompleteAsseStatAgent, self).__init__(parent)
        options = [
            {"name": "fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "table", "type": "infile", "format": "sequence.profile_table"}, #质粒注释文件
            {"name": "output", "type": "outfile", "format": "meta.otu.otu_table"},##输出质粒基因组与基因前缀对应关系
            {"name": "chr", "type": "outfile", "format": "sequence.fasta"},
            {"name": "pla", "type": "outfile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"},
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fa").is_set:
            raise OptionError("必须添加fa的序列文件！")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(CompleteAsseStatAgent, self).end()

class CompleteAsseStatTool(Tool):
    def __init__(self, config):
        super(CompleteAsseStatTool, self).__init__(config)
        self.fa = self.option("fa").prop['path']
        self.table = self.option("table").prop['path']
        self.ample = self.option('sample_name')
        self.perl_path = "/miniconda2/bin/python"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"

    def run_asse_stat(self):
        cmd ='{} {}complete_assemble.py -f {} -d {} -p {} -dir {}'.format(self.perl_path,self.perl_script,self.fa,self.table, self.ample, self.output_dir)
        self.logger.info(cmd)
        command = self.add_command("run_asse_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_asse_stat运行完成")
        else:
            self.set_error("run_asse_stat运行出错!", code="31401001")

    def set_output(self):
        if os.path.getsize(self.work_dir + '/plasmid.type.xls') >0:
           self.option('output').set_path(self.work_dir + '/plasmid.type.xls')
        if os.path.exists(self.work_dir + '/plasmid.fasta') and os.path.getsize(self.work_dir + '/plasmid.fasta') != 0:
            self.option('pla').set_path(self.work_dir + '/plasmid.fasta')
        if os.path.exists(self.work_dir + '/chromosome.fasta'):
            self.option('chr').set_path(self.work_dir + '/chromosome.fasta')
        if os.path.exists(self.output_dir + '/' + self.ample + '_assembly_details.xls'):
            os.remove(self.output_dir + '/' + self.ample + '_assembly_details.xls')
        os.link(self.work_dir + '/' + self.ample + '_assembly_details.xls',self.output_dir + '/' + self.ample + '_assembly_details.xls')
        if os.path.exists(self.output_dir + '/' + self.ample + '_assembly_summary.xls'):
            os.remove(self.output_dir + '/' + self.ample + '_assembly_summary.xls')
        os.link(self.work_dir + '/' + self.ample + '_assembly_summary.xls', self.output_dir + '/' + self.ample + '_assembly_summary.xls')
        self.add_depth()  #zouguanqing 20190327

    def run(self):
        super(CompleteAsseStatTool, self).run()
        self.run_asse_stat()
        self.set_output()
        self.end()

    def linkdir(self, dirpath, dirname):
        """
		link一个文件夹下的所有文件到本module的output目录
		:param dirpath: 传入文件夹路径
		:param dirname: 新的文件夹名称
		:return:
		"""
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))

    def add_depth(self):
        clean_data_summary = os.path.join(self.work_dir, '../BacgenomeQc/PacbioClean/output/{}.PacBio_statistics.xls'.format(self.option('sample_name')))
        if os.path.exists(clean_data_summary):
            with open(clean_data_summary) as f:
                line = f.readline()
                line = f.readline()
                base = float(line.split('\t')[1])
            summary_file = self.output_dir + '/' + self.ample + '_assembly_summary.xls'
            with open(summary_file) as f1, open(summary_file+'_1','w') as fw:
                line = f1.readline()
                fw.write(line.strip()+'\tDepth\n')
                line = f1.readline()
                genome_size = float(line.split('\t')[2])
                fw.write(line.strip()+'\t'+str(round(base/genome_size,2))+'\n')
            os.remove(summary_file)
            os.rename(summary_file+'_1', summary_file)