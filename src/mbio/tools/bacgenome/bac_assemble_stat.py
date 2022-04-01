# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.3.26
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
class BacAssembleStatAgent(Agent):
    """
    细菌基因组扫描图组装结果评估
    """
    def __init__(self, parent):
        super(BacAssembleStatAgent, self).__init__(parent)
        options = [
            {"name": "scaf_seq", "type": "infile","format": "sequence.fasta"},  # 参考序列文件夹
            {"name": "cont_seq", "type": "infile","format": "sequence.fasta"},  # 输入列表文件
            {'name': 'sample_name', "type": "string"},  # 样本名
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("scaf_seq").is_set:
            raise OptionError("必须添加scaf_seq的序列文件！", code="31400301")
        if not self.option("cont_seq").is_set:
            raise OptionError("必须添加cont_seq的序列的文件！", code="31400302")
        if not self.option('sample_name'):
            raise OptionError('必须输入sample_name样品名称！', code="31400303")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(BacAssembleStatAgent, self).end()

class BacAssembleStatTool(Tool):
    def __init__(self, config):
        super(BacAssembleStatTool, self).__init__(config)
        self.scaf_fa = self.option("scaf_seq").prop['path']
        self.cont_fa = self.option("cont_seq").prop['path']
        self.sample_name = self.option("sample_name")
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"


    def run_asse_stat(self):
        output = self.work_dir + '/' + self.sample_name + '.summary.xls'
        cmd ='{} {}assemble_summary.pl {} 1000 {}'.format(self.perl_path,self.perl_script,self.scaf_fa,output)
        self.logger.info(cmd)
        command = self.add_command("assemble_summary", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("assemble_summary运行完成")
        else:
            self.set_error("assemble_summary运行出错!", code="31400301")


    def run_asse_seq(self):
        cmd = '{} {}bac_genome_stat.pl {} {}'.format(self.perl_path,self.perl_script,self.scaf_fa,self.cont_fa)
        self.logger.info(cmd)
        command = self.add_command("asse_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("asse_seq运行完成")
        else:
            self.set_error("asse_seq运行出错!", code="31400302")

    def run_asse_graph(self):
        for i in [1000,2000,5000]:
            output =self.sample_name + '.' + str(i)
            cmd = '{} {}seq_distribution.pl {} {} {} {} {}'.format(self.perl_path, self.perl_script, self.scaf_fa,
                                                                   self.cont_fa,i,i,output)
            self.logger.info(cmd)
            path ='asse_seq' + str(i)
            command = self.add_command(path, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("%s 运行完成" % path)
            else:
                self.set_error("%s 运行出错!" , variables=( path), code="31400303")

    def set_output(self):
        if not os.path.exists(self.output_dir + '/summary'):
            os.mkdir(self.output_dir + '/summary')
        if os.path.exists(self.output_dir + '/summary/' + self.sample_name +'_assembly_scaffold_details.xls'):
            os.remove(self.output_dir + '/summary/' + self.sample_name +'_assembly_scaffold_details.xls')
        os.link(self.work_dir + '/scaffolds.stat.xls',self.output_dir + '/summary/' + self.sample_name +'_assembly_scaffold_details.xls')
        if os.path.exists(self.output_dir + '/summary/'  + self.sample_name +'_assembly_contig_details.xls'):
            os.remove(self.output_dir + '/summary/'  + self.sample_name +'_assembly_contig_details.xls')
        os.link(self.work_dir + '/contigs.stat.xls', self.output_dir + '/summary/'  + self.sample_name +'_assembly_contig_details.xls')
        if os.path.exists(self.output_dir + '/summary/' + self.sample_name  + '_assembly_summary.xls'):
            os.remove(self.output_dir + '/summary/' + self.sample_name  + '_assembly_summary.xls')
        os.link(self.work_dir + '/' + self.sample_name  + '.summary.xls', self.output_dir + '/summary/' + self.sample_name  + '_assembly_summary.xls')
        self.add_depth()  #zouguanqing  20190327
        if not os.path.exists(self.output_dir + '/len'):
            os.mkdir(self.output_dir + '/len')
        for i in [1000,2000,5000]:
            path1 = self.output_dir + '/len/' + self.sample_name  + '.' + str(i) + '.scaffolds.len.xls'
            path2 = self.output_dir + '/len/' + self.sample_name + '.' + str(i) + '.contigs.len.xls'
            path3= self.work_dir + '/' + self.sample_name + '.' + str(i) + '.scaffolds.len.xls'
            path4 = self.work_dir + '/' + self.sample_name + '.' + str(i) + '.contigs.len.xls'
            if os.path.exists(path1):
                os.remove(path1)
            #else:
            #   os.link(path3,path1)
            os.link(path3,path1)  #guanqing.zou 20180823
            if os.path.exists(path2):
                os.remove(path2)
            #else:
            #    os.link(path4,path2)
            os.link(path4,path2)   #guanqing.zou
        self.end()

    def run(self):
        super(BacAssembleStatTool, self).run()
        self.run_asse_stat()
        self.run_asse_seq()
        self.run_asse_graph()
        self.set_output()

    def add_depth(self):
        clean_data_summary = os.path.join(self.work_dir, '../../BacgenomeQc/CleanStat/output/data_QC/{}_Illumina_statistics.xls'.format(self.option('sample_name')))
        if os.path.exists(clean_data_summary):
            with open(clean_data_summary) as f:
                line = f.readline()
                line = f.readline()
                raw_base = float(line.split('\t')[4])
            summary_file = self.output_dir + '/summary/' + self.sample_name  + '_assembly_summary.xls'
            with open(summary_file) as f1, open(summary_file+'_1','w') as fw:
                line = f1.readline()
                fw.write(line.strip()+'\tDepth\n')
                line = f1.readline()
                genome_size = float(line.split('\t')[1])
                fw.write(line.strip()+'\t'+str(round(raw_base/genome_size,2))+'\n')
            os.remove(summary_file)
            os.rename(summary_file+'_1', summary_file)

