#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
import time
import subprocess


class BacGbkAgent(Agent):
    """
    用于扫描图的gbk文件生成
    version 1.0
    author: qingchen.zhang
    fix tool by bacgenome bac_gbk.py
    """

    def __init__(self, parent):
        super(BacGbkAgent, self).__init__(parent)
        options = [
            {"name": "gen_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "trna_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "pro_fa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "anno", "type": "infile", "format": "sequence.profile_table"},  #
            {"name": "sample_name", "type": "string"},
            {"name": "analysis", "type": "string", "default": "uncomplete"},  ###流程分析模式complete，uncomplete
            {"name": "gbk_dir", "type": "outfile", "format": "gene_structure.gbk_dir"},  # 输出文件
            {"name": "all_gbk", "type": "outfile", "format": "gene_structure.gbk"},  # 输出文件
        ]
        self.add_option(options)
        self.list =[]

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('gen_gff').is_set:
            raise OptionError("请设置基因组基因预测基因gff文件！")
        if not self.option('rrna_gff').is_set:
            raise OptionError("请设置基因组rRNA预测基因gff文件！")
        if not self.option('trna_gff').is_set:
            raise OptionError("请设置基因组tRNA预测基因gff文件！")
        if not self.option('pro_fa').is_set:
            raise OptionError("请设置基因组基因蛋白序列文件！")
        if not self.option('genome_fa').is_set:
            raise OptionError("请设置基因组组装序列文件！")
        if not self.option('anno').is_set:
            raise OptionError("请设置基因组基因注释汇总表！")
        if not self.option('analysis'):
            raise OptionError("请提供分析流程类型！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./fastq_stat.xls", "xls", "fastq信息统计表"]
        ])
        super(BacGbkAgent, self).end()


class BacGbkTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(BacGbkTool, self).__init__(config)
        self.genegff =self.option('gen_gff').prop['path']
        self.rrnagff = self.option('rrna_gff').prop['path']
        self.trnagff = self.option('trna_gff').prop['path']
        self.genomefa = self.option('genome_fa').prop['path']
        self.pro = self.option('pro_fa').prop['path']
        self.anno = self.option('anno').prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.perl_script2 = self.config.PACKAGE_DIR + "/fungi_genome/"
        self.list = []
        self.python = self.config.SOFTWARE_DIR + '/program/Python/bin/python'
        self.python_scrit = self.config.PACKAGE_DIR + "/bacgenome/gbk_rm_ty_tbl.py"
        self.tbl2asn = "bioinfo/Genomic/Sofware/tbl2asn/linux64.tbl2asn"
        self.template_sbt = self.config.PACKAGE_DIR + "/bacgenome/template.sbt"
        self.all_gbk = self.work_dir + '/' +self.option('sample_name') + '.gbk'


    def run_tbl(self):
        """
        生成tbl文件
        """
        cmd = "{} {}tbl_generation.pl uncomplete {} {} {} {} {} {} ".format(self.perl_path, self.perl_script2,self.genegff,self.trnagff,self.rrnagff,self.pro ,self.genomefa,self.anno)
        self.logger.info(cmd)
        self.logger.info("开始运行run_tbl")
        command = self.add_command("run_tbl", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_tbl完成")
        else:
            self.set_error("运行run_tbl运行出错!")

    def run_gbk(self):
        """
        根据tbl、fsa文件生成gbk文件
        """
        self.get_list()
        for file in self.list:
            fsa =self.work_dir + '/' + file + '/' + file + '.fsa'
            tbl = self.work_dir + '/' + file + '/' + file + '.tbl'
            ##guanqing.zou 20180904
            with open(tbl) as f:
                lines = f.readlines()
            if len(lines)==1:
                cmd = "{} -t {} -i {}  -V vb".format(self.tbl2asn, self.template_sbt, fsa)
            else:
                cmd = "{} -t {} -i {} -k {} -V vb".format(self.tbl2asn, self.template_sbt, fsa, tbl)

            self.logger.info(cmd)
            path = file.lower() + '_gbk_run'
            self.logger.info("开始运行%s" %path)
            try:
                subprocess.check_output(self.config.SOFTWARE_DIR + '/'+cmd, shell=True)
                self.logger.info("运行%s完成" % path)
            except Exception as e:
                self.logger.info(e)
            os.system("sed  -i '/^LOCUS/s/circular//' %s"%(self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system("sed  -i '/^LOCUS/s/linear//' %s"%(self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system('sed \'s/ACCESSION/ACCESSION   %s/g\' %s -i' % (file, self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system('sed \':a;N;$!ba;s/-\\n[ ]*/M/g\' %s -i' % ( self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system('%s %s %s %s %s' % (self.python, self.python_scrit, self.work_dir + '/' + file + '/' + file + '.tbl', self.work_dir + '/' + file + '/' + file + '.gbf',
                                              self.work_dir + '/' + file + '/' + file + '.gbk'))
            os.system('cat %s >> %s' % (self.work_dir + '/' + file + '/' + file + '.gbk', self.all_gbk))
            self.logger.info("完成运行%s" %path)

    def run_bio_gbk_un(self):
        """
        标准化gbk文件
        """
        if os.path.exists(self.work_dir + 'all.antismash.gbk'):
            os.remove(self.work_dir + 'all.antismash.gbk')
        cmd = "{} {}GBK_generation_un.pl {} {} {} {} {} {} ".format(self.perl_path, self.perl_script,
                                                                            self.genegff, self.trnagff, self.rrnagff,
                                                                            self.pro, self.genomefa, self.anno)
        self.logger.info(cmd)
        self.logger.info("开始运行run_bio_gbk")
        command = self.add_command("run_bio_gbk", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_bio_gbk完成")
        else:
            self.set_error("运行run_bio_gbk运行出错!")

    def get_list(self):
        with open(self.genomefa,'r') as f:
            lines=f.readlines()
            for line in lines:
                if re.search(r'^>',line):
                    line = line.replace('>', '')
                    tmp = line.rstrip('\n').split(' ')[0]
                    print(tmp)
                    self.list.append(tmp)


    def set_output(self):
        """
        设置结果文件目录
        """
        self.logger.info("set output")
        if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '.gbk'):
            os.remove(self.output_dir + '/gbk/Scaffold/' + self.option('sample_name') + '.gbk')
        os.link(self.all_gbk, self.output_dir + '/' + self.option('sample_name') + '.gbk')

        ##链接gff和ptt的结果
        if os.path.exists(self.output_dir + '/'+  self.option('sample_name') + '.gff'):
            os.remove(self.output_dir + '/'+  self.option('sample_name') + '.gff')
        os.link(self.work_dir + '/gff/out.gff',self.output_dir + '/'+  self.option('sample_name') + '.gff')

    def run_gff_ptt(self):
        """
        生成合并的gff文件
        """
        if self.option("analysis") == 'uncomplete':
            gff_src = "generation_gff.pl"
            cmd = "{} {}{} {} {} {} {} {}".format(self.perl_path, self.perl_script2, gff_src,
                                                                self.genegff, self.trnagff, self.rrnagff,self.anno, "gff")
            self.logger.info(cmd)
            self.logger.info("开始生成gff")
            command = self.add_command("run_bio_gff", cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("gff生成完成")
            else:
                self.set_error("gff生成出错!")

    def run(self):
        """
        运行
        """
        super(BacGbkTool, self).run()
        if self.option('analysis') in ['uncomplete']:
            self.run_tbl()
            self.run_gbk()
            # self.run_bio_gbk_un()
            self.run_gff_ptt()
            time.sleep(5)
            self.set_output()
            self.end()

