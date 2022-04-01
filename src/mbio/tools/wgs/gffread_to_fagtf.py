# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180813

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class GffreadToFagtfAgent(Agent):
    """
    gffread转ref_gff,fa文件为*.gtf,*.protein.fa,*.mRNA.fa,*.cds.fa
    lasted modified by hongdong@20190408 增加ref.new.fa
    """
    def __init__(self, parent):
        super(GffreadToFagtfAgent, self).__init__(parent)
        options = [
            {"name": "fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "gff_path", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('GffreadToFagtf')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.GffreadToFagtf.start()
        self.step.update()

    def step_end(self):
        self.step.GffreadToFagtf.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fa"):
            raise OptionError("请设置fa")
        if not self.option("gff_path"):
            raise OptionError("请设置gff_path参数")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(GffreadToFagtfAgent, self).end()


class GffreadToFagtfTool(Tool):
    def __init__(self, config):
        super(GffreadToFagtfTool, self).__init__(config)
        # self.gffread_path = "/bioinfo/rna/cufflinks-2.2.1/gffread"   # 需要下载最新版。崔青美20180813
        self.gffread_path = '/bioinfo/WGS/gffread/gffread-master/gffread'   # cufflinks中gffread存在bug，换软件
        self.gene_path = self.config.PACKAGE_DIR + "/wgs/enrichFastaTitle.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '

    def gffreadtofagtf(self):
        """
        cd /mnt/ilustre/users/qingmei.cui/newmdt/Genome;
        gffread -T -o out.gtf -x out.cds.fa -w out.mRNA.fa -y out.protein.fa
        ./Arabidopsis_thaliana/NCBI/TAIR10_2011-05-11_Columbia/01.newref/ref.gff
        -g ./Arabidopsis_thaliana/NCBI/TAIR10_2011-05-11_Columbia/01.newref/ref.fa
        gffread -J -W -w ref.mRNA.fa -x ref.cds.fa -y 01.newref/ref.pep.fa ref.gff -g ref.fa
        :return:
        """
        cmd = "{} -J -W -T {} -x {}".format(self.gffread_path, self.option("gff_path"), self.output_dir + "/ref.cds.fa")
        cmd += " -w {}".format(self.output_dir + "/ref.mRNA.fa")
        cmd += " -y {}".format(self.output_dir + "/ref.protein.fa")
        cmd += " -g {}".format(self.option('fa').prop['path'])
        cmd += " -o {}".format(self.work_dir + "/ref.gtf")
        self.logger.info(cmd)
        self.logger.info("开始进行GffreadToFagtf")
        command = self.add_command("gffreadtofagtf", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("GffreadToFagtf完成！")
        else:
            self.set_error("GffreadToFagtf出错！")

    def enrichfastatitle(self):
        """
        perl enrichFastaTitle.pl -i ref.mRNA.fa -g ref.gff -o 01.newref/ref.new.mRNA.fa
        :return:
        """
        cmd = "{}{} -i {} -g {} -o {}" \
            .format(self.perl_path, self.gene_path, self.output_dir + "/ref.mRNA.fa", self.option("gff_path"),
                    self.output_dir + "/ref.new.mRNA.fa")
        self.logger.info(cmd)
        self.logger.info("开始进行Getgenefasta")
        command1 = self.add_command("enrichfastatitle", cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("enrichfastatitle完成！")
        else:
            self.set_error("enrichfastatitle出错！")

    def check_mRNA_faid(self):
        """
        检查fa中id的格式是否满足条件， "transcriptid:genename|--|--|--|--:seqid：start：end"
        :return:
        """
        with open(self.output_dir + "/ref.new.mRNA.fa", 'r') as r:
            for line in r:
                if re.match('>.*', line):
                    if len(line.split('>')[-1].split(':')) != 5:
                        self.set_error("ref.new.mRNA.fa文件中的id格式不正确，请检查！")

    def run(self):
        super(GffreadToFagtfTool, self).run()
        self.gffreadtofagtf()
        self.enrichfastatitle()
        self.end()
