# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# version 1.0
# last_modify: 2018.06.013

import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna.trans_step import step_count


class GffFaaFnnFormatAgent(Agent):
    """
    gff文件提取maker预测信息，重命名基因名称， faa fnn 文件基因名称随之改变
    """

    def __init__(self, parent):
        super(GffFaaFnnFormatAgent, self).__init__(parent)
        options = [
            {"name": "gff", "type": "string", "default": ""},
            {"name": "faa", "type": "infile", "format": "sequence.fasta"},
            {"name": "ffn","type":"infile","format": "sequence.fasta"},
            {"name": "pre_name","type":"string","default":"gene"},
            {"name": "genome","type":"infile","format": "sequence.fasta"},
            {"name": "sample", "type": "string","default":"sample"}
        ]

        self.add_option(options)

    def check_options(self):

        if not self.option("gff"):
            raise OptionError("必须设置参gff", code="32101101")
        if not self.option("faa").is_set:
            raise OptionError("必须设置参数faa", code="32101102")
        if not self.option("ffn").is_set:
            raise OptionError("必须设置参数ffn", code="32101103")
        if not self.option("genome").is_set:
            raise OptionError("必须设置参数genome", code="32101104")

    def set_resource(self):
        self._cpu = 1
        self._memory = '1G'

    def end(self):
        super(GffFaaFnnFormatAgent, self).end()


class GffFaaFnnFormatTool(Tool):
    def __init__(self, config):
        super(GffFaaFnnFormatTool, self).__init__(config)
        self.maker_format = os.path.join(self.config.PACKAGE_DIR, "fungi_genome/gff_maker_format.pl")
        self.gff_add =os.path.join(self.config.PACKAGE_DIR, "fungi_genome/gff_add_info.py")
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.python_path = "/program/Python/bin/python"


    def run_cmd1(self):
        cmd1 = "{} {} -gff {} -transcripts {} -proteins {} -genename {}".format(self.perl_path, self.maker_format, self.option("gff"), self.option("ffn").prop['path'], self.option("faa").prop['path'], self.option("pre_name"))
        self.logger.info("开始运行{}".format(cmd1))
        command = self.add_command("run_cmd1", cmd1).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd1))
            self.pos_dic = os.path.join(self.work_dir, "pos.dic")
            self.proten_info = os.path.join(self.work_dir, "protein.info")
            self.gene_info = os.path.join(self.work_dir, "gene.info")
            self.gff_format = os.path.join(self.work_dir,"gff.format")
            self.faa = os.path.join(self.work_dir,'proteins.format')
            self.ffn = os.path.join(self.work_dir,'transcripts.format')
        else:
            self.set_error("%s运行出错!", variables=(cmd1), code="32101101")

    def run_cmd2(self):
        cmd2 = "{} {} {} {} {} {} {} {} {}".format(self.python_path, self.gff_add, self.option("genome").prop["path"], self.proten_info, self.gene_info, self.gff_format, self.ffn, self.faa, self.option('sample'))
        self.logger.info("开始运行{}".format(cmd2))
        command = self.add_command("run_cmd2", cmd2).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd2))
            self.new_gff = os.path.join(self.work_dir,"new.gff")
            self.new_faa = os.path.join(self.work_dir,'new.format.faa')
            self.new_ffn = os.path.join(self.work_dir,"new.format.ffn")
        else:
            self.set_error("%s运行出错!", variables=(cmd2), code="32101102")

    def run_get_information(self):
        perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        port_seq = self.work_dir + '/new.format.faa'
        #nul_seq = self.work_dir + '/new.format.ffn'
        gene_seq = self.work_dir + '/' + self.option('sample')
        genome_fasta = self.option("genome").prop['path']
        cmd = '{} {}trim_gene_info.pl {} {} {} {} {}'.format(self.perl_path, perl_script, genome_fasta,
                                                             gene_seq, port_seq, "ORF",
                                                             self.output_dir)
        command = self.add_command("get_information", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("统计预测结果运行完成")
        else:
            self.set_error("统计预测结果运行出错!", code="32101103")

        self.logger.info("开始计算基因分布")
        #self.logger.info(step_count(self.option('ffn').prop['path'], self.work_dir + "/fnn_stat.xls", 11, 200,
        #                                self.output_dir + "/{}_CDS_length.xls".format(self.option('sample'))))
        #step_count(self.option('ffn').prop['path'], self.work_dir + "/fnn_stat.xls", 11, 200,
        #               self.output_dir + "/{}_CDS_length.xls".format(self.option('sample')))

        self.logger.info(step_count(gene_seq, self.work_dir + "/fnn_stat.xls", 11, 200,
                                       self.output_dir + "/{}_CDS_length.xls".format(self.option('sample'))))
        step_count(gene_seq, self.work_dir + "/fnn_stat.xls", 11, 200,
                      self.output_dir + "/{}_CDS_length.xls".format(self.option('sample')))

    def set_output(self):
        new_paths ="{1}/{0}_CDS.gff,{1}/{0}_CDS.faa,{1}/{0}_CDS.fnn".format(self.option("sample"),self.output_dir).split(',')
        for new_path in new_paths:
            if os.path.exists(new_path):
                os.remove(new_path)
        os.link(self.new_gff,new_paths[0])
        os.link(self.new_faa,new_paths[1])
        os.link(self.new_ffn,new_paths[2])
        cds_xls = self.output_dir+"/{}_CDS_statistics.xls".format(self.option("sample"))
        #len_xls = "{0}/{1}_CDS_length.xls".format(self.output_dir,self.option("sample"))
        os.system("mv {} {}".format(self.output_dir+"/gene_statistics.xls", cds_xls))
        #os.system('rm {0}/*proteins.fasta {0}/*transcripts.fasta'.format(self.output_dir))
        os.system('rm {0}/new.format.faa {0}/{1}'.format(self.output_dir, self.option("sample")))
        os.system('head -n 13 {0}/{1}_CDS_length.xls > {0}/{1}_CDS_length.xls_old'.format(self.output_dir,self.option("sample")))
        os.system('mv {0}/{1}_CDS_length.xls_old {0}/{1}_CDS_length.xls'.format(self.output_dir,self.option("sample")))
        #os.system("sed -i 's/new.format/{0}/' {1} ".format(self.option("sample"),len_xls))
        #os.system("sed -i 's/new.format.faa/{0}/' {1} ".format(self.option("sample"),cds_xls))


    def run(self):
        super(GffFaaFnnFormatTool, self).run()
        self.run_cmd1()
        self.run_cmd2()
        self.run_get_information()
        self.set_output()
        self.end()
