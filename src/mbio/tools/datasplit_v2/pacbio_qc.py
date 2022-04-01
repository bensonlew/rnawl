# -*_ coding: utf-8 -*-
# __author__ = 'XueQinwen'
# last_modified: 20210907
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os

class PacbioQcAgent(Agent):
    def __init__(self,parent):
        super(PacbioQcAgent,self).__init__(parent)
        options = [
            {"name":"input_bam","type":"infile","format":"align.bwa.bam"},
            {"name":"sn_name","type":"string"},
            {"name":"primer_type","type":"string"}
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_option(self):
        """
        参数检查
        """
        if not self.option("input_bam"):
            raise OptionError("没有输入文件")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(PacbioQcAgent, self).end()

class PacbioQcTool(Tool):
    def __init__(self, config):
        super(PacbioQcTool, self).__init__(config)
        self.perl = "program/perl-5.24.0/bin/perl"
        self.python3 = "program/Python35/bin"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.python3, LD_LIBRARY_PATH=self.lib)
        self.qc_script = self.config.PACKAGE_DIR + '/datasplit/trim_fqSeq_pacbio.pl'
        self.bam2fastq = "program/SmrtLink/smrtlink/smrtcmds/bin/bam2fastq"
        self.fq_rename = self.config.PACKAGE_DIR + "/datasplit/fastqRename.py"


    def run(self):
        super(PacbioQcTool, self).run()
        self.run_bam2fastq()
        self.run_qc()
        self.set_output()
        self.end()


    def run_bam2fastq(self):
        bamfile = self.option("input_bam").prop['path']
        self.sample_name = os.path.basename(bamfile).split(".")[1]
        cmd = "{} -o {}.ccs {}".format(self.bam2fastq,self.sample_name, bamfile)
        command = self.add_command("bam2fastq",cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fastq生成完成")
        else:
            self.set_error("fastq出现错误")
        
    def run_qc(self):
        self.fastqfile = os.path.join(self.work_dir,"{}.ccs.fastq.gz".format(self.sample_name))
        self.trim_fastqfile = os.path.join(self.work_dir,"{}_trim.fastq".format(self.sample_name))
        self.out_fastqfile = os.path.join(self.output_dir,"{}_value.fastq".format(self.sample_name))
        self.logger.info("正在给样本{}进行长度指控，引物{}".format(self.sample_name,self.option("primer_type"))) 
        cmd = "{} {} -i {} -o {}".format(self.perl,self.qc_script,self.fastqfile, self.trim_fastqfile)
        if self.option("primer_type")=="27F_1492R" or "20F_1492R":
            cmd += " -m 1000 -x 1800"
        elif self.option("primer_type") == "ITS1F_ITS4R":
            cmd += " -m 300 -x 900"
        else:
            self.set_error("引物不在三代拆分QC范围内")
        self.logger.info(cmd)
        command = self.add_command("qc_smrt", cmd).run()
        self.wait(command)
        if command.return_code  == 0:
            self.logger.info("质控完成")
        else:
            self.set_error("质控发生错误")
        cmd = "{}/python {} -i {} -o {} -n {}".format(self.python3, self.fq_rename,
                  self.trim_fastqfile, self.out_fastqfile, self.option("sn_name"))
        command = self.add_command("fastq_rename", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("fastqRenamefq重命名运行成功")
        else:
            self.set_error("fastqRenamefq重命名运行失败")
        os.system("gzip {}".format(self.out_fastqfile))

    def set_output(self):
        try:
            os.link(self.fastqfile,
                        os.path.join(self.output_dir,"{}.ccs.fastq.gz".format(self.sample_name)))
        except Exception as e:
            self.set_error("设置结果目录失败{}".format(e))


    
    

                 


    