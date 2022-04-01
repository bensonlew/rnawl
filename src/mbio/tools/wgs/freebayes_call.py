# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180404

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class FreebayesCallAgent(Agent):
    """
    author = HONGDONG
    version = 1.0
    last modify = 20180404
    SNP工具,使用Freebayes进行call snp与indel
    """
    def __init__(self, parent):
        super(FreebayesCallAgent, self).__init__(parent)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "bed_file", "type": "string"},  # 切割后的bed文件
            {"name": "bam_list", "type": "string"}  # 所有样本的bam文件列表
        ]
        self.add_option(options)
        self.step.add_steps('call_snp')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.call_snp.start()
        self.step.update()

    def step_end(self):
        self.step.call_snp.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("ref_fasta").is_set:
            raise OptionError("缺少ref_fasta参数", code="34502701")
        if not self.option("bed_file"):
            raise OptionError("缺少bed_file参数", code="34502702")
        if not self.option("bam_list"):
            raise OptionError("缺少bam_list参数", code="34502703")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".")[0]
        if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="34502704")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '100G'
        
    def end(self):
        super(FreebayesCallAgent, self).end()


class FreebayesCallTool(Tool):
    def __init__(self, config):
        super(FreebayesCallTool, self).__init__(config)
        self._version = '1.0.1'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.freebayes = 'bioinfo/WGS/freebayes/bin/'
        
    def freebayes_call(self):
        """
        注意这里加上-=是为了获取GQ值，用于后面的统计
        freebayes  -f /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.3.2.var/02.ref-config/ref.fa
        -t /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.3.2.var/step07.variant-c
        all-freebayes/1.bed -L /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.3.2.var/st
        ep07.variant-call-freebayes/bam.list
        -v /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.3.2.var/step07.variant-
        call-freebayes/1.vcf
        """
        cmd = "{}freebayes -f {} -t {} -L {} -v {}.vcf -= -X -u"\
            .format(self.freebayes, self.option("ref_fasta").prop['path'], self.option("bed_file"),
                    self.option("bam_list"), os.path.join(self.output_dir, os.path.basename(self.option('bed_file'))))
        self.logger.info(cmd)
        self.logger.info("开始进行freebayes_call")
        command = self.add_command("cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("freebayes_call完成！")
        else:
            self.set_error("freebayes_call出错！", code="34502701")
            self.set_error("freebayes_call出错！", code="34502704")

    def run(self):
        super(FreebayesCallTool, self).run()
        self.freebayes_call()
        self.end()
