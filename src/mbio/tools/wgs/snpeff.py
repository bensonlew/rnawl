# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180411

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import random
import os
import re


class SnpeffAgent(Agent):
    """
    snp/indel比较分析接口中用于样本比较分析
    """
    def __init__(self, parent):
        super(SnpeffAgent, self).__init__(parent)
        options = [
            {"name": "ref_name", "type": "string"},
            {"name": "snpEff_config", "type": "string"},   # snpEff.config
            {"name": "filter_recode_vcf", "type": "string"},   # pop.snp.filter.recode.vcf/pop.indel.filter.recode.vcf
            {"name": "types", "type": "string", "default": "snp"}
        ]
        self.add_option(options)
        self.step.add_steps('snpeff')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.snpeff.start()
        self.step.update()

    def step_end(self):
        self.step.snpeff.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("snpEff_config"):
            raise OptionError("缺少snpEff_config参数", code="34505701")
        if not self.option("filter_recode_vcf"):
            raise OptionError("缺少filter_recode_vcf参数", code="34505702")
        if not self.option("types"):
            raise OptionError("缺少variation_type参数", code="34505703")
        else:
            if self.option("types") not in ["snp", "indel"]:
                raise OptionError("分析类型%s不合法！",variables=(self.option("types")), code="34505704")
        if not self.option("ref_name"):
            raise OptionError("缺少ref_fasta参数", code="34505705")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '55G'
        
    def end(self):
        super(SnpeffAgent, self).end()


class SnpeffTool(Tool):
    def __init__(self, config):
        super(SnpeffTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.snpeff_path = self.config.SOFTWARE_DIR + '/bioinfo/annotation/snpeff/'
        self.snpEff_config_base = self.config.SOFTWARE_DIR
        
    def snp_eff(self):
        """
        java -Xmx50G -jar snpEff.jar -v ref -csvStats snp.anno.csv -c snpEff.config
        pop.snp.filter.recode.vcf > snp.anno.primary.vcf
        && java -Xmx50G -jar snpEff.jar -v ref -csvStats
        indel.anno.csv -c snpEff.config pop.indel.filter.recode.vcf > indel.anno.primary.vcf
        :return:
        """
        self.make_new_snpEff_config()
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "script_{}.sh".format(now_time)
        cmd1 = "java -Xmx50G -jar {}snpEff.jar -v {} " \
            .format(self.snpeff_path, self.option("ref_name"))
        if self.option("types") == "snp":
            cmd1 += "-csvStats {} -c {} {} > {}"\
                .format(os.path.join(self.output_dir, "snp.anno.csv"), self.work_dir + "/snpEff.config",
                        self.option("filter_recode_vcf"), os.path.join(self.output_dir, "snp.anno.primary.vcf"))
        elif self.option("types") == "sv":
            cmd1 += "-csvStats {} -c {} {} > {}" \
                .format(os.path.join(self.output_dir, "pop.sv.anno.csv"), self.work_dir + "/snpEff.config",
                        self.option("filter_recode_vcf"), os.path.join(self.output_dir, "pop.sv.anno.vcf"))
        else:
            cmd1 += "-csvStats {} -c {} {} > {}" \
                .format(os.path.join(self.output_dir, "indel.anno.csv"), self.work_dir + "/snpEff.config",
                        self.option("filter_recode_vcf"), os.path.join(self.output_dir, "indel.anno.primary.vcf"))
        self.logger.info(cmd1)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd1)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！",variables=(file_path), code="34505701")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始进行snp_filter")
        command1 = self.add_command("snp_eff", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("snp_eff完成！")
        else:
            self.set_error("snp_eff出错！", code="34505702")
            self.set_error("snp_eff出错！", code="34505705")
        os.system('rm {}'.format(file_path))

    def run(self):
        super(SnpeffTool, self).run()
        self.snp_eff()
        self.end()

    def make_new_snpEff_config(self):
        """
        转换snpeff中的文件路径
        data.dir =/mnt/lustre/users/sanger/app/database/dna_geneome/Oryza_sativa/NCBI/GCF_001433935.1/2015.10.10/
        ref.genome : ref
        :return:
        """
        self.logger.info("开始转换snpEff_config")
        path = ''
        with open(self.option("snpEff_config"), 'r') as r, open(self.work_dir + "/snpEff.config", 'w') as w1:
            for line in r:
                if re.match('data.*', line):
                    path = os.path.join(self.snpEff_config_base, line.strip().split('=')[1].split('app')[1].lstrip('/'))
                    break
                else:
                    pass
            self.logger.info("path:{}".format(path))
            w1.write("data.dir ={}\n".format(path))
            w1.write("ref.genome : ref\n")
        self.logger.info("转换snpEff_config完成")
