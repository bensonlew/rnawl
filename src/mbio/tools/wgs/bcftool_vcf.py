# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180408

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class BcftoolVcfAgent(Agent):
    """
    author = HONGDONG
    version = 1.0
    last modify = 20180408
    该tool用于对vcf进行处理，首先将分割后的vcf文件使用concat进行合并，然后annotate，最后在进行sort排序生成最后需要的vcf文件
    """
    def __init__(self, parent):
        super(BcftoolVcfAgent, self).__init__(parent)
        options = [
            {"name": "vcf_list", "type": "string"}  # 所有分割后的vcf文件列表
        ]
        self.add_option(options)
        self.step.add_steps('vcf_make')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.vcf_make.start()
        self.step.update()

    def step_end(self):
        self.step.vcf_make.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("vcf_list"):
            raise OptionError("缺少vcf_list参数", code="34500501")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '200G'
        
    def end(self):
        super(BcftoolVcfAgent, self).end()


class BcftoolVcfTool(Tool):
    def __init__(self, config):
        super(BcftoolVcfTool, self).__init__(config)
        self._version = '1.0.1'
        self.bcftools = "bioinfo/seq/bcftools-1.7/"
        
    def bcftools_concat(self):
        """
        进行vcf文件合并
        bcftools concat  1.noid.vcf  2.noid.vcf  3.noid.vcf  4.noid.vcf 5.noid.vcf 6.noid.vcf 7.noid.vcf -o pop.noid.vcf
            -O v
        """
        vcf_list = self.get_vcf_list()
        cmd1 = "{}bcftools concat {}-o {}/pop.noid.vcf -O v"\
            .format(self.bcftools, vcf_list, self.work_dir)
        self.logger.info(cmd1)
        self.logger.info("开始进行bcftools_concat")
        command1 = self.add_command("cmd1", cmd1).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("bcftools_concat完成！")
        else:
            self.set_error("bcftools_concat出错！", code="34500501")
            self.set_error("bcftools_concat出错！", code="34500510")

    def bcftools_annotate(self):
        """
        进行vcf文件注释
        bcftools annotate --set-id +'%CHROM\_%POS' pop.noid.vcf -o pop.nosort.vcf
        :return:
        """
        cmd2 = "{}bcftools annotate --set-id +'%CHROM\_%POS' {} -o {}/pop.nosort.vcf"\
            .format(self.bcftools, self.work_dir + '/pop.noid.vcf', self.work_dir)
        self.logger.info(cmd2)
        self.logger.info("开始进行bcftools_annotate")
        command2 = self.add_command("cmd2", cmd2).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("bcftools_annotate完成！")
        else:
            self.set_error("bcftools_annotate出错！", code="34500502")
            self.set_error("bcftools_annotate出错！", code="34500511")

    def bcftools_sort(self):
        """
        进行vcf文件排序
        bcftools sort -m 100G -T temp/ -o pop.variant.vcf pop.nosort.vcf
        :return:
        """
        cmd3 = "{}bcftools sort -m 100G -T {}/temp/ -o {}/pop.variant.vcf {}"\
            .format(self.bcftools, self.work_dir, self.output_dir, self.work_dir + "/pop.nosort.vcf")
        self.logger.info(cmd3)
        self.logger.info("开始进行bcftools_sort")
        command3 = self.add_command("cmd3", cmd3).run()
        self.wait(command3)
        if command3.return_code == 0:
            self.logger.info("bcftools_sort完成！")
        else:
            self.set_error("bcftools_sort出错！", code="34500503")
            self.set_error("bcftools_sort出错！", code="34500512")

    def get_vcf_list(self):
        """
        获取vcf列表参数文档
        :return:
        """
        vcf_list = ""
        with open(self.option("vcf_list"), 'r') as r:
            data = r.readlines()
            for line in data:
                line = line.strip()
                vcf_list += line
                vcf_list += " "
        return vcf_list

    def run(self):
        super(BcftoolVcfTool, self).run()
        self.bcftools_concat()
        self.bcftools_annotate()
        self.bcftools_sort()
        self.end()
