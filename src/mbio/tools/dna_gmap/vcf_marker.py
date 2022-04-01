# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180605

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class VcfMarkerAgent(Agent):
    """
    vcf2marker.pl 将vcf(snp.indel)文件的标记转换成type:
    type:aaxaa aaxbb abxcc abxcd ccxab efxeg hkxhk lmxll nxxnp pmdep
    /mnt/ilustre/users/qingmei.cui/newmdt/sanger/6.gmap/new_pip_data_20180707/vcf2marker_data/vcf2marker.pl
    cd /mnt/ilustre/users/qingmei.cui/newmdt/sanger/6.gmap/new_pip_data_20180707/vcf2marker_data;
    perl /mnt/ilustre/users/long.huang/vcf2marker/vcf2marker.pl -vcf 10000_test_vcf -out pop -PID 14w_27 -MID 16S15
    """
    def __init__(self, parent):
        super(VcfMarkerAgent, self).__init__(parent)
        options = [
            {"name": "vcf_path", "type": "infile", "format": "dna_gmap.vcf"},  # 可以传入.vcf或.vcf.gz
            {"name": "pid", "type": "string"},
            {"name": "mid", "type": "string"},
            {"name": "primary_marker", "type": "outfile", "format": "dna_gmap.marker"}
        ]
        self.add_option(options)
        self.step.add_steps('VcfMarker')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.VcfMarker.start()
        self.step.update()

    def step_end(self):
        self.step.VcfMarker.finish()
        self.step.update()

    def check_options(self):
        if not self.option("vcf_path"):
            raise OptionError("请设置vcf文件", code="34801401")  # 必须有
        if not self.option("pid"):
            raise OptionError("请设置PID参数：父本名称", code="34801402")
        if not self.option("mid"):
            raise OptionError("请设置MID参数：母本名称", code="34801403")

    def set_resource(self):
        """
        运行所需资源
        vcf在边存边释放
        """
        self._cpu = 2
        self._memory = '30G'

    def end(self):
        super(VcfMarkerAgent, self).end()


class VcfMarkerTool(Tool):
    def __init__(self, config):
        super(VcfMarkerTool, self).__init__(config)
        # self.vcf2marker_path = self.config.PACKAGE_DIR + "/dna_gmap/vcf2markercui.pl"
        self.vcf2marker_path = self.config.PACKAGE_DIR + "/dna_gmap/vcf2marker.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl'

    def vcfmarker(self):
        """
        :return:
        """
        cmd = "{} {} -vcf {} -PID {} -MID {} -out {}".format(
            self.perl_path, self.vcf2marker_path, self.option("vcf_path").prop['path'],
            self.option("pid"), self.option("mid"),
            self.output_dir + "/pop.primary.marker")
        self.logger.info("^^^^^cmd:{}".format(cmd))
        self.logger.info("^^^^^开始进行VcfMarker")
        command = self.add_command("vcfmarker", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("VcfMarker完成！")
        else:
            self.set_error("VcfMarker出错！", code="34801401")
            raise Exception("VcfMarker出错！", code="34801402")
        self.option("primary_marker").set_path(self.output_dir + "/pop.primary.marker")

    def run(self):
        super(VcfMarkerTool, self).run()
        self.vcfmarker()
        self.end()
