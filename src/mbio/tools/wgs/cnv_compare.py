# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180408

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class CnvCompareAgent(Agent):
    """
    CnvDiff用于交互分析中cnv变异检测--cnv比较分析
    """
    def __init__(self, parent):
        super(CnvCompareAgent, self).__init__(parent)
        options = [
            {"name": "sample1", "type": "string"},
            {"name": "sample2", "type": "string"},
            {"name": "is_same", "type": "string", "default": "true"},  # 页面样本间拷贝数变异是否相同
            {"name": "variation_type", "type": "string", "default": "deletion,duplication"},  # 变异位点类型
            {"name": "variation_len", "type": "string"}  # 变异区域长度,页面不选传到这里就是0-0，但是传给脚本的时候该参数不传
        ]
        self.add_option(options)
        self.step.add_steps('cnvdiff')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cnvdiff.start()
        self.step.update()

    def step_end(self):
        self.step.cnvdiff.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("sample1"):
            raise OptionError("缺少样本1参数", code="34501501")
        if not self.option("sample2"):
            raise OptionError("缺少样本2参数", code="34501502")
        if not self.option("is_same"):
            raise OptionError("缺少is_same参数", code="34501503")
        if not self.option("variation_type"):
            raise OptionError("缺少variation_type参数", code="34501504")
        else:
            if self.option("variation_type") not in ["deletion,duplication", "deletion", "duplication"]:
                raise OptionError("变异位点类型%s不合法！",variables=(self.option("variation_type")), code="34501505")
        if not self.option("variation_len"):
            raise OptionError("缺少variation_len参数", code="34501506")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '10G'
        
    def end(self):
        super(CnvCompareAgent, self).end()


class CnvCompareTool(Tool):
    def __init__(self, config):
        super(CnvCompareTool, self).__init__(config)
        self.cnv_diff_path = self.config.PACKAGE_DIR + "/wgs/cnv.diff.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '
        
    def cnv_diff(self):
        """
        要重新写下
        :return:
        """
        cmd = "{}{} -cnv1 {} -cnv2 {} -b {} -o {}"\
            .format(self.perl_path, self.cnv_diff_path, self.option("sample1"), self.option("sample2"),
                    self.option("is_same"), self.output_dir)
        if self.option("variation_len") != "-":
            cmd += " -l {}".format(self.option("variation_len"))
        if self.option("variation_type") == "":
            cmd += " -t {}".format("deletion,duplication")
        else:
            cmd += " -t {}".format(self.option("variation_type"))
        self.logger.info(cmd)
        self.logger.info("开始进行cnv_diff")
        command = self.add_command("cnv_diff", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cnv_diff完成！")
        else:
            self.set_error("cnv_diff出错！", code="34501501")
            self.set_error("cnv_diff出错！", code="34501504")

    def run(self):
        super(CnvCompareTool, self).run()
        self.cnv_diff()
        self.end()
