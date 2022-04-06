# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180411

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SampleCompareAgent(Agent):
    """
    snp/indel比较分析接口中用于样本比较分析, 所有参数多选的时候，都以英文逗号分隔
    """
    def __init__(self, parent):
        super(SampleCompareAgent, self).__init__(parent)
        options = [
            {"name": "sample1", "type": "string"},
            {"name": "sample2", "type": "string"},
            {"name": "is_same", "type": "string", "default": "true"},  # 页面样本间拷贝数变异是否相同
            {"name": "vcf_file", "type": "string"},
            {"name": "funtype", "type": "string"},   # HIGH,LOW,MODIFIER,MODERATE，没有值的时候前端传空
            {"name": "efftype", "type": "string"},   # efftype 没有值的时候前端传空
            {"name": "len1", "type": "string", "default": "1"},  # snp的时候默认len1与len2是1, 页面上变异长度开始与结束一致
            {"name": "len2", "type": "string", "default": "1"},
            {"name": "dep1", "type": "string"},  # 2,5
            {"name": "dep2", "type": "string"},  # 2,6
            {"name": "location", "type": "string"}  # chr1,1,200
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
            raise OptionError("缺少样本1参数", code="34504701")
        if not self.option("sample2"):
            raise OptionError("缺少样本2参数", code="34504702")
        if not self.option("is_same"):
            raise OptionError("缺少is_same参数", code="34504703")
        if not self.option("vcf_file"):
            raise OptionError("缺少vcf_file参数", code="34504704")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '10G'

    def end(self):
        super(SampleCompareAgent, self).end()


class SampleCompareTool(Tool):
    def __init__(self, config):
        super(SampleCompareTool, self).__init__(config)
        self.script_path = self.config.PACKAGE_DIR + "/wgs/sample_vcf.pl"
        self.script_path1 = self.config.PACKAGE_DIR + "/wgs/single_vcf.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '
        self.python = "miniconda2/bin/python"
        self.distribution_graph_path = self.config.PACKAGE_DIR + "/wgs/distribution_graph.py"

    def sampletosample(self):
        """
        脚本中的len1指的是样本1的突变位点长度格式1,1； len2是指的是样本2的突变长度格式1,1，该接口中样本1与样
        本2的长度一样
        :return:
        """
        cmd = "{}{} -sample1 {} -sample2 {} -vcf {} -out {}"\
            .format(self.perl_path, self.script_path, self.option("sample1"), self.option("sample2"),
                    self.option("vcf_file"), self.output_dir)
        if self.option("is_same") == "true":
            cmd += " -equal {}".format(1)
        else:
            cmd += " -equal {}".format(0)
        if self.option("len1") == '' and self.option("len2") != '':
            cmd += " -len1 {} -len2 {}".format(",".join(["-1", self.option("len2")]),
                                               ",".join(["-1", self.option("len2")]))
        elif self.option("len1") != '' and self.option("len2") == '':
            cmd += " -len1 {} -len2 {}".format(",".join([self.option("len1"), "1000000000000000000000000000000"]),
                                               ",".join([self.option("len1"), "1000000000000000000000000000000"]))
        elif self.option("len1") != '' and self.option("len2") != '':
            cmd += " -len1 {} -len2 {}".format(",".join([self.option("len1"), self.option("len2")]),
                                               ",".join([self.option("len1"), self.option("len2")]))
        if self.option("funtype") == "":
            cmd += " -funtype {}".format("all")
        else:
            cmd += " -funtype {}".format(self.option("funtype"))
        if self.option("efftype") == "":
            cmd += " -efftype {}".format("all")
        else:
            cmd += " -efftype {}".format(self.option("efftype"))
        if self.option("dep1") != ",":   # 当为，的时候前端没有传值
            temp = self.option("dep1").split(',')
            if temp[0] == "" and temp[1] != "":
                cmd += " -dep1 {}".format(",".join(["-1", temp[1]]))
            elif temp[0] != "" and temp[1] == "":
                cmd += " -dep1 {}".format(",".join([temp[0], "100000000"]))
        if self.option("dep2") != ",":
            temp = self.option("dep2").split(',')
            if temp[0] == "" and temp[1] != "":
                cmd += " -dep2 {}".format(",".join(["-1", temp[1]]))
            elif temp[0] != "" and temp[1] == "":
                cmd += " -dep2 {}".format(",".join([temp[0], "100000000"]))
        if self.option("location") != ",,":
            temp = self.option("location").split(',')
            if temp[0] == "" or temp[1] == "" or temp[2] == "":
                pass
            else:
                cmd += " -pos {}".format(self.option("location"))
        self.logger.info(cmd)
        self.logger.info("开始进行sample_compare")
        command = self.add_command("sample_compare", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("sample_compare完成！")
        else:
            self.set_error("sample_compare出错！", code="34504701")
            self.set_error("sample_compare出错！", code="34504710")

    def sample2ref(self):
        """
        针对样本与ref的比较，原本是想直接在sample2sample中进行判断的，但是上述函数已经写的很乱了，还是算了，重新定义一个函数
        :return:
        """
        if self.option("sample1") == "reference":
            sample = self.option("sample2")
        else:
            sample = self.option("sample1")
        cmd = "{}{} -sample {} -vcf {} -out {}" \
            .format(self.perl_path, self.script_path1, sample, self.option("vcf_file"), self.output_dir)
        if self.option("is_same") == "true":
            cmd += " -equal {}".format(1)
        else:
            cmd += " -equal {}".format(0)
        if self.option("len1") == '' and self.option("len2") != '':
            cmd += " -len1 {}".format(",".join(["-1", self.option("len2")]))
        elif self.option("len1") != '' and self.option("len2") == '':
            cmd += " -len1 {}".format(",".join([self.option("len1"), "1000000000000000000000000000000"]))
        elif self.option("len1") != '' and self.option("len2") != '':
            cmd += " -len1 {}".format(",".join([self.option("len1"), self.option("len2")]))
        if self.option("funtype") == "":
            cmd += " -funtype {}".format("all")
        else:
            cmd += " -funtype {}".format(self.option("funtype"))
        if self.option("efftype") == "":
            cmd += " -efftype {}".format("all")
        else:
            cmd += " -efftype {}".format(self.option("efftype"))
        if self.option("sample1") == "reference":
            if self.option("dep1") != ",":  # 当为，的时候前端没有传值
                temp = self.option("dep1").split(',')
                if temp[0] == "" and temp[1] != "":
                    cmd += " -dep1 {}".format(",".join(["-1", temp[1]]))
                elif temp[0] != "" and temp[1] == "":
                    cmd += " -dep1 {}".format(",".join([temp[0], "100000000"]))
        elif self.option("sample2") == "reference":
            if self.option("dep2") != ",":  # 当为，的时候前端没有传值
                temp = self.option("dep2").split(',')
                if temp[0] == "" and temp[1] != "":
                    cmd += " -dep1 {}".format(",".join(["-1", temp[1]]))
                elif temp[0] != "" and temp[1] == "":
                    cmd += " -dep1 {}".format(",".join([temp[0], "100000000"]))
        if self.option("location") != ",,":
            temp = self.option("location").split(',')
            if temp[0] == "" or temp[1] == "" or temp[2] == "":
                pass
            else:
                cmd += " -pos {}".format(self.option("location"))
        self.logger.info(cmd)
        self.logger.info("开始进行sample_compare")
        command = self.add_command("sample_compare", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("sample_compare完成！")
        else:
            self.set_error("sample_compare出错！", code="34504702")
            self.set_error("sample_compare出错！", code="34504711")

    def distribution_graph(self):
        """
        将win.stat进行滑窗，用于画染色体分布图
        """
        win_stat = self.output_dir + "/win.stat"
        step_win_stat = self.work_dir + "/win.stat.xls"
        cmd = "{} {} -step {} -i {} -o {}".format(self.python, self.distribution_graph_path, 100000,
                                                  win_stat, step_win_stat)
        self.logger.info("开始进行染色体分布图滑窗")
        command = self.add_command("distribution_graph", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("distribution_graph完成！")
        else:
            self.set_error("distribution_graph出错！", code="34504703")
            self.set_error("distribution_graph出错！", code="34504712")

    def run(self):
        super(SampleCompareTool, self).run()
        if self.option("sample1") == "reference" or self.option("sample2") == "reference":
            self.sample2ref()
        else:
            self.sampletosample()
        self.distribution_graph()
        self.end()
