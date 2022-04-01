## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdongxuan"

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re


class OtuAssociationAgent(Agent):
    """
    调用otu2shared.pl，将otu表转化为shared.txt，然后使用mothur计算物种之间的相似性，构建出物种相似性网络
    version v1.0
    author: hongdongxuan
    last_modify: 2016.12.07
    """
    def __init__(self, parent):
        super(OtuAssociationAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "lable", "type": "float", "default": 0.03}, # 设定lable（距离水平）的值
            {"name": "method", "type": "string", "default": "pearson"}  # 设定计算相似性的算法
        ]
        self.add_option(options)
        self.step.add_steps("OtuAssociation")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.OtuAssociation.start()
        self.step.update()

    def stepfinish(self):
        self.step.OtuAssociation.finish()
        self.step.update()


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("otutable"):
            raise OptionError("必须输入otutable文件", code="32704901")
        if self.option('lable') < 0 or self.option('lable') >= 1:
            raise OptionError('lable值必须是在[0,1)范围之内的值', code="32704902")
        if self.option('method') not in ['pearson', 'spearman', 'kendall']:
            raise OptionError('错误的物种相似性计算方式：%s', variables=(self.option('method')), code="32704903")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["shared.txt", "txt", "shared文件"]
                    ])
        result_dir.add_regexp_rules([
            [r"*\.otu\.corr", "corr", "物种相似性结果文件"],

        ])
        super(OtuAssociationAgent, self).end()


class OtuAssociationTool(Tool):
    """
    将otu转为shared.txt，使用mothur就算物种间相似性
    perl otu2shared_xhd.pl -i otu.new.xls -l 0.03 -o shared.txt
    mothur "#otu.association(shared=shared.txt, method=pearson)"
    """
    def __init__(self, config):
        super(OtuAssociationTool, self).__init__(config)
        self._version = '1.0.1'
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '
        self.mothur_path = 'bioinfo/meta/mothur-1.30/mothur.1.30'
        self.script_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/meta/scripts/otu2shared_xhd.pl")

    def run_otu2shared(self):
        one_cmd = self.perl_path + self.script_path + " -i %s -l %s -o %s" % (self.option('otutable').prop["path"], self.option('lable'), "shared.txt")
        self.logger.info(one_cmd)
        self.logger.info("开始运行one_cmd")
        cmd = self.add_command("one_cmd", one_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行one_cmd成功")
        else:
            self.logger.info("运行one_cmd出错")

    def run_OtuAssociation(self):
        two_cmd = self.mothur_path + " \"#otu.association(shared=shared.txt, method=%s)\"" % (self.option('method'))
        self.logger.info(two_cmd)
        self.logger.info("开始运行two_cmd")
        cmd = self.add_command("two_cmd", two_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行two_cmd成功")
        else:
            self.logger.info("运行two_cmd失败")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir + '/')
        for f in results:
            if re.search(r'.*\.corr$', f):  # modified by hongdongxuan 20170324
                # os.link(self.work_dir + '/' + f, self.output_dir + "/" + f)
                #os.link(self.work_dir + '/' + f, self.output_dir + "/shared.0.03." + self.option('method') + ".corr")
                os.link(self.work_dir + '/' + f, self.output_dir + "/shared." + str(self.option('lable')) + "." + self.option('method') + ".corr")  #modified by zhengyuan 20171121
            elif str(f) == "shared.txt":
                os.link(self.work_dir + '/' + f, self.output_dir + "/" + f)
            else:
                pass
        self.logger.info('设置文件夹路径成功')


    def run(self):
        super(OtuAssociationTool, self).run()
        self.run_otu2shared()
        self.run_OtuAssociation()
        self.set_output()
        self.end()
