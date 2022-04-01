# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.25

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SsrMisaAgent(Agent):
    """
    软件：misa.pl
    """
    def __init__(self, parent):
        super(SsrMisaAgent, self).__init__(parent)
        options = [
            {"name": "scafseq_file", "type": "string"},  # scafSeq文件
            {"name": "needini", "type": "bool", "default": False},  # 做misa的时候是否要配置ini文件，默认不配
            {"name": "rept_1", "type": "int", "default": 10},
            {"name": "rept_2", "type": "int", "default": 6},
            {"name": "rept_3", "type": "int", "default": 5},
            {"name": "rept_4", "type": "int", "default": 5},
            {"name": "rept_5", "type": "int", "default": 5},
            {"name": "rept_6", "type": "int", "default": 5},
            {"name": "ssr_distance", "type": "int", "default": 100},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("scafseq_file"):
            raise OptionError("请设置scafseq_file文件", code="34506301")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(SsrMisaAgent, self).end()


class SsrMisaTool(Tool):
    def __init__(self, config):
        super(SsrMisaTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.misa_path = self.config.PACKAGE_DIR + "/wgs/misa.pl"
        self.misa_path2 = self.config.PACKAGE_DIR + "/wgs_v2/misa.pl"
        self.ssr_type_stat = self.config.PACKAGE_DIR + "/wgs/SSR.type.stat.pl"
        self.misa_ini = self.config.SOFTWARE_DIR + "/bioinfo/WGS/misa.ini"

    def run_misa(self):
        """
        misa.pl
        """
        if self.option("needini"):
            self.misa_ini = self.work_dir + '/misa.ini'
            misa_path = self.work_dir + "/misa.pl"
            os.link(self.misa_path2, misa_path)
        else:
            misa_path = self.misa_path
        scafseq_path = self.work_dir + "/" + os.path.basename(self.option("scafseq_file"))
        self.checklink(scafseq_path)
        os.link(self.option("scafseq_file"), scafseq_path)
        cmd = "{} {} {} {}".format(self.perl_path, misa_path, scafseq_path, self.misa_ini)
        command = self.add_command("ssr_misa", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ssr_misa完成")
        else:
            self.set_error("ssr_misa失败", code="34506301")
        misa_file = os.path.basename(scafseq_path) + ".misa"
        self.checklink(os.path.join(self.output_dir, misa_file))
        os.link(os.path.join(self.work_dir, misa_file), os.path.join(self.output_dir, misa_file))
        self.checklink(os.path.join(self.output_dir, os.path.basename(self.option("scafseq_file"))))
        os.link(self.option("scafseq_file"), os.path.join(self.output_dir, os.path.basename(self.option("scafseq_file"))))

    def make_misa_ini(self):
        """
        解析页面传过来的参数，然后写成标准格式的misa.ini文件,并且要确保它和misa.pl在一个文件夹里
        字典解析
        :param kwargs:
        :return:
        """
        with open(self.work_dir + '/misa.ini','w') as ini:
            param1 = '1' + '-' + str(self.option("rept_1"))
            param2 = '2' + '-' + str(self.option("rept_2"))
            param3 = '3' + '-' + str(self.option("rept_3"))
            param4 = '4' + '-' + str(self.option("rept_4"))
            param5 = '5' + '-' + str(self.option("rept_5"))
            param6 = '6' + '-' + str(self.option("rept_6"))
            param7 = str(self.option("ssr_distance"))
            line1_ideography = 'definition(unit_size,min_repeats):'
            line2_ideography = 'interruptions(max_difference_between_2_SSRs):'
            line1 = [line1_ideography, param1, param2, param3, param4, param5, param6]
            line2 = [line2_ideography, param7]
            ini.write(" ".join(tuple(line1)) + '\n')
            ini.write(" ".join(tuple(line2)) + '\n')

    def checklink(self,filepath):
        if os.path.exists(filepath):
            os.remove(filepath)
        return True

    def run(self):
        super(SsrMisaTool, self).run()
        if self.option("needini"):
            self.make_misa_ini()
        self.run_misa()
        self.end()

# # -*- coding: utf-8 -*-
# # __author__ = 'zengjing' # qingmei
# # modified 2018.04.25

# from biocluster.agent import Agent
# from biocluster.tool import Tool
# from biocluster.core.exceptions import OptionError
# import os


# class SsrMisaAgent(Agent):
#     """
#     软件：misa.pl
#     """
#     def __init__(self, parent):
#         super(SsrMisaAgent, self).__init__(parent)
#         options = [
#             {"name": "scafseq_file", "type": "string"},  # scafSeq文件
#         ]
#         self.add_option(options)

#     def check_options(self):
#         if not self.option("scafseq_file"):
#             raise OptionError("请设置scafseq_file文件")

#     def set_resource(self):
#         self._cpu = 2
#         self._memory = "5G"

#     def end(self):
#         super(SsrMisaAgent, self).end()


# class SsrMisaTool(Tool):
#     def __init__(self, config):
#         super(SsrMisaTool, self).__init__(config)
#         self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
#         self.misa_path = self.config.PACKAGE_DIR + "/wgs/misa.pl"
#         self.ssr_type_stat = self.config.PACKAGE_DIR + "/wgs/SSR.type.stat.pl"
#         self.misa_ini = self.config.SOFTWARE_DIR + "/bioinfo/WGS/misa.ini"

#     def run_misa(self):
#         """
#         misa.pl
#         """
#         scafseq_path = self.work_dir + "/" + os.path.basename(self.option("scafseq_file"))
#         if os.path.exists(scafseq_path):
#             os.remove(scafseq_path)
#         os.link(self.option("scafseq_file"), scafseq_path)  # 连接输入文件到work_dir下

#         cmd = "{} {} {} {}".format(self.perl_path, self.misa_path, scafseq_path, self.misa_ini)     # 执行程序cmd
#         command = self.add_command("ssr_misa", cmd).run()   # 跑程序
#         self.wait() # 等待程序运行结束
#         if command.return_code == 0:
#             self.logger.info("ssr_misa完成")
#         else:
#             self.set_error("ssr_misa失败") # 提示信息
#     misa_file = os.path.basename(self.option("scafseq_file")) + ".misa"
#     j = os.path.join(self.output_dir, misa_file)
#     if os.path.exists(j):
#             os.remove(j)
#         os.link(os.path.join(self.work_dir, misa_file), j)
#     i = os.path.join(self.output_dir, os.path.basename(self.option("scafseq_file")))
#     if os.path.exists(i):
#             os.remove(i)
#         os.link(self.option("scafseq_file"), i)

#     def run(self):
#         super(SsrMisaTool, self).run()
#         self.run_misa()
#         self.end()
