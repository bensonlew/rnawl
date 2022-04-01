# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0517


from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import subprocess   # 79行！


class SsrRefPrimerDesignAgent(Agent):
    """
    软件：primer3_core，脚本：sample.ssr.p3in.pl、sample.ssr.p3out.pl
    对样本参考基因组的SSR进行引物设计
    """
    def __init__(self, parent):
        super(SsrRefPrimerDesignAgent, self).__init__(parent)
        options = [
            {"name": "misa_file", "type": "string"},  # scafSeq.misa文件
            {"name": "scafseq_file", "type": "string"},  # scafSeq文件  # ref.fa.misa 
            {"name": "tm1", "type": "float", "default": 57.0},  # float, Tm1 (℃)
            {"name": "tm2", "type": "float", "default": 63.0},  # float, Tm2 (℃),要大于tm1
            {"name": "product_size", "type": "string", "default": "300-500"},  # Product Size(bp),数值要为int,范围间用-分隔，多个用,分隔,如：100-300
            {"name": "primer_num", "type": "int", "default": 3},  # Max Pairs Primer Number,范围:[1,5]
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("misa_file"):
            raise OptionError("请设置misa_file文件", code="34506501")
        if not self.option("scafseq_file"): # &&
            raise OptionError("scafseq_file待设计引物fa文件不存在", code="34506502") # &&

    def set_resource(self):
        self._cpu = 5
        self._memory = "30G"

    def end(self):
        super(SsrRefPrimerDesignAgent, self).end()


class SsrRefPrimerDesignTool(Tool):
    def __init__(self, config):
        super(SsrRefPrimerDesignTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64")
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.ssr_p3in_path = self.config.PACKAGE_DIR + "/wgs/ref.ssr.p3in.pl" # 样本基因组组装用sample.ssr.p3in.pl
        self.ssr_p3out_path = self.config.PACKAGE_DIR + "/wgs/ref.ssr.p3out.pl"
        self.ssr_primersort_path = self.config.PACKAGE_DIR + "/wgs/sort.ssr.pl"        
        self.primer3_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/primer3/src/"

    def run_ssr_primer_in(self):
        """
        sample.ssr.p3in.pl
        """
        cmd = "{} {} -d {} ".format(self.perl_path, self.ssr_p3in_path, self.option("misa_file"))
        cmd += "-r {} -T1 {} -T2 {}".format(self.option("product_size"), self.option("tm1"), self.option("tm2"))
        cmd += " -p {} -ref {} ".format(self.option("primer_num"), self.option("scafseq_file"))
        cmd += " -o {}".format(self.work_dir)  # self.output_dir 缓存文件输出地方，output输出上层！！
        command = self.add_command("ssr_primer_in", cmd).run()
        self.wait() # 
        if command.return_code == 0:
            self.logger.info("ssr_primer_in完成，.p3in生成成功！")
        else:
            self.set_error("ssr_primer_in失败，.p3in生成失败！", code="34506501")

    def run_primer3_core(self):
        """
        ./primer3_core
        """
        # self.sample_name = os.path.basename(self.option("misa_file")).split(".")[0]
        self.sample_name = "ssr"
        p3in_file = os.path.join(self.work_dir, self.sample_name + ".p3in")
        cmd = "cd {} && ./primer3_core --output {} {}".format(self.primer3_path, self.work_dir + "/variation.p3out", 
                                                              p3in_file) # 生成SSR的p3out文件。
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("运行primer3_core完成")
        except subprocess.CalledProcessError:
            self.logger.info("运行primer3_core出错")
            self.set_error("运行primer3_core出错", code="34506502")
        os.system("cd {}".format(self.work_dir))   # 切换回路径

    def run_ssr_primer_out(self):
        """
        sample.ssr.p3out.pl
        """
        cmd = "{} {} {}".format(self.perl_path, self.ssr_p3out_path, os.path.join(self.work_dir, "variation.p3out"))
        command = self.add_command("ssr_primer_out", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ssr_primer_out完成")
        else:
            self.set_error("ssr_primer_out失败", code="34506503")

    def run_ssr_primer_sort(self):
        """
        sort.ssr.pl
        """
        cmd = "{} {} -r {} -s {}".format(self.perl_path, self.ssr_primersort_path, os.path.join(self.work_dir, "variation.result"), 
                                        os.path.join(self.output_dir, "ssr.ref.result.xls"))  # 补上.ssr
        command = self.add_command("ssr_primer_sort", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ssr_primer_sort完成")
        else:
            self.set_error("ssr_primer_sort失败", code="34506504")

    # def set_output(self):
    #     os.link(os.path.join(self.work_dir, "variation.result"), os.path.join(self.output_dir,  self.sample_name + ".result"))

    def run(self):
        super(SsrRefPrimerDesignTool, self).run()
        self.run_ssr_primer_in()
        self.run_primer3_core()
        self.run_ssr_primer_out()
        self.run_ssr_primer_sort()
        # self.set_output()
        self.end()
