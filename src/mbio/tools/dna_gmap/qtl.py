# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.26

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class QtlAgent(Agent):
    """
    遗传图谱：QTL定位
    """
    def __init__(self, parent):
        super(QtlAgent, self).__init__(parent)
        options = [
            {"name": "trait_file", "type": "infile", "format": "dna_gmap.trait"},  # 性状文件
            {"name": "total_map", "type": "infile", "format": "dna_gmap.lg"},  # pop_type为F2时是total.csv，CP的时候是total.male.map/total.female.map/total.sexAver.map
            {"name": "total_phase", "type": "infile", "format": "dna_gmap.lg"},  # pop_type为CP时才有，total.male.phase /total.female.phase /total.sexAver.phase
            {"name": "pop_type", "type": "string", "default": "F2"},  # pop type, CP/F2
            {"name": "location_method", "type": "string", "default": "cim"},  # 定位方法,cim/scanone
            # {"name": "trait_type", "type": "bool", "default": False},  # 性状文件的格式, 非二进制/二进制
            {"name": "pm_num", "type": "int", "default": 1000},  # 置换检验次数，[0,+∞]
            {"name": "p_value", "type": "float"},  # 阈值，p value,默认值0.05
            {"name": "lod_value", "type": "int"},  # 阈值，LOD,默认值3
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("trait_file").is_set:
            raise OptionError("请设置性状文件trait_file", code="34801501")
        if not self.option("total_map").is_set:
            raise OptionError("请设置map的结果文件", code="34801502")
        if self.option("pop_type") not in ["F1", "F2", "DH", "BC", "Ril"]:
            raise OptionError("pop_type类型只能为F1/F2/DH/BC/Ril,而不能是%s", variables=(self.option("pop_type")), code = "34801503")
        if self.option("pop_type") == "F1":
            if not self.option("total_phase").is_set:
                raise OptionError("请设置map对应的phase", code="34801504")
        if self.option("location_method") not in ["cim", "scanone"]:
            raise OptionError("定位方法出错，只能为cim/scanone", code="34801505")
        if self.option("pm_num") < 0:
            raise OptionError("置换检验次数只能是大于等于0的正整数", code="34801506")
        if not self.option("p_value") and not self.option("lod_value"):
            raise OptionError("必须设置阈值P_value或LOD", code="34801507")
        if self.option("p_value"):
            if self.option("p_value") > 1 or self.option("p_value") < 0:
                raise OptionError("阈值p_value的范围是[0,1]", code="34801508")
        if self.option("lod_value"):
            if self.option("lod_value") < 0:
                raise OptionError("阈值lod_value的范围是[0,+∞]", code="34801509")
        if self.option("pop_type") == "F1" and self.option("location_method") == "cim":
            raise OptionError("群体类型为F1的时候定位方法不能是cim,请检查", code="34801510")

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(QtlAgent, self).end()


class QtlTool(Tool):
    def __init__(self, config):
        super(QtlTool, self).__init__(config)
        self.rscript = "program/R-3.3.3/bin/Rscript"
        self.qtl_nocp = self.config.PACKAGE_DIR + "/dna_gmap/qtl-NOCP.R"
        self.qtl_cp = self.config.PACKAGE_DIR + "/dna_gmap/qtl-CP.R"
        self.btl_nocp = self.config.PACKAGE_DIR + "/dna_gmap/btl-NOCP.R"
        self.btl_cp = self.config.PACKAGE_DIR + "/dna_gmap/btl-CP.R"
        if self.option("pop_type") == "DH" or self.option("pop_type") == "BC":
            self.pop_type = "bc"
        elif self.option("pop_type") == "F2":
            self.pop_type = "f2"
        elif self.option("pop_type") == "Ril":
            self.pop_type = "riself"
        else:
            self.pop_type = "F1"

    def run_qtl_cp(self):
        """
        pop_type是CP的情况
        """
        if self.option("trait_file").prop["trait_type"] == "btl":
            self.cp = self.btl_cp
        else:
            self.cp = self.qtl_cp
        cmd = "{} {} --map {}".format(self.rscript, self.cp, self.option("total_map").prop["path"])
        cmd += " --trt {} --loc {}".format(self.option("trait_file").prop["path"], self.option("total_phase").prop["path"])
        cmd += " --out {} --num {} --method {}".format(self.work_dir, self.option("pm_num"), self.option("location_method"))
        if self.option("p_value"):
            cmd += " --pvalue {}".format(self.option("p_value"))
        if self.option("lod_value"):
            cmd += " --lod {}".format(self.option("lod_value"))
        self.run_cmd(cmd, "qtl_cp")

    def run_qtl_nocp(self):
        """
        pop_type是F2的情况
        """
        if self.option("trait_file").prop["trait_type"] == "btl":
            self.nocp = self.btl_nocp
        else:
            self.nocp = self.qtl_nocp
        cmd = "{} {} --mark {}".format(self.rscript, self.nocp, self.option("total_map").prop["path"])
        cmd += " --trt {} --out {}".format(self.option("trait_file").prop["path"], self.work_dir)
        cmd += " --method {} --num {}".format(self.option("location_method"), self.option("pm_num"))
        cmd += " --pop {}".format(self.pop_type)
        # cmd += " --pop bcsft --bc 2 --f 2"
        if self.option("p_value"):
            cmd += " --pvalue {}".format(self.option("p_value"))
        if self.option("lod_value"):
            cmd += " --lod {}".format(self.option("lod_value"))
        self.run_cmd(cmd, "qtl_nocp")

    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("qtl定位失败", code="34801501")
            self.set_error("qtl运行失败", code="34801504")

    def set_output(self):
        no_qtl = True
        trait_name = os.path.basename(self.option("trait_file").prop["path"]).split(".txt")[0]
        for f in os.listdir(self.work_dir):
            old = os.path.join(self.work_dir, f)
            if f.endswith("qtl.csv"):
                new = os.path.join(self.output_dir, trait_name + ".qtl.csv")
                if os.path.exists(new):
                    os.remove(new)
                os.link(old, new)
                no_qtl = False
            if f.endswith("scan.csv"):
                new = os.path.join(self.output_dir, trait_name + ".scan.csv")
                if os.path.exists(new):
                    os.remove(new)
                os.link(old, new)
        if no_qtl:
            self.logger.info("qtl定位失败")

    def run(self):
        super(QtlTool, self).run()
        if self.pop_type == "F1":
            self.run_qtl_cp()
        else:
            self.run_qtl_nocp()
        self.set_output()
        self.end()
