# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190222

import os
import json
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SvCompareAgent(Agent):
    """
    sv比较分析
    """
    def __init__(self, parent):
        self._sheet = parent
        super(SvCompareAgent, self).__init__(parent)
        options = [
            {"name": "sv_vcf", "type": "infile", "format": "wgs_v2.vcf", "required": True},  # sv的vcf文件
            {"name": "params_config", "type": "string", "required": True},  # sv比较分析的参数
            {"name": "analysis_model", "type": "string", "required": True},  # 分析模式single or multiple
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sv_vcf").is_set:
            raise OptionError("请设置sv的vcf文件")
        if not self.option("params_config"):
            raise OptionError("请设置参数params_config")
        if self.option("analysis_model") not in ["single", "multiple"]:
            raise OptionError("分析模式只能是single or multiple,请检查")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SvCompareAgent, self).end()


class SvCompareTool(Tool):
    def __init__(self, config):
        super(SvCompareTool, self).__init__(config)
        self.python_path = "program/Python/bin/python"
        self.sv_compare_path = self.config.PACKAGE_DIR + "/wgs_v2/sv_compare.py"

    def get_sv_config(self):
        """
        根据参数组建参数的config文件
        """
        self.sv_comfig = os.path.join(self.work_dir, "params.json")
        with open(self.sv_comfig, "w") as w:
            w.write(self.option("params_config") + "\n")

    def run_sv_compare(self):
        cmd = "{} {} -sv {} ".format(self.python_path, self.sv_compare_path, self.option("sv_vcf").prop["path"])
        cmd += "-c {} -o {} -m {}".format(self.sv_comfig, self.output_dir, self.option("analysis_model"))
        command = self.add_command("sv_compare", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("SV比较分析运行完成")
        else:
            self.set_error("SV比较分析运行失败，请检查")

    def set_db(self):
        self.logger.info("开始进行SV比较分析导表")
        sv_api = self.api.api("wgs_v2.sv_compare")
        for f in os.listdir(self.output_dir):
            if f.endswith("summary.xls"):
                sv_api.add_sg_sv_compare_stat(compare_id=self.option("main_id"), summary_path=os.path.join(self.output_dir, f))
            elif f.endswith("detail.xls"):
                sv_api.add_sg_sv_compare_detail(compare_id=self.option("main_id"), detail_path=os.path.join(self.output_dir, f), type=self.option("analysis_model"))
        data_json = os.path.dirname(self.work_dir) + "/data.json"
        s3_output_dir = json.loads(open(data_json).read())["output"]
        sv_api.update_vcf_path(compare_id=self.option("main_id"), output_dir=self.output_dir, s3_output_dir=s3_output_dir)

    def run(self):
        super(SvCompareTool, self).run()
        self.get_sv_config()
        self.run_sv_compare()
        if self.option("main_id"):
            self.set_db()
        self.end()
