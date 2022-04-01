# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190402

import os
import json
from bson.objectid import ObjectId
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SsrCompareAgent(Agent):
    """
    ssr比较分析
    单条件批量比较参数：'{"allele_num": 2, "compare_type": "same", "sample_info": ["BDZ_VS_AA"],
      "region_select": {"chr1": "1-500000"}, "depth": "2-100", "gene_type": "same"}'
     多条件组合参数：'{"allele_num": 2, "sample_params": {"BDZ_VS_AA": {"AA": {"depth": "-100",
      "gene_type": "same"}, "compare_type": "diff", "BDZ": {"depth": "-100", "gene_type": "same"}}},
      "group_params": {"group1": {"sample": ["BDZ", "AA"], "depth": "-", "gene_fre": "-0.7", "miss_fre": "-"}},
      "region_select": "all"}'
    """
    def __init__(self, parent):
        super(SsrCompareAgent, self).__init__(parent)
        options = [
            {"name": "ssr_vcf", "type": "infile", "format": "wgs_v2.vcf", "required": True},
            {"name": "params_config", "type": "string"},
            {"name": "analysis_model", "type": "string", "required": True},  # 分析模式single or multiple
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ssr_vcf").is_set:
            raise OptionError("请设置SSR的vcf文件")
        if not self.option("params_config"):
            raise OptionError("请设置SSR比较的参数")
        if self.option("analysis_model") not in ["single", "multiple"]:
            raise OptionError("分析类型只能是single/multiple")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SsrCompareAgent, self).end()


class SsrCompareTool(Tool):
    def __init__(self, config):
        super(SsrCompareTool, self).__init__(config)
        self.python_path = "program/Python/bin/python"
        self.ssr_compare_path = self.config.PACKAGE_DIR + "/wgs_v2/ssr_compare.py"

    def get_params_config(self):
        """
        根据参数组建参数的config文件
        """
        self.params_config = os.path.join(self.work_dir, "params.json")
        with open(self.params_config, "w") as w:
            w.write(self.option("params_config") + "\n")

    def run_ssr_compare(self):
        """
        根据参数进行SSR比较分析
        """
        cmd = "{} {} -ssr {}".format(self.python_path, self.ssr_compare_path, self.option("ssr_vcf").prop["path"])
        cmd += " -c {} -m {} -o {}".format(self.params_config, self.option("analysis_model"), self.output_dir)
        command = self.add_command("ssr_compare", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("SSR比较分析运行完成")
        else:
            self.set_error("SSR比较分析运行失败")

    def set_db(self):
        self.logger.info("开始进行SSR比较分析结果导表")
        ssr_api = self.api.api("wgs_v2.ssr_specimen")
        compare_id = self.option("main_id")
        data_json = os.path.dirname(self.work_dir) + "/data.json"
        s3_output_dir = json.loads(open(data_json).read())["output"]
        subname, vcf_path = [], []
        for f in os.listdir(self.output_dir):
            if f.endswith("stat.xls"):
                name = f.split(".stat.xls")[0] if self.option("analysis_model") == "single" else None
                ssr_api.add_sg_ssr_compare_stat(compare_id=compare_id, ssr_stat_path=os.path.join(self.output_dir, f), name=name)
            if f.endswith("detail.xls"):
                name = f.split(".detail.xls")[0] if self.option("analysis_model") == "single" else None
                subname.append(name)
                vcf_path.append(os.path.join(s3_output_dir, f))
                ssr_api.add_sg_ssr_compare_detail(compare_id=compare_id, ssr_detail_path=os.path.join(self.output_dir, f), name=name)
        if self.option("analysis_model") == "single":
            update_dict = {"subname": subname, "vcf_path": vcf_path}
            ssr_api.update_info(coll="sg_ssr_compare", query_dict={"main_id": ObjectId(compare_id)}, update_dict=update_dict)
        else:
            ssr_api.update_info(coll="sg_ssr_compare", query_dict={"main_id": ObjectId(compare_id)}, update_dict={"vcf_path": vcf_path[0]})


    def run(self):
        super(SsrCompareTool, self).run()
        self.get_params_config()
        self.run_ssr_compare()
        if self.option("main_id"):
            self.set_db()
        self.end()
