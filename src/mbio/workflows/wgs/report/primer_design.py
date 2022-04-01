# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180608

import os
import re
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class PrimerDesignWorkflow(Workflow):
    """
    交互分析：样本基因组SSR分析及引物设计
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PrimerDesignWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "ref_fa", "type": "string"},
            {"name": "diff_variant", "type": "infile", 'format': 'bsa.vcf'},  # 比较分析的结果diff.variant
            {"name": "tm1", "type": "float", "default": 57.0},  # float, Tm1 (℃)
            {"name": "tm2", "type": "float", "default": 63.0},  # float, Tm2 (℃),要大于tm1
            {"name": "product_size", "type": "string"},  # Product Size(bp),数值要为int,范围间用-分隔，多个用,分隔,如：100-300
            {"name": "primer_num", "type": "int"},  # Max Pairs Primer Number,范围:[1,5]
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("请设置ref_fa文件", code="14500701")
        if not self.option("diff_variant"):
            raise OptionError("请设置比较分析的结果diff_variant文件", code="14500702")
        if self.option("tm1") >= self.option("tm2"):
            raise OptionError("tm1要小于tm2", code="14500703")
        if not self.option("product_size"):
            raise OptionError("请设置product_size", code="14500704")
        if not re.search(r"(\d+)-(\d+).*", self.option("product_size")):
            raise OptionError("%s product_size格式不正确,用-分隔", variables=(self.option("product_size")), code="14500705")
            # raise OptionError("{} product_size格式不正确,用-分隔".format(self.option("product_size")))
        if not self.option("primer_num"):
            raise OptionError("请设置primer_num", code="14500706")
        if self.option("primer_num") > 5 or self.option("primer_num") < 1:
            raise OptionError("primer_num:%s范围不为[1,5]", variables=(self.option("primer_num")), code="14500707")
            # raise OptionError("primer_num:{}范围不为[1,5]".format(self.option("primer_num")))

    def run_primer_design(self):
        options = {
            "ref_fa": self.option("ref_fa"),
            "diff_variant": self.option("diff_variant").prop['path'],
            "tm1": self.option("tm1"),
            "tm2": self.option("tm2"),
            "product_size": self.option("product_size"),
            "primer_num": self.option("primer_num")
        }
        self.primer_design = self.add_tool("wgs.primer_design")
        self.primer_design.set_options(options)
        self.primer_design.on("end", self.set_output)
        self.primer_design.run()

    def set_output(self):
        for f in os.listdir(self.primer_design.output_dir):
            f1 = os.path.join(self.primer_design.output_dir, f)
            f2 = os.path.join(self.output_dir, f)
            if os.path.exists(f2):
                os.remove(f2)
            os.link(f1, f2)
        self.set_db()
        self.end()

    def set_db(self):
        """
        将结果保存到mongo
        """
        self.logger.info("将结果保存到mongo")
        primer_id = self.option("main_id")
        primer_api = self.api.api("wgs.primer_design")
        if self.option("project_type"):
            primer_api._project_type = self.option("project_type")
        primer_result = os.path.join(self.work_dir, "variation.result")
        primer_result = os.path.join(self.primer_design.work_dir, "variation.result")
        download_file = self._sheet.output.rstrip('/') + "/variant_result.xls"
        self.logger.info(download_file)
        primer_api.add_sg_primer_detail(primer_id, primer_result, download_file)

    def run(self):
        self.run_primer_design()
        super(PrimerDesignWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(PrimerDesignWorkflow, self).end()
