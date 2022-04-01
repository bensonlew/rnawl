# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20200609

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
import gevent.subprocess as subprocess
import json
from collections import OrderedDict
from biocluster.config import Config


class DiffGenesetPrepareModule(Module):
    """
    该Module用于基因融合分析，默认使用方法
    """
    def __init__(self, work_id):
        super(DiffGenesetPrepareModule, self).__init__(work_id)
        options = [
            # {"name": "geneset_names", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "geneset_names", "type": "string", "default": None},
            {"name": "diff_id", "type": "string"},
            {"name": "task_id", "type": "string", "default": None},
            {"name": "level", "type": "string", "default": None},
            {"name": "species", "type": "string", "default": "Homo_sapiens"}
        ]
        self.add_option(options)
        self.geneset_prepare_tools=[]
        self.genest_infos = OrderedDict()
        self.common = self.add_tool("medical_transcriptome.diff_geneset.annot_prepare")
        self.total_dict = OrderedDict()


    def check_options(self):
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def get_genesets_info(self):
        genesets_info = json.loads(self.option("geneset_names"))
        for gi in genesets_info:
            regulate = gi
            for geneset in genesets_info[regulate]:
                geneset_name = geneset
                geneset_details = geneset_name.split("_")
                ctrl,test ="_".join(geneset_details[:-2]).split("_vs_")
                compare = ctrl+"|"+test
                self.genest_infos[geneset_name]={
                    "geneset_name": geneset_name,
                    "compare" : compare,
                    "regulate":regulate
                }
        self.prepare_common_file()

    def prepare_common_file(self):
        self.logger.info("开始准备差异数据挖掘公共文件")
        opts = {
            "task_id": self.option("task_id"),
            "level": self.option("level"),
            "species": self.option("species"),
        }
        self.common.set_options(opts)
        self.common.on("end", self.prepare_geneset_infos)
        self.common.run()

    def prepare_geneset_infos(self):
        for key in self.genest_infos:
            diff_geneset = self.add_tool("medical_transcriptome.diff_geneset.diff_geneset_prepare")
            opts = {
                "diff_id": self.option("diff_id"),
                "task_id": self.option("task_id"),
                "compare": self.genest_infos[key]["compare"],
                "regulate": self.genest_infos[key]["regulate"],
                "geneset_name": self.genest_infos[key]["geneset_name"],
                "level": self.option("level"),
                "species": self.option("species"),
            }
            diff_geneset.set_options(opts)
            self.geneset_prepare_tools.append(diff_geneset)
        if self.geneset_prepare_tools:
            if len(self.geneset_prepare_tools) > 1:
                self.on_rely(self.geneset_prepare_tools, self.set_output)
            elif len(self.geneset_prepare_tools) == 1:
                self.geneset_prepare_tools[0].on('end', self.set_output)
        else:
            self.set_error("geneset_prepare_tools列表为空！")
        for tool in self.geneset_prepare_tools:
            gevent.sleep(1)
            tool.run()

    def set_output(self):
        common_file_path = self.common.option("common_file_json").prop["path"]
        common_file_dict = json.load(open(self.common.option("common_file_json").prop["path"]))
        self.total_dict["common_file"] = common_file_dict
        self.total_dict["genesets"] =OrderedDict()
        for geneset_prepare_tools in self.geneset_prepare_tools:
            geneset_json = geneset_prepare_tools.option("geneset_json").prop["path"]
            geneset_dict = json.load(open(geneset_json))
            geneset_name = geneset_dict["geneset_name"]
            self.total_dict["genesets"][geneset_name] = geneset_dict
        with open(os.path.join(self.output_dir,"prepare_json"),"w") as f:
            json.dump(self.total_dict,f,indent=2)
        self.end()


    def run(self):
        super(DiffGenesetPrepareModule, self).run()
        self.logger.info("开始运行差异基因集数据挖掘数据准备")
        self.logger.info("首先解析基因集列表文件")
        self.get_genesets_info()

    def end(self):
        super(DiffGenesetPrepareModule, self).end()



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "geneset_prepare" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "medical_transcriptome.diff_geneset.diff_geneset_prepare",
            "instant": False,
            "options": dict(
                geneset_names='{"down": ["H1_vs_H3_down_12", "S1_vs_S3_down_12", "H1_vs_H2_down_12", "S1_vs_S2_down_12"], "all": ["H1_vs_H3_all_12", "S1_vs_S3_all_12", "H1_vs_H2_all_12", "S1_vs_S2_all_12"], "up": ["H1_vs_H3_up_12", "S1_vs_S3_up_12", "H1_vs_H2_up_12", "S1_vs_S2_up_12"]}',
                diff_id="5f45c6d117b2bf78d9c9c16d",
                task_id ="medical_transcriptome",
                level = "G",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()