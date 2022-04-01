# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class KeggEnrichBubbleMultWorkflow(Workflow):
    """
    富集气泡图小工具
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(KeggEnrichBubbleMultWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table_file1", "type": "infile", "format": "small_rna.common"},
            {"name": "table_file2", "type": "infile", "format": "small_rna.common"},
            {"name": "table_file3", "type": "infile", "format": "small_rna.common"},
            {"name": "table_file4", "type": "infile", "format": "small_rna.common"},
            {"name": "table_file5", "type": "infile", "format": "small_rna.common"},
            {"name": "group_num", "type":"int"},
            {"name": "p_value", "type": "float"},
            {"name": "result_num", "type": "int"},
            {"name": "output_form", "type":"string"},
            {"name": "ids","type":"string"},
            {"name": "show_order","type":"string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.enrich_tool = self.add_tool("tool_lab.enrich_bubble")

    def check_option(self):
        """
        参数检查
        """
        if not self.option("group_num"):
            if not self.option("table_file1"):
                raise OptionError("必须设置作图数据")
            if not self.option("p_value"):
                raise OptionError("必须设置P值")
            if not self.option("result_num"):
                raise OptionError("result_num")
        else:
            for i in range(1,self.option("group_num")):
                if not self.option("table_file{}".format(i)) :
                    raise OptionError("缺少输入表格{}".format(i))
                elif not self.option("table_file{}".format(i)).prop["path"].endwith("xls") and not self.option("table_file{}".format(i)).prop["path"].endwith("xlsx"):
                    raise OptionError("输入文件格式有误，为Excel文件")
            if self.option("output_form") not in ["union","intersection","ids"]:
                raise OptionError("需要提供输出形式")
            if self.option("output_form") == "interation":
                if not self.option("show_order"):
                    raise OptionError("当结果以交集展示时，需要指定显示方式")
    
    def run_tools(self):
        if not self.option("group_num"):
            options = {
                "table_file1":self.option("table_file1"),
                "p_value":self.option("p_value"),
                "result_num":self.option("result_num"),
            }
            self.major = os.path.basename(self.option("table_file1").prop["path"]).split(".")[0]
        else:
            options = {
                "output_form":self.option("output_form"),
                "group_num":self.option("group_num"),
                "p_value":self.option("p_value"),
            }   
            self.major = ""
            for i in range(1,self.option("group_num")+1):
                options["table_file{}".format(i)] = self.option("table_file{}".format(i))
            if self.option("output_form") == "intersection":
                options["show_order"] = self.option("show_order")
                self.major = os.path.basename(self.option(self.option("show_order")).prop["path"]).split(".")[0]
            else:
                self.major = os.path.basename(self.option("table_file1").prop["path"]).split(".")[0]
            if self.option("output_form") == "ids":
                options["ids"]= self.option("ids")
            else:
                options["result_num"] = self.option("result_num")
        self.enrich_tool.set_options(options)
        self.enrich_tool.on("end", self.set_output)
        self.enrich_tool.run()

    def set_output(self):
        self.logger.info('strat set_output as {}'.format(
            self.__class__.__name__))
        if not self.option("group_num"):
            input_name = os.path.basename(self.option("table_file1").prop["path"]).split(".")[0]
            os.link(os.path.join(self.enrich_tool.output_dir,"{}_bubble.xls".format(input_name))
                ,os.path.join(self.output_dir,"{}_bubble.xls".format(input_name)))
        else:    
            for i in range(1,self.option("group_num")+1):
                table_name = "table_file{}".format(i)
                input_name = os.path.basename(self.option(table_name).prop["path"]).split(".")[0]
                os.link(os.path.join(self.enrich_tool.output_dir,"{}_bubble.xls".format(input_name))
                        ,os.path.join(self.output_dir,"{}_bubble.xls".format(input_name)))
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_bubble = self.api.api("tool_lab.enrich_bubble")
        api_bubble.add_detail(self.option("main_id"),self.output_dir,self.major)
        api_bubble.add_tooltip(self.option("main_id"),os.path.join(self.enrich_tool.output_dir,"term.txt"),"kegg")    
        self.logger.info("导表结束")
        self.end()

    def run(self):
        """
        运行
        """
        self.run_tools()
        super(KeggEnrichBubbleMultWorkflow, self).run()

    def end(self):
        super(KeggEnrichBubbleMultWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    # data = {
    #     'name': 'test_enrichbubble',
    #     'id': 'enrichbubble_' +  str(random.randint(1, 10000)),
    #     'type': 'workflow',
    #     'options': {
    #         "table_file1": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/go_pvalue.xls",
    #         "table_file2": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/go_pvalue1.xls",
    #         "group_num":2,
    #         "output_form":"intersection",
    #         "p_value" : 0.05, 
    #         "result_num":20,
    #         "show_order":"table_file1",
    #         "main_id" : "5e9e6a6017b2bf2049a81be3"
    #     }
    # }
    data = {
        'name': 'test_enrichbubble',
        'id': 'enrichbubble_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        'options': {
        "table_file1": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/GO-1.xls",
        "table_file2": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/GO-2.xls",
        "table_file3": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/GO-3.xls",
        "output_form":"ids",
        "p_value" : 0.05, 
        "result_num":20,
        "show_order":"table_file1",
        "ids":"map00942;map00902;map00592;map00945;map04075;map00941;map00900;map00520;map00563;map00940;map00062;map00908;map03022;map00920;map04130;map04016;map00944;map04070;",
        "main_id" : "5e9e6a6017b2bf2049a81be1"
        }
    }
    wsheet = Sheet(data=data)
    wf = KeggEnrichBubbleMultWorkflow(wsheet)
    wf.run()

            
