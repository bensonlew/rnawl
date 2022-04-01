# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'binbin.zhao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.itraq_and_tmt.go_graph import Terms
import os


class GoBarAgent(Agent):
    """
    CorHeatmapAgent:用于生成之间的correlation
    """
    def __init__(self, parent):
        super(GoBarAgent, self).__init__(parent)
        options = [
            {"name": "go_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "main_id", "type": "string"},

        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("go_table"):
            raise OptionError("必须输入snp_table文件")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        """
        计算结束
        """
        super(GoBarAgent, self).end()


class GoBarTool(Tool):
    """
    用于相关性分析
    """
    def __init__(self, config):
        super(GoBarTool, self).__init__(config)

    def calculate_go(self):
        """
        开始计算绘图
        """
        a = Terms(obo_fp=self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome/obo_file/go-basic.obo")
        self.logger.info("开始计算go")
        self.logger.info(a.get_term("GO:0000016"))
        outfile = self.output_dir + "/go_table"
        with open(outfile, "w") as w, open(self.option("go_table").prop["path"]) as f:
            lines = f.readlines()
            go_num = {}
            list_id = []
            for line in lines:
                go_id = line.strip().split('\t')[1:]
                if ';' in go_id[0] or ',' in go_id[0] or '、' in go_id[0] or ' ' in go_id[0]:
                    self.logger.info("GO号必须以制表符分隔，请修改！")
                    self.set_error("GO号必须以制表符分隔，请修改！")
                    break
                else:
                    for id in go_id:
                        if id not in go_num.keys():
                            go_num[id] = 1
                            list_id.append(id)
                        else:
                            go_num[id] += 1
            if len(list_id) > 0:
                w.write("go_id\tdescription\tterm_type\tnumber\tpercent\n")
                for go_id in list_id:
                    try:
                        a.get_term(go_id).name
                    except:
                        pass
                    else:
                        w.write(go_id + "\t" + a.get_term(go_id).name + "\t" + a.get_term(go_id).namespace + "\t"
                            + str(go_num[go_id]) + "\t" +str(go_num[go_id] * 100/len(list_id)) + "%" + "\n")
                    # w.write(go_id + "\t" + a.get_term(go_id).name + "\t" + a.get_term(go_id).namespace + "\t"
                    #         + str(go_num[go_id]) + "\t" +str(go_num[go_id] * 100/len(list_id)) + "%" + "\n")

    def run(self):
        """
        运行
        """
        super(GoBarTool, self).run()
        self.calculate_go() # 对输入数据进行处理，得到new_corr_table.csv
        self.set_db()
        self.end()

    def set_db(self):
        self.logger.info("tool中导表结束")
        api_gobar = self.api.api('tool_lab.go_bar')
        api_gobar.add_sg_gobar(self.option('main_id'), self.output_dir + "/go_table")
        self.logger.info("tool中导表结束")