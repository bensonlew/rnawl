# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 2018.06.29

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os


class GetMarkersAgent(Agent):
    """
    图谱标记详情中进行生成新的连锁分群
    逻辑根据图谱评估的结果，按照miss_ratio，signif以及useless_marker过滤得到有用的marker
    lasted modified by hongdong @ 20180629
    """
    def __init__(self, parent):
        super(GetMarkersAgent, self).__init__(parent)
        options = [
            {"name": "total_lg", "type": "infile", "format": "dna_gmap.lg"},
            {"name": "miss_ratio_start", "type": "float", "default": 30},
            {"name": "miss_ratio_end", "type": "float", "default": 100},
            {"name": "signif_start", "type": "float", "default": 0.05},
            {"name": "signif_end", "type": "float", "default": 1},
            {"name": "marker_info_path", "type": "string"},
            {"name": "useless_marker", "type": "string"},  # "marker01;marker02"
            {"name": "filter_lg", "type": "outfile", "format": "dna_gmap.lg"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("total_lg").is_set:
            raise OptionError("请设置total_lg文件", code="34800301")
        if not self.option("marker_info_path"):
            raise OptionError("请设置marker_info_path文件", code="34800302")

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(GetMarkersAgent, self).end()


class GetMarkersTool(Tool):
    def __init__(self, config):
        super(GetMarkersTool, self).__init__(config)

    def set_marker_file(self):
        useful_marker = self.get_useful_marker()
        self.logger.info("开始过滤markers")
        out_file = os.path.join(self.output_dir + "/total.lg")
        with open(self.option("total_lg").prop['path'], "r") as r, open(out_file, "w") as w:
            for line in r:
                if re.match(">.*", line):
                    w.write(line)
                else:
                    temp = line.strip().split("\t")
                    new_markfer = []
                    for i in range(0, len(temp)):
                        if temp[i] in useful_marker:
                            new_markfer.append(temp[i])
                    w.write("\t".join(new_markfer) + '\n')
        self.logger.info("开始过滤markers成功")
        self.option("filter_lg").set_path(out_file)

    def get_useful_marker(self):
        """
        这一步原本想放到package或者api中，因为marker文件比较大所以就写了这个tool
        total.marker.info, total.female.info
        :return:
        """
        markers = []
        uless_marker = self.useless_marker()
        for m in os.listdir(self.option("marker_info_path")):
            if re.match(r'total\..*\.info$', m):
                with open(os.path.join(self.option("marker_info_path"), m), "r") as r:
                    for line in r:
                        temp = line.strip().split('\t')
                        if re.match(r'^#.*', temp[0]):
                            pass
                        else:
                            if len(temp) == 10:
                                miss_ratio = (1 - float(temp[5]) / float(temp[4])) * 100
                                signif = float(temp[8])
                                if self.option("miss_ratio_end") >= miss_ratio >= self.option("miss_ratio_start") \
                                        and self.option("signif_end") >= signif >= self.option("signif_start") \
                                        and temp[0] not in uless_marker:
                                    markers.append(temp[0])
        self.logger.info("获取有用的markers成功")
        return markers

    def useless_marker(self):
        markers = []
        for m in self.option("useless_marker").split(";"):
            markers.append(m)
        return markers

    def run(self):
        super(GetMarkersTool, self).run()
        self.set_marker_file()
        self.end()
