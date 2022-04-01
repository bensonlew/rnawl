# -*- coding: utf-8 -*-
# __author__ = 'wentianliu'
# last modify 20181226

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import shutil


class TagStructureAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(TagStructureAgent, self).__init__(parent)
        options = [
            {"name": "populations_tag", "type": "infile", "format": "noref_wgs.populations_tag"},   # populations.tag
            {"name": "tag_id", "type": "string"},  # tag_id
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('populations_tag'):
            raise OptionError('必须输入:populations_tag', code="35501007")
        if not self.option('tag_id'):
            raise OptionError('必须输入:tag_id', code="35501008")
        if not self.option('main_id'):
            raise OptionError('必须输入:main_id', code="35501009")

    def set_resource(self):
        """
        运行所需资源
        vcftools一般不耗费内存资源，这里暂时不设置动态内存
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(TagStructureAgent, self).end()


class TagStructureTool(Tool):
    def __init__(self, config):
        super(TagStructureTool, self).__init__(config)

    def set_db(self):
        """
        导表
        """
        with open(self.option('populations_tag').prop['path'], "r")as fr:
            lines = fr.readlines()
            for i in range(len(lines)):
                if lines[i].strip("\n") == (">" + self.option('tag_id')):
                    ref = lines[i+1].strip()
                    snp = []
                    tmp = lines[i+2].strip("\n").split(" ")
                    for x in range(len(tmp)):
                        if tmp[x] == '*':
                            snp.append(x)
                            ref_list = list(ref)
                            ref_list[x] = "N"
                            ref = ''.join(ref_list)
                    sv_api = self.api.api("noref_wgs.tag_detail")
                    sv_api.add_sg_tag_detail(main_id=self.option('main_id'), tag_id=self.option('tag_id'), ref=ref, snp=snp)

    def run(self):
        super(TagStructureTool, self).run()
        self.set_db()
        self.end()
