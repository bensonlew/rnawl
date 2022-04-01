# -*- coding: utf-8 -*-
# __author__ = 'zhengyuan'


from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class HmdbMapModule(Module):
    def __init__(self, work_id):
        super(HmdbMapModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "method", "type": "string", "default": "global"}
        ]
        self.add_option(options)
        # self.step.add_steps("map", "annotation")
        self.map = self.add_tool("hmdb.map")
        self.annotation = self.add_tool("taxon.cgc_taxon")

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须输入fasta文件")
        if self.option("method") not in ["global", "local"]:
            raise OptioError("比对方法必须为global或local")

    def run_map(self):
        self.map.set_options({
            "fa_query": self.option("query"),
            "method": self.option("method")
        })
        self.map.on("end", self.run_annotation)
        self.map.run()

    def run_annotation(self):
        infile = self.map.option("abund")
        self.annotation.set_options({
            "genetable": infile
        })
        self.output_dir = self.annotation.output_dir
        self.annotation.on("end", self.end)
        self.annotation.run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种组成分析结果目录"],
            ['annotation_result.xls', 'xls', 'gene对应注释文件']
        ])
        super(HmdbMapModule, self).end()


    def run(self):
        super(HmdbMapModule, self).run()
        self.run_map()
