# -*- coding: utf-8 -*-
# __author__ = 'zhengyuan'


from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class HmdbBlastModule(Module):
    def __init__(self, work_id):
        super(HmdbBlastModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "evalue", "type": "float", "default": 1e-5}
        ]
        self.add_option(options)
        # self.step.add_steps("blast", "annotation")
        self.blast = self.add_tool("align.blast")
        self.annotation = self.add_tool("taxon.cgc_taxon")

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须输入fasta文件")

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_blast(self):
        self.blast.set_options({
            "query": self.option("query"),
            "query_type": "nucl",
            "database": "cgc",
            "outfmt": 6,
            "blast": "blastn",
            "evalue": self.option("evalue")
        })
        self.blast.on("end", self.run_annotation)
        self.blast.run()

    def run_annotation(self):
        infile = self.blast.option("outtable")
        self.annotation.set_options({
            "blastoutput": infile
        })
        self.output_dir = self.annotation.output_dir
        self.annotation.run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种组成分析结果目录"],
            ['annotation_result.xls', 'xls', 'gene对应注释文件']
        ])
        super(HmdbBlastModule, self).end()


    def run(self):
        super(HmdbBlastModule, self).run()
        self.run_blast()
        self.annotation.on("end", self.end)
