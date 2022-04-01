# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class HumannAnnoModule(Module):

    def __init__(self, work_id):
        super(HumannAnnoModule, self).__init__(work_id)
        self.step.add_steps('align', 'humann_anno')
        options = [
            {"name": "reads_fasta", "type": "infile",
             "format": "sequence.fasta_dir"},
            {"name": "lines", "type": "int", "default": 2000},
            {"name": "align_software", 'type': "string", "default": "diamond"},  # 或者是blast
            {"name": "evalue", "type": "float", "default": 1e-5}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("reads_fasta"):
            raise OptionError("必须reads序列fasta文件夹")
        if self.option("align_software") not in ['blast', 'diamond']:
            raise OptionError("比对软件必须是blast和diamond之一: {}".format(self.option("align_software")))
        return True

    def run(self):
        super(HumannAnnoModule, self).run()
        self.step.align.start()
        self.step.update()
        self.run_align()

    def run_align(self):
        """"""
        if self.option("align_software") == "diamond":
            align = "align.diamond"
        else:
            align = "align.blast"
        self.align_all = []
        self.logger.warn("此处需要进行fasta文件检查，可能很慢，造成问题")
        for fasta in self.option("reads_fasta").prop['fasta_fullname']:
            align_module = self.add_module(align)
            align_module.set_options(
                {
                    "query": fasta,
                    "lines": self.option("lines"),
                    "query_type": "nucl",
                    "database": "kegg_microbe",
                    "outfmt": 5,
                    "blast": "blastx",
                    "evalue": self.option("evalue")
                }
            )
            self.align_all.append(align_module)
        if len(self.align_all) == 1:
            self.align_all[0].on("end", self.run_humann)
        elif len(self.align_all) > 1:
            self.on_rely(self.align_all, self.run_humann)
        else:
            raise Exception("fasta文件夹对象没有任何fasta文件")
        for i in self.align_all:
            i.run()

    def run_humann(self):
        self.step.align.finish()
        self.step.humann_anno.start()
        self.step.update()
        align_out = os.path.join(self.work_dir, "output")
        if not os.path.exists(align_out):
            os.makedirs(align_out)
        for i in self.align_all:
            file_name = os.listdir(i.work_dir + "/blast_tmp")
            name = os.path.basename(i.option('query').prop['path']).split('.')[0]
            for f in file_name:
                os.link(i.work_dir + "/blast_tmp/" + f, align_out + "/" + name + '-' + f)
            # os.link(i.option("outxml").path, align_out + '/' + os.path.basename(i.option("query").path) + '.xml')
        tool = self.add_tool("annotation.kegg.human_anno")
        tool.set_options({"blastout": align_out})
        tool.on('end', self.step_end)
        tool.run()

    def step_end(self):
        self.step.humann_anno.finish()
        self.step.update()
        self.end()

    def end(self):
        repaths = [
            [".", "", "huamnn注释结果文件夹"],
        ]
        regexps = [
            [r'.*\.xml$', 'xml', '比对结果文件'],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(HumannAnnoModule, self).end()
