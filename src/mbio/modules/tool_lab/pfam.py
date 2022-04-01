#-*- coding: utf-8 -*-

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class PfamModule(Module):
    def __init__(self, work_id):
        super(PfamModule, self).__init__(work_id)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},    #非冗余基因集的输入
            {"name": "evalue", "type": "float", "default": 1e-5},   #默认的evalue值
            {"name": "database", "type": "string"},  #pfam注释的数据库必须要加此参数--pfam
        ]
        self.add_option(options)    #检查option是否是list格式，体重每个option是否是字典类型

    def check_options(self):
        if not self.option('fasta').is_set:
            raise OptionError("必须设置参数fasta")
        return True

    def run_align(self):
        self.align_hmmscan = self.add_tool('align.mg_hmmscan')
        self.align_hmmscan.set_options({
            "faa_file": self.option('fasta'),
            "database": "pfam_v33.1",
        })
        self.align_hmmscan.on('end', self.run_anno)
        self.align_hmmscan.run()

    def run_anno(self):
        self.anno.set_options({
            'hmmscan_result': os.path.join(self.cat_out.output_dir, 'hmmscan_result'),
        })
        self.anno.on('end', self.set_output)
        self.anno.run()

    def set_output(self):
        self.option('pfam_result_dir', self.output_dir)
        anno_allfiles = os.listdir(self.anno.output_dir)
        stat_allfiles = os.listdir(self.pfam_anno_stat_tool.output_dir)
        out_oldfiles = [os.path.join(self.anno.output_dir, i) for i in anno_allfiles]
        stat_oldfiles = [os.path.join(self.pfam_anno_stat_tool.output_dir, i) for i in stat_allfiles]
        out_oldfiles.extend(stat_oldfiles)
        out_name = anno_allfiles + stat_allfiles
        output_newfiles = [os.path.join(self.output_dir, i) for i in out_name]
        for i in range(len(output_newfiles)):
            if os.path.exists(output_newfiles[i]):
                os.remove(output_newfiles[i])
            os.link(out_oldfiles[i], output_newfiles[i])
        self.end()

    def run(self):
        super(PfamModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(PfamModule, self).end()

