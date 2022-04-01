# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:2017.08.28

from biocluster.module import Module
import os
import shutil
import re
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.file_sample import FileSampleFile


class MapGenesetModule(Module):
    def __init__(self, work_id):
        super(MapGenesetModule, self).__init__(work_id)
        options = [
            {"name": "fafile", "type": "infile", "format": "sequence.fasta"},  # 非冗余基因集
            {"name": "insertsize", "type": "infile", "format": "sample.insertsize_table"},  # 插入片段文件
            {"name": "QC_dir", "type": "infile", "format": "sequence.fastq_dir"},  # qc后fastq文件夹
            {"name": "reads_abundance", "type": "outfile", "format": "sequence.profile_table"},  # reads_abundance
            {"name": "rpkm_abundance", "type": "outfile", "format": "sequence.profile_table"},  # rpkm_abundance
            {"name": "seed", "type": "int", "default": 35},
            # align the initial n bps as a seed means whole lengths of read
            {"name": "mode", "type": "int", "default": 4},
            # match mode for each read or the seed part of read, which shouldn't contain more than 2 mismatches: 0 for exact mathc only; 1 for 1 mismatch; 2 for 2 mismatch; 4 for find the best hits
            {"name": "processors", "type": "int", "default": 6},
            {"name": "mismatch", "type": "int", "default": 20},  # maximum number of mismatches allowed on a read
            {"name": "repeat", "type": "int", "default": 1},  # how to report repeat hits, 0=none, 1=random one, 2=all
            {"name": "identity", "type": "float", "default": 0.95}  # identity
        ]
        self.add_option(options)
        self.build = self.add_tool("align.bwt_builder")
        self.aligner = []
        self.step.add_steps('bwt_build')
        self.samples = {}

    def check_options(self):
        if not self.option("repeat") in [0, 1, 2]:
            raise OptionError("repeat必须为0,1,或2", code="21100401")
        if not self.option("mode") in [0, 1, 2, 4]:
            raise OptionError("repeat必须为0,1,2,或4", code="21100402")
        if not 0 < self.option("seed") <= 256:
            raise OptionError('seed参数必须设置在1-256之间:%s', variables=(self.option('seed')), code="21100403")
        if not 0 < self.option("identity") <= 1:
            raise OptionError("identity必须在0，1之间", code="21100404")
        if not self.option("fafile").is_set:
            raise OptionError("必须提供非冗余基因集", code="21100405")
        if not self.option("QC_dir").is_set:
            raise OptionError("必须提供质控后的fq文件夹", code="21100406")
        if self.option("QC_dir").is_set:
            list_path = os.path.join(self.option("QC_dir").prop["path"], "list.txt")
            if not os.path.exists(list_path):
                OptionError("缺少list文件", code="21100407")
            row_num = len(open(list_path, "r").readline().split())
            if row_num != 3:
                raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code="21100408")
        if not self.option("insertsize").is_set:
            raise OptionError("必须提供插入片段文件", code="21100409")

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def build_index(self):
        self.build.set_options({"fafile": self.option("fafile")})
        self.logger.info("build index")
        self.build.on("start", self.set_step, {'start': self.step.bwt_build})
        self.build.on("end", self.set_step, {'end': self.step.bwt_build})
        self.build.on("end", self.soap_align)
        self.build.run()

    def get_column(self, filename, splitregex='\t'):
        with open(filename, 'rt') as handle:
            for ln in handle:
                items = re.split(splitregex, ln)
                yield items[0], items[1]

    def get_list(self):
        list_path = os.path.join(self.option("QC_dir").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        return samples

    def soap_align(self):
        self.samples = self.get_list()
        n = 1
        for x, y in self.get_column(self.option("insertsize").prop['path']):
            align = self.add_tool("align.soap_aligner")
            self.step.add_steps('align_{}'.format(n))
            options = {
                "sample": x,
                "insertSize": y,
                "index": self.build.option("build_dir"),
                "fq_r": os.path.join(self.option("QC_dir").prop["path"], self.samples[x]["r"]),
                "fq_l": os.path.join(self.option("QC_dir").prop["path"], self.samples[x]["l"]),
                #"fq_s": os.path.join(self.option("QC_dir").prop["path"], self.samples[x]["s"]),
                "repeat": self.option("repeat"),
                "seed": self.option("seed"),
                "mode": self.option("mode"),
                "processors": self.option("processors"),
                "mismatch": self.option("mismatch"),
                "identity": self.option("identity"),
                "fafile": self.option("fafile")
            }
            ### 兼容无s端的
            if self.samples[x].has_key("s"):
                options["fq_s"] = os.path.join(self.option("QC_dir").prop["path"], self.samples[x]["s"])
            align.set_options(options)
            step = getattr(self.step, 'align_{}'.format(n))
            step.start()
            align.on("end", self.finish_update, "align_{}".format(n))
            n += 1
            self.aligner.append(align)
            self.logger.info(x + str(y))
            #self.logger.info(self.samples[x]["r"] + self.samples[x]["l"] + self.samples[x]["s"])
        if len(self.aligner) == 1:
            self.aligner[0].on("end", self.set_output)
        else:
            self.on_rely(self.aligner, self.set_output)
        for tool in self.aligner:
            tool.run()

    # def profile(self):
    #     if not os.path.exists(os.path.join(self.output_dir, "map_dir")):
    #         os.mkdir(os.path.join(self.output_dir, "map_dir"))
    #     for tool in self.aligner:
    #         out_files = os.listdir(tool.option("map_dir").prop['path'])
    #         for f in out_files:
    #             f_path = os.path.join(tool.option("map_dir").prop['path'], f)
    #             self.linkdir(f_path, "map_dir/" + f)
    #     opts = {
    #         "fafile": self.option("fafile"),
    #         "insertsize": self.option("insertsize"),
    #         "map_dir": os.path.join(self.output_dir, 'map_dir')
    #     }
    #     self.unigene_profile.set_options(opts)
    #     self.unigene_profile.on("end", self.set_output)
    #     self.unigene_profile.run()
    #
    def set_output(self):
        if not os.path.exists(os.path.join(self.output_dir, "map_dir")):
            os.mkdir(os.path.join(self.output_dir, "map_dir"))
        for tool in self.aligner:
            out_files = os.listdir(tool.option("map_dir").prop['path'])
            for f in out_files:
                f_path = os.path.join(tool.option("map_dir").prop['path'], f)
                self.linkdir(f_path, "map_dir/" + f)
        self.end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(MapGenesetModule, self).run()
        self.build_index()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(MapGenesetModule, self).end()
