# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180208

from biocluster.module import Module
import os
import shutil
import re
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.file_sample import FileSampleFile


class SplitMapModule(Module):
    def __init__(self, work_id):
        super(SplitMapModule, self).__init__(work_id)
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
        self.step.add_steps('split')
        self.split = self.add_tool("sequence.cdhit_split_fasta")
        self.unigene_profile = self.add_module("statistical.gene_profile")
        self.map_list = []

    def check_options(self):
        if not self.option("repeat") in [0, 1, 2]:
            raise OptionError("repeat必须为0,1,或2", code="21100601")
        if not self.option("mode") in [0, 1, 2, 4]:
            raise OptionError("repeat必须为0,1,2,或4", code="21100602")
        if not 0 < self.option("seed") <= 256:
            raise OptionError('seed参数必须设置在1-256之间:%s', variables=(self.option('seed')), code="21100603")
        if not 0 < self.option("identity") <= 1:
            raise OptionError("identity必须在0，1之间", code="21100604")
        if not self.option("fafile").is_set:
            raise OptionError("必须提供非冗余基因集", code="21100605")
        if not self.option("QC_dir").is_set:
            raise OptionError("必须提供质控后的fq文件夹", code="21100606")
        if self.option("QC_dir").is_set:
            list_path = os.path.join(self.option("QC_dir").prop["path"], "list.txt")
            if not os.path.exists(list_path):
                OptionError("缺少list文件", code="21100607")
            row_num = len(open(list_path, "r").readline().split())
            if row_num != 3:
                raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code="21100608")
        if not self.option("insertsize").is_set:
            raise OptionError("必须提供插入片段文件", code="21100609")

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

    def div_fasta(self):
        self.number = os.path.getsize(self.option("fafile").prop['path']) / 4000000000 + 1
        # self.number = 3
        self.split.set_options({"gene_tmp_fa": self.option("fafile"),
                                "ou_dir": self.work_dir + '/fasta_div',
                                "number": self.number,
                                "pre": "geneset.div",
                                "order": 0
                                })
        self.logger.info(self.split)
        self.split.on("start", self.set_step, {'start': self.step.split})
        self.split.on("end", self.set_step, {'end': self.step.split})
        self.split.on("end", self.map)
        self.split.run()

    def map(self):
        n = 1
        for i in range(0, self.number):
            map = self.add_module("align.map_geneset")
            self.step.add_steps('map_{}'.format(n))
            fa = self.work_dir + '/' + self.split.option("pre") + '-' + str(i) + '.fa'
            if os.path.exists(fa):
                os.remove(fa)
            os.link(self.split.option("ou_dir") + '/' + self.split.option("pre") + '-' + str(i),fa)
            map.set_options({
                "fafile": fa,
                "insertsize": self.option("insertsize"),
                "QC_dir": self.option("QC_dir"),
                "seed": self.option("seed"),
                "mode": self.option("mode"),
                "processors": self.option("processors"),
                "mismatch": self.option("mismatch"),
                "repeat": self.option("repeat"),
                "identity": self.option("identity")
            })
            step = getattr(self.step, 'map_{}'.format(n))
            step.start()
            map.on("end", self.finish_update, "map_{}".format(n))
            n += 1
            self.map_list.append(map)
        if len(self.map_list) == 1:
            self.map_list[0].on("end", self.profile)
        else:
            self.on_rely(self.map_list, self.profile)
        for module in self.map_list:
            module.run()

    def profile(self):
        dir_list = []
        for module in self.map_list:
            dir_list.append(module.output_dir + '/map_dir')
        map_dir = ','.join(dir_list)
        opts = {
            "fafile": self.option("fafile"),
            "insertsize": self.option("insertsize"),
            "map_dir": map_dir
        }
        self.unigene_profile.set_options(opts)
        self.unigene_profile.on("end", self.set_output)
        self.unigene_profile.run()

    def set_output(self):
        self.linkdir(self.unigene_profile.output_dir, "gene_profile")
        self.option('reads_abundance', self.unigene_profile.option("reads_abundance"))
        self.option('rpkm_abundance', self.unigene_profile.option("rpkm_abundance"))
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
        super(SplitMapModule, self).run()
        self.div_fasta()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(SplitMapModule, self).end()
