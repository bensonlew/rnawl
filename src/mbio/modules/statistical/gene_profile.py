# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:2018.0228

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class GeneProfileModule(Module):
    def __init__(self, parent):
        super(GeneProfileModule, self).__init__(parent)
        options = [
            {"name": "map_dir", "type": "string", "default": ""},  # map结果
            {"name": "fafile", "type": "infile", "format": "sequence.fasta"},  # 非冗余基因集fasta文件
            {"name": "insertsize", "type": "infile", "format": "sample.insertsize_table"},  # 插入片段文件
            {"name": "rpkm_abundance", "type": "outfile", "format": "sequence.profile_table"},  # RPKM丰度
            {"name": "reads_abundance", "type": "outfile", "format": "sequence.profile_table"},  # reads丰度
        ]
        self.add_option(options)
        self.merge_profile = self.add_tool("statistical.merge_profile")
        self.soap_info = self.add_tool("statistical.soap_info")
        self.sample_profile = []
        self.step.add_steps("pre", "merge")

    def check_options(self):
        if not self.option("fafile").is_set:
            raise OptionError("必须提供非冗余基因集", code="24100101")
        if not self.option("map_dir"):
            raise OptionError("必须提供map结果", code="24100102")
        if not self.option("insertsize").is_set:
            raise OptionError("必须插入片段文件", code="24100103")

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

    def pre_run(self):
        self.soap_info.set_options({
            "insertsize": self.option("insertsize"),
            "map_dir": self.option("map_dir"),
        })
        self.soap_info.on("start", self.set_step, {"start": self.step.pre})
        self.soap_info.on("end", self.set_step, {"end": self.step.pre})
        self.soap_info.run()

    def profile_run(self):
        f = open(self.soap_info.option("soap_info").prop['path'],'r')
        lines = f.readlines()
        n = 1
        sample = []
        for line in lines[1:]:
            line = line.strip()
            line_info = line.split("\t")
            profile = self.add_tool("statistical.sample_unigene_profile")
            self.step.add_steps('profile_{}'.format(n))
            opts = {
                "map_result": line_info[2],
                "fafile": self.option("fafile"),
                "insertsize": line_info[1],
                "sample": line_info[0],
                "out_dir": self.work_dir
            }
            sample.append(line_info[0])
            profile.set_options(opts)
            step = getattr(self.step, 'profile_{}'.format(n))
            step.start()
            n += 1
            self.sample_profile.append(profile)
        f.close()
        self.sample = ','.join(sample)
        if len(self.sample_profile) == 1:
            self.sample_profile[0].on('end', self.merge_run)
        else:
            self.on_rely(self.sample_profile, self.merge_run)
        for tool in self.sample_profile:
            tool.run()

    def merge_run(self):
        self.merge_profile.set_options({
            "profile_dir": self.work_dir ,
            "fafile": self.option("fafile"),
            "samples": self.sample,
        })
        self.merge_profile.on("start", self.set_step, {"start": self.step.merge})
        self.merge_profile.on("end", self.set_step, {"end": self.step.merge})
        self.merge_profile.run()

    def set_output(self):
        # self.linkdir(self.merge.output_dir, os.path.join(self.output_dir, 'uniGeneset'))
        # self.linkdir(self.length_tool.output_dir, os.path.join(self.output_dir, 'length_distribute'))
        self.linkdir(self.merge_profile.output_dir, "")
        self.option('reads_abundance').set_path(self.output_dir + '/reads_number.xls')
        self.option('rpkm_abundance').set_path(self.output_dir + '/RPKM.xls')
        self.end()

    def run(self):
        super(GeneProfileModule, self).run()
        self.soap_info.on("end", self.profile_run)
        self.merge_profile.on("end", self.set_output)
        self.pre_run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["reads_number.xls", "xls", "reads丰度表"],
            ["reads_percent.xls", "xls", "reads相对丰度表"],
            ["base_percent.xls", "xls", "reads数除以基因长度丰度表"],
            ["base_number.xls", "xls", "reads数除以基因长度相对丰度表"],
            ["RPKM.xls", "xls", "RPKM丰度表"],
            ["RPKM_percent.xls", "xls", "RPKM相对丰度表"]
        ])
        super(GeneProfileModule, self).end()

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
