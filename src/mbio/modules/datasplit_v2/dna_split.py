# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190124

import os
import re
import gevent
from collections import defaultdict
from biocluster.module import Module
from biocluster.file import download
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import MultiFileTransfer


class DnaSplitModule(Module):
    """
    动植物基因组多个文库进行二次拆分，拆分出每个文库对应的样本
    """
    def __init__(self, work_id):
        super(DnaSplitModule, self).__init__(work_id)
        options = [
            {"name": "lib_path", "type": "infile", "format": "datasplit.path"},  # list,存放文库文件夹对应的路径信息,第一列路径
            {"name": "library_info", "type": "infile", "format": "datasplit.dna_lib_info"},  # 文库信息，包含文库中样本酶的信息，文库类型等
            {"name": "ziplevel", "type": "string", "default": "6"},  # 压缩级别
            {"name": "combinatorial", "type": "string", "default": "2"},  # Use combinatorial barcode matching
            {"name": "mismatch", "type": "string", "default": "0"},  # mismatch
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("lib_path").is_set:
            raise OptionError("必须输入文库文件夹对应的路径信息")
        if not self.option("library_info").is_set:
            raise OptionError("必须输入文库的样本信息表")

    def get_info(self):
        """
        按照文库，将对应的文库序列及文库信息分开；
        文库类型为RAD及GBS的需要二次拆分，PE不需要跳过
        """
        pe_sample, lib_list = {}, []
        self.lib_type, self.lib_path, self.pe_info = {}, {}, {}
        sample_info = defaultdict(list)
        with open(self.option("library_info").prop["path"])as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                if re.search("PE", tmp[2]):
                    pe_sample[tmp[1]] = tmp[3]
                elif re.search("WGS", tmp[2]):
                    pe_sample[tmp[1]] = tmp[3]
                else:
                    if tmp[1] not in lib_list:
                        lib_list.append(tmp[1])
                    if re.search("RAD", tmp[2]):
                        self.lib_type[tmp[1]] = "RAD"
                    else:
                        self.lib_type[tmp[1]] = "GBS"
                    sample_info[tmp[1]].append(tmp)
        for lib in self.lib_type.keys():
            barcode_path = os.path.join(self.work_dir, lib + ".barcode.config")
            with open(barcode_path, "w+") as w:
                if self.lib_type[lib] == "RAD":
                    for tmp in sample_info[lib]:
                        # w.write(tmp[4] + "\t" + tmp[3] + "\n")  # enzyme1,sample
                        w.write(tmp[4] + "\t" + tmp[-1]+"--"+tmp[3] + "\n")  # enzyme1,sample
                else:
                    for tmp in sample_info[lib]:
                        # w.write(tmp[4] + "\t" + tmp[5] + "\t" + tmp[3] + "\n")  # enzyme1,enzyme2,sample
                        w.write(tmp[4] + "\t" + tmp[5] + "\t" + tmp[-1]+"--"+tmp[3] + "\n")  # enzyme1,enzyme2,sample
        with open(self.option("lib_path").prop["path"])as f:  # 文库对应序列路径
            for line in f:
                tmp = line.strip().split("\t")
                lib = tmp[0]
                if lib in pe_sample.keys():  # 如果是PE，单独以样本为键存路径
                    self.pe_info[lib + "--" + pe_sample[lib]] = tmp[1].strip().split(";")  # 换行符要去掉
                elif lib in lib_list:
                    self.lib_path[lib] = tmp[1].strip().split(";")  # 换行符要去掉
                else:
                    raise Exception("文库路径中的文库名{}在样本信息表中不存在，请核实！".format(lib))

    def run_axe_demux(self):
        self.tools = []
        for lib in self.lib_path:
            opts = {
                "fq1": self.lib_path[lib][0],
                "fq2": self.lib_path[lib][1],
                "lib_name": lib,
                "lib_type": self.lib_type[lib],
                "barcode_info": os.path.join(self.work_dir, lib + ".barcode.config"),
                "mismatch": self.option("mismatch"),
                "ziplevel": self.option("ziplevel"),
                "combinatorial": self.option("combinatorial")
            }
            self.axe_demux = self.add_tool("datasplit_v2.axe_demux")
            self.axe_demux.set_options(opts)
            self.tools.append(self.axe_demux)
        if len(self.tools) == 0:
            self.set_output()
        elif len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
            self.tools[0].run()
        else:
            self.on_rely(self.tools, self.set_output)
            for tool in self.tools:
                tool.run()

    def link(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        if old.startswith("s3:"):
            download(old, new)
        else:
            os.link(old, new)

    def download_pe_file(self):
        for sample in self.pe_info.keys():
            r1_path = self.output_dir + "/" + sample + ".R1.raw.fastq.gz"
            r2_path = self.output_dir + "/" + sample + ".R2.raw.fastq.gz"
            self.link(self.pe_info[sample][0], r1_path)
            self.link(self.pe_info[sample][1], r2_path)

    def set_output(self):
        self.logger.info("设置结果目录")
        if self.tools:
            for tool in self.tools:
                for f in os.listdir((tool.output_dir)):
                    old = os.path.join(tool.output_dir, f)
                    new = os.path.join(self.output_dir, f)
                    self.link(old, new)
            self.end()
        else:
            gevent.spawn_later(5, self.end)
        # self.end()

    def run(self):
        super(DnaSplitModule, self).run()
        self.get_info()
        # if len(self.pe_info) != 0:
        #     self.download_pe_file()
        self.run_axe_demux()

    def end(self):
        super(DnaSplitModule, self).end()
