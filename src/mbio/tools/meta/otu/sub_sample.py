# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import os
import re
import subprocess
import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.files.meta.otu.otu_table import OtuTableFile


class SubSampleAgent(Agent):
    """
    version 1.0
    author: xuting
    last_modify: 2015.11.19
    需要mothur 版本1.30
    需要shared2otu.pl
    需要otu2shared.pl
    """
    def __init__(self, parent):
        super(SubSampleAgent, self).__init__(parent)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table,meta.otu.tax_summary_dir"},  # 输入的OTU文件
            {"name": "out_otu_table", "type": "outfile", "format": "meta.otu.otu_table"},  # 输出的OTU文件
            {"name": "level", "type": "string", "default": "otu"},
            {"name": "size", "type": "string", "default": "min"}
        ]
        self.add_option(options)
        self.step.add_steps("sub_sample")
        self.on('start', self.start_sub_sample)
        self.on('end', self.end_sub_sample)

    def start_sub_sample(self):
        self.step.sub_sample.start()
        self.step.update()

    def end_sub_sample(self):
        self.step.sub_sample.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option("in_otu_table").is_set:
            raise OptionError("输入的OTU文件不能为空", code="32705501")
        if self.option("level") not in ['otu', 'domain', 'kindom', 'phylum', 'class',
                                        'order', 'family', 'genus', 'species']:
            raise OptionError("请选择正确的分类水平", code="32705502")
        return True

    def set_resource(self):
        """
        设置所需的资源
        """
        self._cpu = 2
        self._memory = "20G"
        self._memory_increase_step = 20

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ['\.subsample\.', 'meta.otu.otu_table', "抽平后的otu表格"]
        ])
        super(SubSampleAgent, self).end()


class SubSampleTool(Tool):
    def __init__(self, config):
        super(SubSampleTool, self).__init__(config)
        self.mothur_path = "bioinfo/meta/mothur-1.30/mothur.1.30"
        self.shared2otu_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/meta/scripts/shared2otu.pl")
        self.otu_tax = dict()
        self.has_tax = False
        self.basename = ""

    def sub_sample(self):
        """
        运行mothur的subsample，进行抽平
        """
        new2old = dict()  # mothur进行抽平的时候会对OTU进行重命名， 我们需要找回原有的OTU名字
        if self.option("in_otu_table").format == "meta.otu.tax_summary_dir":
            otu_table = os.path.basename(self.option("otu_table").get_table(self.option("level")))
            self.basename = otu_table
        else:
            otu_table = os.path.basename(self.option("in_otu_table").prop["path"])
            self.basename = otu_table
            self.option("in_otu_table").get_info()
            if self.option("in_otu_table").prop["metadata"] == "taxonomy":
                no_taxnomy_otu_table_path = os.path.join(self.work_dir, otu_table + ".no_tax")
                tax_dic = self.option("in_otu_table").del_taxonomy(self.option("in_otu_table").prop['path'], no_taxnomy_otu_table_path)
                self.has_tax = True
        shared_path = os.path.join(self.work_dir, otu_table + ".shared")
        mothur_dir = os.path.join(self.work_dir, "mothur")
        if not os.path.exists(mothur_dir):
            os.mkdir(mothur_dir)
        my_table = OtuTableFile()
        if self.option("in_otu_table").format == "meta.otu.tax_summary_dir":
            my_table.set_path(self.option("in_otu_table").get_table(self.option("level")))
        else:
            if self.has_tax is True:
                my_table.set_path(no_taxnomy_otu_table_path)
            else:
                my_table.set_path(self.option("in_otu_table").prop['path'])
        my_table.get_info()
        (min_sample, min_num) = my_table.get_min_sample_num()
        if self.option("size") != "min":
            if min_num < int(self.option("size")):
                self.set_error("自定义抽平序列数目大于最小的样本序列数目！", code="32705501")
                raise Exception("自定义抽平序列数目大于最小的样本序列数目！")

        with open(my_table.prop["path"], "rb") as r:
            r.next()
            c = 1
            for line in r:
                line = line.split("\t")
                new2old["OTU" + str(c)] = line[0]
                c += 1

        my_table.convert_to_shared(shared_path)
        if self.option("size") in [0, "0", "min"]:
            cmd = self.mothur_path + " \"#set.dir(output=" + mothur_dir\
                + ");sub.sample(shared=" + shared_path + ")\""
        else:
            cmd = self.mothur_path + " \"#set.dir(output=" + mothur_dir\
                + ");sub.sample(shared=" + shared_path + ", size=" + str(self.option("size")) + ")\""
        sub_sample_cmd = self.add_command("sub_sample_cmd", cmd)
        self.logger.info("开始运行sub.sample")
        sub_sample_cmd.run()
        self.wait(sub_sample_cmd)
        if sub_sample_cmd.return_code == 0:
            self.logger.info("运行sub.sample完成")
        else:
            self.set_error("运行sub.sample出错", code="32705502")
        self.logger.info("运行share2otu,将shared文件转化为otu")
        dir_ = os.listdir(mothur_dir)
        sub_sampled_shared = ""
        for file_ in dir_:
            if re.search(r'subsample', file_):
                sub_sampled_shared = os.path.join(mothur_dir, file_)
                break
        if sub_sampled_shared == "":
            self.set_error("mothur未能成功的产生subsample结果", code="32705503")
        match = re.search(r"(^.+)(\..+$)", self.basename)
        prefix = match.group(1)
        suffix = match.group(2)
        sub_sampled_otu = os.path.join(self.work_dir, "output", prefix + ".subsample" + suffix)
        sub_sampled_otu_renamed = os.path.join(self.work_dir, prefix + ".subsample" + suffix + ".rename")
        cmd = self.shared2otu_path + " -l 0.97 -i " + sub_sampled_shared + " -o " + sub_sampled_otu_renamed
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError:
            self.set_error("shared2otu.pl 运行出错", code="32705504")

        # 将新的OTU表用就有的OTU表替换掉
        with open(sub_sampled_otu_renamed, "rb") as r, open(sub_sampled_otu, "wb") as w:
            line = r.next()
            w.write(line)
            for line in r:
                line = line.split("\t")
                line[0] = new2old[line[0]]
                w.write("\t".join(line))
        self.option("out_otu_table").set_path(sub_sampled_otu)
        self.option("out_otu_table").check()
        if self.has_tax is True:
            sub_sampled_otu_obj = OtuTableFile()
            sub_sampled_otu_obj.set_path(sub_sampled_otu)
            sub_sampled_otu_obj.get_info()
            tmp_name = sub_sampled_otu + ".tmp"
            sub_sampled_otu_obj.add_taxonomy(tax_dic, sub_sampled_otu, tmp_name)
            os.remove(sub_sampled_otu)
            shutil.copy2(tmp_name, sub_sampled_otu)
            os.remove(tmp_name)

    def run(self):
        """
        运行
        """
        super(SubSampleTool, self).run()
        self.sub_sample()
        self.end()
