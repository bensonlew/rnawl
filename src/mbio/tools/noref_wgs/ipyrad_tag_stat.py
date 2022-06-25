# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20190110

import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class IpyradTagStatAgent(Agent):
    """
    ipyrad 统计tag信息统计和Consensus信息统计
    """
    def __init__(self, parent=None):
        super(IpyradTagStatAgent, self).__init__(parent)
        options = [
            {"name": "cluster_dir", "type": "infile", "format": "bsa.dir", "required": True},  # data_clust_0.85文件夹
            {"name": "data_loci", "type": "infile", "format": "noref_wgs.list_file", "required": True},  # data.loci
            {"name": "total_sample_num", "type": "int", "required": True},  #总的样本数
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("cluster_dir").is_set:
            raise OptionError("请设置聚类文件夹data_clust_0.85", code="35500707")
        if not self.option("data_loci").is_set:
            raise OptionError("请设置ipyrad data_outfiles里的data_loci文件", code="35500708")
        if not self.option("total_sample_num"):
            raise OptionError("请设置总样本数", code="35500709")

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(IpyradTagStatAgent, self).end()


class IpyradTagStatTool(Tool):
    def __init__(self, config):
        super(IpyradTagStatTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.python_path = "miniconda2/bin/python"
        self.tag_stat = self.config.PACKAGE_DIR + "/noref_wgs/tagsinfo.pl"
        self.consensus_path = self.config.PACKAGE_DIR + "/noref_wgs/consensus_stat_ipyrad.py"

    def run_tag_stat(self):
        """
        统计tag信息表：tag_stat.xls和每个样本每个tag的深度表：sample_tag_depth.xls
        """
        tag_out = os.path.join(self.output_dir, "tag_stat.xls")
        self.tag_depth = os.path.join(self.work_dir, "sample_tag_depth.xls")
        cmd = "{} {} -in {} ".format(self.perl_path, self.tag_stat, self.option("cluster_dir").prop["path"])
        cmd += "-out {} -dep {}".format(tag_out, self.tag_depth)
        command = self.add_command("tag_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("tag信息表统计成功")
        else:
            self.set_error("tag信息表统计失败，请检查", code="35500705")

    def get_tag_depth_distribution(self):
        """
        统计每个样本的tag深度分布表
        """
        self.logger.info("开始进行每个样本的tag深度分布统计")
        depth_info = {}
        sample_list, depth_list = [], []
        with open(self.tag_depth, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                sample = item[0]
                depth = int(item[1])
                if sample not in sample_list:
                    depth_info[sample] = {}
                    sample_list.append(sample)
                if depth not in depth_info[sample].keys():
                    depth_info[sample][depth] = 0
                depth_info[sample][depth] += 1
                if depth not in depth_list:
                    depth_list.append(depth)
        depth_list.sort()
        with open(os.path.join(self.output_dir, "tag_depth.xls"), "w") as w:
            w.write("#Depth\t{}\n".format("\t".join(sample_list)))
            for depth in depth_list:
                s_depth = []
                for sample in sample_list:
                    try:
                        s_depth.append(str(depth_info[sample][depth]))
                    except:
                        s_depth.append(str(None))
                w.write(str(depth) + "\t" + "\t".join(s_depth) + "\n")

    def consensus_stat(self):
        """
        统计consensus结果和聚类分布结果
        """
        consen_stat = os.path.join(self.output_dir, "consensus_stat.xls")
        cover_path = os.path.join(self.output_dir, "consensus_coverage.xls")
        pop_tag = os.path.join(self.output_dir, "populations.tag")
        cmd = "{} {} -i {} ".format(self.python_path, self.consensus_path, self.option("data_loci").prop["path"])
        cmd += "-s {} -c {} -p {} -t {}".format(consen_stat, cover_path, pop_tag, str(self.option("total_sample_num")))
        command = self.add_command("consensus_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("consensus信息表统计成功")
        else:
            self.set_error("consensus信息表统计失败，请检查", code="35500706")

    def run(self):
        super(IpyradTagStatTool, self).run()
        self.run_tag_stat()
        self.get_tag_depth_distribution()
        self.consensus_stat()
        self.end()
