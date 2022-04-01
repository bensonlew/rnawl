# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
import pandas as pd
from mbio.packages.metagenomic.common import link_file


class AssembleAssess2Module(Module):
    """
    微生物基因组组装结果评估
    author: guhaidong
    last_modify: 2019.04.26
    """
    def __init__(self, work_id):
        super(AssembleAssess2Module, self).__init__(work_id)
        options = [
            {"name": "seq_scaf", "type": "infile", "format": "sequence.fasta"},  # 最佳Gapcloser的scaffold文件
            {'name': 'sample_name', "type": "string"},  # 样本名
            {'name': 'base_size', "type": "float"}, # 理论基因组大小
            {'name': "abund_table", "type": "infile", "format": "sequence.profile_table"},
            {'name': "fq1", "type": "infile", "format": "sequence.fastq"},
            {"name": "fq2", "type": "infile", "format": "sequence.fastq"},
            {"name": "scaffold", "type": "outfile", "format": "sequence.fasta"},  # 输出文件
        ]
        self.step.add_steps('scaf_agp_contig', 'bac_stat')
        self.scaf_agp_contig = self.add_tool('assemble.scaf_agp_contig')
        self.bac_stat = self.add_tool('bacgenome.bac_assemble_stat')
        self.coverage = self.add_module("bacgenome.mapping")
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('seq_scaf').is_set:
            raise OptionError('必须输入seq_scaf序列文件', code="21400301")
        if not self.option('sample_name'):
            raise OptionError('必须输入sample_name样品名称！', code="21400302")

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run_gaploser_scaf(self):
        opts = {
            "seq_scaf": self.option('seq_scaf'),
            "sample_name": self.option('sample_name'),
        }
        if not self.option("fq1").is_set:  # 做三代数据结果统计，scaffold事先已排好顺序
            opts["sort_scaf"] = False
        self.scaf_agp_contig.set_options(opts)
        self.scaf_agp_contig.on('end', self.set_output, 'scaf_agp_contig')
        self.scaf_agp_contig.run()
        self.step.scaf_agp_contig.finish()
        self.step.update()

    def run_bac_stat(self):
        opts = {
            "scaf_seq": self.scaf_agp_contig.option('scaffold'),
            "cont_seq": self.scaf_agp_contig.option('contig'),
            "sample_name": self.option('sample_name'),
        }
        self.bac_stat.set_options(opts)
        self.bac_stat.on('end', self.set_output, 'bac_stat')
        self.bac_stat.run()
        self.step.bac_stat.finish()
        self.step.update()

    def run_coverage(self):
        if self.option("fq1").is_set:
            self.coverage.set_options({
                "ref_fa": self.scaf_agp_contig.option('scaffold'),
                "fq1": self.option("fq1"),
                "fq2": self.option("fq2")
            })
            self.coverage.on("end", self.set_output, "coverage")
            self.coverage.run()

    def run(self):
        """
        运行
        :return:
        """
        super(AssembleAssess2Module, self).run()
        self.scaf_agp_contig.on('end',self.run_bac_stat)
        self.bac_stat.on("end", self.run_coverage)
        self.run_gaploser_scaf()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        allfiles = os.listdir(dirpath)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if event["data"] == "scaf_agp_contig":
            self.linkdir(self.scaf_agp_contig.output_dir , self.output_dir + '/assembly')
            self.option('scaffold',self.output_dir + '/assembly/' + self.option('sample_name') + '_scaf.fna')
        if event["data"] == "bac_stat":
            self.linkdir(self.bac_stat.output_dir + '/summary', self.output_dir + '/assembly')
            self.linkdir(self.bac_stat.output_dir + '/len', self.output_dir + '/len')
            # self.linkdir(self.bac_stat.work_dir + '/scaffold',self.output_dir + '/scaffold')
            # self.linkdir(self.bac_stat.work_dir + '/contig', self.output_dir + '/contig')
            if not self.option("fq1").is_set:
                self.set_output({"data": "chr_coverage"})
        if event["data"] == "coverage":
            data = pd.read_table(self.coverage.option("cov_out").prop["path"])
            cov = data[["contigName", "totalAvgDepth"]]
            detail_data = pd.read_table(self.bac_stat.output_dir + '/summary/' + self.option("sample_name") + "_assembly_scaffold_details.xls", header=None, index_col=0)
            detail_data = detail_data[range(1, 4)]
            new_detail = pd.merge(detail_data, cov, left_index=True, right_on="contigName", how="left")
            new_detail.set_index("contigName", drop=True, inplace=True)
            # new_detail.to_csv(self.output_dir + "/" + self.option("sample_name") + "_assembly_scaffold_details.xls", sep="\t", header=False, index=True)
            new_detail.to_csv(self.output_dir + "/assembly/" + self.option("sample_name") + "_assembly_scaffold_details.xls", sep="\t", header=False, index=True)
            stat_data = pd.read_table(self.bac_stat.output_dir + "/summary/" + self.option("sample_name") + "_assembly_summary.xls")
            sample_cov = (data["totalAvgDepth"] * data["contigLen"]).sum() / data["contigLen"].sum()
            if "Coverage" in stat_data.columns:
                stat_data.loc[:, "Coverage"] = sample_cov
            else:
                stat_data.insert(2, "Coverage", [sample_cov])
            stat_data.to_csv(self.output_dir + "/assembly/" + self.option("sample_name") + "_assembly_summary.xls", sep="\t", header=True, index=False)
            link_file(self.coverage.coverage.option("maxbin_depth").prop["path"], self.output_dir + '/assembly/' + self.option("sample_name") + ".abund")
            link_file(self.coverage.option("cov_out2").prop["path"], self.output_dir + "/assembly/" + self.option("sample_name") + ".depth")
            self.end()
        if event["data"] == "chr_coverage":
            stat_data = pd.read_table(self.bac_stat.output_dir + "/summary/" + self.option("sample_name") + "_assembly_summary.xls", header=0)
            sample_cov = self.option("base_size") / float(stat_data.iloc[0,1])
            if "Coverage" in stat_data.columns:
                stat_data.loc[:, "Coverage"] = sample_cov
            else:
                stat_data.insert(2, "Coverage", [sample_cov])
            self.logger.info("=====Coverage")
            self.logger.info(stat_data)
            stat_data.to_csv(self.output_dir + "/assembly/" + self.option("sample_name") + "_assembly_summary.xls", sep="\t", header=True, index=False)
            detail_data = pd.read_table(self.bac_stat.output_dir + '/summary/' + self.option("sample_name") + "_assembly_scaffold_details.xls", header=None, index_col=0)
            detail_data = detail_data[range(1, 4)]
            self.logger.info(">>>>>>>>>>")
            if self.option("abund_table").is_set:
                data = pd.read_table(self.option("abund_table").prop["path"])
                new_detail = detail_data.reset_index()
                new_detail = pd.concat([new_detail, data["totalAvgDepth"]], axis=1)
                new_detail = new_detail.set_index(0)
                self.logger.info(new_detail)
                new_detail.to_csv(self.output_dir + "/assembly/" + self.option("sample_name") + "_assembly_scaffold_details.xls", sep="\t", header=False, index=True)
                new_detail.index.name = "contigName"
                new_detail["totalAvgDepth"].to_csv(self.output_dir + '/assembly/' + self.option("sample_name") + ".abund", sep="\t", header=True, index=True)
                mean_depth_path = self.option("abund_table").prop["path"].rstrip("abund") + "depth"
                link_file(mean_depth_path, self.output_dir + "/assembly/" + self.option("sample_name") + ".depth")
            else:
                detail_data["totalAvgDepth"] = sample_cov
                self.logger.info(detail_data)
                detail_data.to_csv(self.output_dir + "/assembly/" + self.option("sample_name") + "_assembly_scaffold_details.xls", sep="\t", header=False, index=True)
                detail_data.index.name = "contigName"
                detail_data["totalAvgDepth"] = sample_cov
                detail_data["totalAvgDepth"].to_csv(self.output_dir + "/assembly/" + self.option("sample_name") + ".abund", sep="\t", header=True, index=True)
            self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AssembleAssess2Module, self).end()