# -*- coding: utf-8 -*-
# __author__ = "wangzhaoyue"

"""fastp质控工具 """
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class FastpAgent(Agent):
    """
    fastp
    """
    def __init__(self, parent=None):
        super(FastpAgent, self).__init__(parent)
        options = [
            {"name": "fq1", "type": "infile", "format": "sequence.fastq"},  # 原始数据
            {"name": "fq2", "type": "infile", "format": "sequence.fastq"},  # 原始数据
            {"name": "qualified_quality_phred", "type": "string", "default": "15"},  # -q,一个碱基合格的质量值,默认15表示phred质量> = Q15是合格的。
            {"name": "length_required", "type": "string", "default": "30"},  # -l,长度过滤参数，比此值短的读取将被丢弃
            {"name": "cut_by_quality5", "type": "string"},  # -5,根据前面(5 ")的质量，允许每个读切割，默认是禁用的
            {"name": "cut_by_quality3", "type": "string"},  # -3,根据后面(3 ")的质量，允许每个读切割，默认是禁用的
            {"name": "cut_mean_quality", "type": "string", "default": "20"},  # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {"name": "n_base_limit", "type": "string", "default": "5"},  # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {"name": "compression", "type": "string", "default": "2"},  # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {"name": "thread", "type": "string", "default": "3"},  # -w,线程数
            {"name": "adapter_sequence", "type": "string"},  # --adapter_sequence,the adapter for read1, default:AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
            {"name": "adapter_sequence_r2", "type": "string"},  # --adapter_sequence_r2,the adapter for read2 (PE data only),default:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
            {"name": "sample_name", "type": "string"},  # 样本名
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检测
        """
        if not self.option("fq1"):
            raise OptionError("必须输入1端序列", code="34502501")
        if not self.option("fq2"):
            raise OptionError("必须输入2端序列", code="34502502")
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 12
        self._memory = "15G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(FastpAgent, self).end()


class FastpTool(Tool):
    def __init__(self, config):
        super(FastpTool, self).__init__(config)
        self.fastp_path = "/bioinfo/seq/fastp "

    def run(self):
        super(FastpTool, self).run()
        self.run_fastp()
        self.set_output()
        self.end()

    def run_fastp(self):
        """
        运行fastp,进行质控
        """
        if self.option("sample_name"):
            sample_name = self.option("sample_name")
        else:
            sample_name = os.path.basename(self.option("fq1").prop["path"]).split(".")[0]
        cmd = self.fastp_path + " -i %s -I %s -o %s -O %s -q %s -l %s -M %s -n %s -j %s -h %s -z %s -w %s" % (
        self.option("fq1").prop["path"], self.option("fq2").prop["path"],
        self.work_dir + "/" + sample_name + ".clean.1.fastq.gz", self.work_dir + "/" + sample_name + ".clean.2.fastq.gz",
        self.option("qualified_quality_phred"),  self.option("length_required"),
        self.option("cut_mean_quality"),self.option("n_base_limit"),
        self.work_dir + "/" + sample_name + ".json", self.work_dir + "/" + sample_name + ".html", self.option("compression"),
        self.option("thread"))
        if self.option("cut_by_quality5"):
            cmd += " -5 %s" % (self.option("cut_by_quality5"))
        if self.option("cut_by_quality3"):
            cmd += " -3 %s" % (self.option("cut_by_quality3"))
        if self.option("adapter_sequence"):
            cmd += " --adapter_sequence %s" % (self.option("adapter_sequence"))
        if self.option("adapter_sequence_r2"):
            cmd += " --adapter_sequence_r2 %s" % (self.option("adapter_sequence_r2"))
        self.logger.info("运行fastp，进行质控")
        command = self.add_command("fastp_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fastp运行完成")
        else:
            self.set_error("fastp运行出错!", code="34502501")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            all_files = os.listdir(self.work_dir)
            for files in all_files:
                if files.endswith(".fastq.gz"):
                    if os.path.exists(self.output_dir + "/" + files):
                        os.remove(self.output_dir + "/" + files)
                    os.link(self.work_dir + "/" + files, self.output_dir + "/" + files)
            self.logger.info("设置fastp分析结果目录成功")
        except Exception as e:
            self.logger.info("设置fastp分析结果目录失败{}".format(e))
            self.set_error("设置fastp分析结果目录失败%s", variables=(e), code="34502502")
