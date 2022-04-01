# -*- coding: utf-8 -*-

# __author__ = "linfang.jin"
# time: 2017/1/25 9:33

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import re


class MapspliceMapAgent(Agent):
    def __init__(self, parent):
        super(MapspliceMapAgent, self).__init__(parent)
        options = [
            # 必须设置参数
            {"name": "seq_type", "type": "string", "default": "paired"},  # ‘paired’ 或者 ‘single’
            {"name": "bowtie_index_dir", "type": "string", "default": ""},
            {"name": "ref_dir", "type": "string", "default": ""},
            {"name": "_1_fq", "type": "infile", "format": "sequence.fastq"},
            {"name": "_2_fq", "type": "infile", "format": "sequence.fastq"},
            {"name": "threads", "type": "int", "default": 10},
            {"name": "out_dir", "type": "string", "default": self.output_dir},
            # 以下为比对相关参数
            # {"name": "ref_gtf", "type": "infile", "format": "ref_rna.assembly.gtf"},
            {"name": "segments_len", "type": "int", "default": 25},  # 限定在[18,25] 之间
            {"name": "min_map_len", "type": "int", "default": 50},
            {"name": "min_intron_len", "type": "int", "default": 50},
            {"name": "max_intron_len", "type": "int", "default": 300000},
            {"name": "non_cano_double_anchor", "type": "int", "default": 1},  # 0:不设定，1：设定
            {"name": "non_cano_single_anchor", "type": "int", "default": 1},  # 0:不设定，1：设定
            {"name": "splice_mismatch_max_num", "type": "int", "default": 1},  # 必须在[0,2]之间
            {"name": "max_append_mismatch", "type": "int", "default": 3},
            {"name": "max_insertion_len", "type": "int", "default": 3},  # 范围为[0,10]
            {"name": "max_deletion_len", "type": "int", "default": 6},
            {"name": "double_anchor_option", "type": "int", "default": 0},
            {"name": "single_anchor_option", "type": "int", "default": 0},
            {"name": "fusion_option", "type": "int", "default": 0},
            {"name": "fusion_non_canonical", "type": "int", "default": 0},
            {"name": "min_fusion_distance", "type": "int", "default": 10000}
        ]
        self.add_option(options)
        self.step.add_steps("mapsplice_map")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def check_options(self):
        """
        重写参数检查
        :return:
        """
        if not self.option("ref_dir"):
            raise OptionError("必须设置参考基因组文件夹 ")
        if not self.option("bowtie_index_dir"):
            raise OptionError("必须设置参考基因组索引文件夹")
        if not self.option("seq_type"):
            raise OptionError("必须设置测序类型： paired 或者 single")
        if not self.option("_1_fq"):
            raise OptionError("必须设置fastq文件")
        if self.option("seq_type") == "paired" and self.option("_2_fq") is None:
            raise OptionError("当测序类型为双端测序时，必须设置左右两端fastq文件")
        if self.option("segments_len") < 18 or self.option("segments_len") > 25:
            raise OptionError("reads 被分为的小片段（segments）长度应在18bp 和25bp 之间")
        if self.option("splice_mismatch_max_num") < 0 or self.option("splice_mismatch_max_num") > 2:
            raise OptionError("splice_mismatch_max_num 长度应在0bp 和2bp 之间")
        if self.option("min_intron_len") > self.option("max_intron_len"):
            raise OptionError("splice junctions 的最小长度不应该大于splice junction的最大长度")
        if self.option("max_insertion_len") > 10 or self.option("max_insertion_len") < 0:
            raise OptionError("max_insertion_len 应介于[0,10]之间")
        if self.option("double_anchor_option") not in (0, 1):
            raise OptionError("double_anchor_option 必须设置为0或1")
        if self.option("single_anchor_option") not in (0, 1):
            raise OptionError("single_anchor_option 必须设置为0或1")
        if self.option("fusion_option") not in (0, 1):
            raise OptionError(" fusion_option必须设置为0或1")
        if self.option("fusion_non_canonical") not in (0, 1):
            raise OptionError(" fusion_non_canonical 必须设置为0或1")

    def step_start(self):
        self.step.mapsplice_map.start()
        self.step.update()

    def step_end(self):
        self.step.mapsplice_map.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 10
        self._memory = "100G"

    def end(self):
        """
        agent结束后一些文件德操作

        :return:
        """
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"],
        # ])
        # result_dir.add_regexp_rules([
        #     ["ASEvents", "文件夹", "AS事件详细信息表文件夹"],
        #     ["MATS_output", "文件夹", "AS事件上reads信息"],
        #     ["commands.txt", "txt", "命令历史"],
        #     ["config.txt", "txt", "软件运行配置信息"],
        #     ["summary.txt", "txt", "rMATS结果总结文件"],
        # ])
        super(MapspliceMapAgent, self).end()


class MapspliceMapTool(Tool):
    def __init__(self, config):
        super(MapspliceMapTool, self).__init__(config)
        self._version = "v2.1.8"
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/MapSplice-v2.1.8/mapsplice.py"
        self.Python_path = "program/Python/bin/python"
        self.double_anchor_dic = {0: "", 1: "--non-canonical-double-anchor"}
        self.single_anchor_dic = {0: "", 1: "--non-canonical-single-anchor"}
        self.fusion_dic = {0: "", 1: "--fusion"}
        self.fusion_non_cano_dic = {0: "", 1: "--fusion-non-canonical"}

    def run_mapsplice_map(self):
        """
        运行mapsplice_map
        :return:
        """
        if self.option("seq_type") == "paired":
            cmd = "{} {} -p {} -o {}  --bam  -s {}  --min-map-len {} -i {} -I {} -m {} --max-append-mis {} --ins {}  --del {} --min-fusion-distance {}  {}  {}  {}  {}  -c {} -x {} -1 {} -2 {}".format(
                self.Python_path, self.script_path, self.option("threads"), self.output_dir,
                self.option("segments_len"), self.option("min_map_len"), self.option("min_intron_len"),
                self.option("max_intron_len"), self.option("splice_mismatch_max_num"),
                self.option("max_append_mismatch"),
                self.option("max_insertion_len"), self.option("max_deletion_len"),
                self.option("min_fusion_distance"),
                self.double_anchor_dic[self.option("double_anchor_option")],
                self.single_anchor_dic[self.option("single_anchor_option")],
                self.fusion_dic[self.option("fusion_option")],
                self.fusion_non_cano_dic[self.option("fusion_non_canonical")],
                self.option("ref_dir"), self.option("bowtie_index_dir"),
                self.option("_1_fq").path,
                self.option("_2_fq").path)
        elif self.option("seq_type") == "single":
            cmd = "{} {} -p {} -o {}  --bam  -s {}  --min-map-len {} -i {} -I {} -m {} --max-append-mis {} --ins {}  --del {} --min-fusion-distance {}  {}  {}  {}  {}  -c {} -x {} -1 {}".format(
                self.Python_path, self.script_path, self.option("threads"), self.output_dir,
                self.option("segments_len"), self.option("min_map_len"), self.option("min_intron_len"),
                self.option("max_intron_len"), self.option("splice_mismatch_max_num"),
                self.option("max_append_mismatch"),
                self.option("max_insertion_len"), self.option("max_deletion_len"),
                self.option("min_fusion_distance"),
                self.double_anchor_dic[self.option("double_anchor_option")],
                self.single_anchor_dic[self.option("single_anchor_option")],
                self.fusion_dic[self.option("fusion_option")],
                self.fusion_non_cano_dic[self.option("fusion_non_canonical")],
                self.option("ref_dir"), self.option("bowtie_index_dir"),
                self.option("_1_fq").path)

        self.logger.info("运行mapsplice_map")
        command = self.add_command("mapsplice_map_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("mapsplice_map运行完成")
        else:
            self.set_error("mapsplice_map运行出错!")

    def run(self):
        """
        运行
         :return:
        """
        super(MapspliceMapTool, self).run()
        self.run_mapsplice_map()
        self.end()
