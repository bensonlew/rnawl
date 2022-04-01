# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/4/23'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import json
from mbio.packages.metagenomic.common import link_file
import shutil


class HgapAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(HgapAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.bam", "required": True},
            {"name": "genome_size", "type": "string", "default": "5000000", "required": True},
            {"name": "sample_name", "type": "string", "required": True},
            {"name": "scaf", "type": "outfile", "format": "sequence.fasta", "required": True},
            {"name": "cpu", "type": "int", "default": 8},
            {"name": "mem", "type": "int", "default": 50}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        # modified check_option
        pass

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = self.option("cpu")
        self._memory = "%sG" % self.option("mem")


class HgapTool(Tool):
    def __init__(self, config):
        super(HgapTool, self).__init__(config)
        self.smrt_tool = "bioinfo/Genomic/Sofware/smrttools/smrtcmds/bin/"
        if self.option("genome_size") == self.option("genome_size").upper():
            self.genome_size = self.option("genome_size")  # eg. "5000000"
        else:
            self.genome_size = self.option("genome_size").upper().rstrip("MBP") + "000000"  # eg. "5Mb" "5m" "5M" "5MB" "5000000bp"

    def run_dataset(self):
        if os.path.isfile("data.xml"):
            return 0
        cmd = "{}dataset create --type SubreadSet --force data.xml {}".format(self.smrt_tool, self.option("bam_file").prop["path"])
        command = self.add_command("dataset", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_dataset运行完成")
        else:
            self.set_error("run_dataset运行出错！")

    def run_xml(self):
        if os.path.isfile("preset_hgap_new.json"):
            return 0
        cmd = "{}pbsmrtpipe show-template-details pbsmrtpipe.pipelines.polished_falcon_fat -j preset_hgap.json".format(self.smrt_tool)
        command = self.add_command("xml", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_xml完成")
            self.change_opt("preset_hgap.json", "preset_hgap_new.json", {"falcon_ns.task_options.HGAP_GenomeLength_str": self.genome_size})
        else:
            self.set_error("run_dataset运行出错！")

    def change_opt(self, input, output, update_dic):
        with open(input, "r") as load_f, open(output, "w") as dump_f:
            load_dict = json.load(load_f)
            load_dict["taskOptions"].update(update_dic)
            json.dump(load_dict, dump_f, indent=4)

    def run_hgap(self):
        if os.path.isdir(os.path.join(self.work_dir, "result")):
            shutil.rmtree(os.path.join(self.work_dir, "result"))
        cmd = "{}pbsmrtpipe pipeline-id -e eid_subread:data.xml -o result --preset-json preset_hgap_new.json --local-only pbsmrtpipe.pipelines.polished_falcon_fat".format(self.smrt_tool)
        command = self.add_command("hgap", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("hgap完成")
        else:
            self.end()
            # self.set_error("run_hgap运行出错！")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        from_f = os.path.join(self.work_dir, "result/tasks/genomic_consensus.tasks.variantcaller-0/consensus.fasta")
        to_f = os.path.join(self.output_dir, self.option("sample_name") + ".scaf.fna")
        link_file(from_f, to_f)
        self.option("scaf").set_path(to_f)

    def run(self):
        super(HgapTool, self).run()
        self.run_dataset()
        self.run_xml()
        self.run_hgap()
        self.set_output()
        self.end()