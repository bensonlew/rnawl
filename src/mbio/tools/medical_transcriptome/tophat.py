# -*- coding: utf-8 -*-
# __author__ = "shijin,shicaiping"

import os
import unittest

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class TophatAgent(Agent):
    """
    tophat  
    version 2.0
    last_modify: by shicaiping at 20180508
    """

    def __init__(self, parent):
        super(TophatAgent, self).__init__(parent)
        options = [
            {"name": "ref_genome", "type": "string"},
            {"name": "genome_version", "type": "string", "default": "Custom"},  # 参考基因组版本
            {"name": "genome_annot_version", "type": "string", "default": "Custom"},  # 参考基因组注释版本
            {"name": "mapping_method", "type": "string"},  # 比对软件，tophat or hisat
            {"name": "seq_method", "type": "string"},
            {"name": "single_end_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "left_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "right_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "bam_output", "type": "outfile", "format": "align.bwa.bam"},
            {"name": "assemble_method", "type": "string", "default": "none"},
            {"name": "sample", "type": "string"},
            {"name": "num_threads", "type": "int", "default": 8},
            {"name": "mate_std", "type": "int", "default": 50},  # 末端配对插入片段长度标准差
            {"name": "mid_dis", "type": "int", "default": 50},  # 两个成对引物间的距离中间值
            {"name": "result_reserved", "type": "int", "default": 5},  # 最多保留的比对结果数目
            {"name": "strand_specific", "type": "bool", "default": False}
        ]
        self.add_option(options)
        self.step.add_steps("tophat")
        self.on("start", self.step_start)
        self.on("end", self.step_end)
        self._memory_increase_step = 40

    def step_start(self):
        self.step.tophat.start()
        self.step.update()

    def step_end(self):
        self.step.tophat.finish()
        self.step.update()

    def check_options(self):
        if self.option("seq_method") == "PE":
            if self.option("single_end_reads").is_set:
                raise OptionError("上传的是单端测序的序列，请上传双端序列")
            elif not (self.option("left_reads").is_set and self.option("right_reads").is_set):
                raise OptionError("缺少某端序列")
        else:
            if not self.option("single_end_reads").is_set:
                raise OptionError("请上传单端序列")
            elif self.option("left_reads").is_set or self.option("right_reads").is_set:
                raise OptionError("只需要有单端的序列")
        if not self.option("assemble_method").lower() in ["cufflinks", "stringtie", "none"]:
            raise OptionError("请选择拼接软件")

    def set_resource(self):
        self._cpu = 1
        if self.option("seq_method") == "PE":
            self._memory = "{}G".format((os.path.getsize(self.option("left_reads").prop["path"]) + os.path.getsize(
                self.option("right_reads").prop["path"])) / 1024 ** 3 * 4 + 40)
        else:
            self._memory = "{}G".format(
                (os.path.getsize(self.option("single_end_reads").prop["path"])) / 1024 ** 3 * 4 + 40)

    def end(self):
        super(TophatAgent, self).end()


class TophatTool(Tool):
    def __init__(self, config):
        super(TophatTool, self).__init__(config)
        self.program = {
            "tophat2": "bioinfo/align/tophat-2.1.1/tophat-2.1.1.Linux_x86_64/tophat2",
            "samtools": "bioinfo/align/samtools-1.8/samtools"
        }
        self.path = {
            "bam": os.path.join(self.output_dir, "{}.bam".format(self.option("sample")))
        }

    def run(self):
        super(TophatTool, self).run()
        self.pre_tophat2()
        self.end()

    def pre_tophat2(self):
        db = Config().get_mongo_client(mtype="ref_rna_v2", dydb_forbid=True)[Config().get_mongo_dbname("ref_rna_v2", dydb_forbid=True)]
        col = db["sg_genome_db"]
        self.logger.debug("name => {}".format(self.option("ref_genome")))
        self.logger.debug("assembly => {}".format(self.option("genome_version")))
        self.logger.debug("annot_version => {}".format(self.option("genome_annot_version")))
        genome_info = col.find_one({
            "name": self.option("ref_genome"),
            "assembly": self.option("genome_version"),
            "annot_version": self.option("genome_annot_version")
        })
        self.path["bowtie_index"] = os.path.join(
            self.config.SOFTWARE_DIR, "database/Genome_DB_finish/{}".format(genome_info["dna_index"]))
        self.run_tophat2()

    def run_tophat2(self):
        cmd = self.program["tophat2"]
        cmd += " --max-multihits {}".format(self.option("result_reserved"))
        if self.option("strand_specific"):
            cmd += " --library-type fr-firststrand"
        cmd += " --num-threads {}".format(self.option("num_threads"))
        cmd += " --mate-inner-dist {}".format(self.option("mid_dis"))
        cmd += " --mate-std-dev {}".format(self.option("mate_std"))
        cmd += " {}".format(self.path["bowtie_index"])
        if self.option("seq_method") == "PE":
            cmd += " {} {}".format(self.option("left_reads").path, self.option("right_reads").path)
        else:
            cmd += " {}".format(self.option("single_end_reads").path)
        self.logger.info("start running tophat2 for alignment")
        command = self.add_command("run_tophat2", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.samtools_sort()
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state("memory_limit", "memory is low!")
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running tophat2")
                self.samtools_sort()
        else:
            self.set_error("fail to run tophat2")

    def samtools_sort(self):
        cmd = "{} sort -m 2G -o {} -@ 4 {}".format(
            self.program["samtools"],
            os.path.join(self.work_dir, "accepted_hits.bam"),
            os.path.join(self.work_dir, "tophat_out/accepted_hits.bam")
        )
        command = self.add_command("samtools_sort", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.set_output()
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state("memory_limit", "memory is low!")
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running samtools")
                self.set_output()
        else:
            self.set_error("fail to run samtools")

    def set_output(self):
        if os.path.exists(self.path["bam"]):
            os.remove(self.path["bam"])
        os.link(os.path.join(self.work_dir, "accepted_hits.bam"), self.path["bam"])
        self.option("bam_output").set_path(self.path["bam"])


class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "tophat_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "ref_rna_v2.tophat",
            "instant": False,
            "options": dict(
                ref_genome="Arabidopsis_thaliana",
                genome_version="TAIR10",
                genome_annot_version="Ensembl_43",
                seq_method="PE",
                left_reads="/mnt/ilustre/users/sanger-dev/workspace/20200313/Refrna_tsg_37039/FastpRna/output/fastq"
                           "/TR2_1.clean.1.fastq",
                right_reads="/mnt/ilustre/users/sanger-dev/workspace/20200313/Refrna_tsg_37039/FastpRna/output/fastq"
                            "/TR2_1.clean.2.fastq",
                assemble_method="stringtie",
                sample="TR2_1"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
