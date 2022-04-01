# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.03

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class BwaMemAgent(Agent):
    """
    软件: bwa
    bwa的meme方法, 注：ref.fa要好建索引
    """
    def __init__(self, parent):
        super(BwaMemAgent, self).__init__(parent)
        options = [
            {"name": "fastq_l", "type": "infile", "format": "sequence.fastq"},  # 左端fastq序列
            {"name": "fastq_r", "type": "infile", "format": "sequence.fastq"},  # 右端fastq序列
            {"name": "sample_name", "type": "string"},  #样本名称
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "num", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fastq_l").is_set:
            raise OptionError("请设置左端fastq序列", code="34501101")
        if not self.option("fastq_r").is_set:
            raise OptionError("请设置右端fastq序列", code="34501102")
        if not self.option("ref_fa").is_set:
            raise OptionError("请设置参考序列", code="34501103")
        else:
            ref_dir = os.path.dirname(self.option("ref_fa").prop["path"])
            base_name = os.path.basename(self.option("ref_fa").prop["path"])
            amb_path = os.path.join(ref_dir, base_name + ".amb")
            ann_path = os.path.join(ref_dir, base_name + ".ann")
            bwt_path = os.path.join(ref_dir, base_name + ".bwt")
            fai_path = os.path.join(ref_dir, base_name + ".fai")
            pac_path = os.path.join(ref_dir, base_name + ".pac")
            sa_path = os.path.join(ref_dir, base_name + ".sa")
            if not os.path.exists(amb_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.amb", code="34501104")
            if not os.path.exists(ann_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.ann", code="34501105")
            if not os.path.exists(bwt_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.bwt", code="34501106")
            if not os.path.exists(fai_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.fai", code="34501107")
            if not os.path.exists(pac_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.pac", code="34501108")
            if not os.path.exists(sa_path):
                raise OptionError("ref_fa的索引不全，缺少ref.fa.sa", code="34501109")

    def set_resource(self):
        self._cpu = 12
        self._memory = "30G"

    def end(self):
        super(BwaMemAgent, self).end()


class BwaMemTool(Tool):
    def __init__(self, config):
        super(BwaMemTool, self).__init__(config)
        self.bwa_path = self.config.SOFTWARE_DIR + "/bioinfo/wes/bwa/bwa"
        self.bwa_sh_path = "bioinfo/WGS/bwa_mapping.sh"

    def run_bwa_mem(self):
        """
        bwa mem,将reads比对到参考序列ref_fa上
        """
        if self.option("sample_name"):
            sample_name = self.option("sample_name")
        else:
            sample_name = os.path.basename(self.option("fastq_l").prop["path"]).split(".fastq")[0]
        name = sample_name.split("-")[0]
        header = "@RG\\tID:{}\\tLG:{}\\tLB:1\\tPL:illumina\\tSM:{}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq".format(self.option("num"), name, name)
        sam_path = self.output_dir + "/" + sample_name + ".sam"
        if os.path.exists(self.output_dir + "/" + sample_name + ".sam"):
            os.remove(self.output_dir + "/" + sample_name + ".sam")
        # header = "@RG\\tID:{}\\tLB:LB1\\tSM:{}\\tPL:ILLUMINA".format(self.option("sample_name"), self.option("sample_name"))
        # cmd = "{} {} \"@RG\\tID:{}\\tLB:LB1\\tSM:{}\\tPL:ILLUMINA\" {} {} {} {} {}".format(self.bwa_sh_path, self.bwa_path,\
        #        self.option("sample_name"), self.option("sample_name"), self.hg19_ref, self.option("fastq_l").prop["path"],\
        #        self.option("fastq_r").prop["path"], self.work_dir + "/aln-pe.sam",
        #        8)  # -t 8；使用8个threads
        cmd = "{} {} \"@RG\\tID:{}\\tLG:{}\\tLB:1\\tPL:illumina\\tSM:{}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq\" {} {} {} {} {}"\
               .format(self.bwa_sh_path, self.bwa_path, self.option("num"), name, name, self.option("ref_fa").prop["path"],\
               self.option("fastq_l").prop["path"], self.option("fastq_r").prop["path"], sam_path, 8)
        self.logger.info(cmd)
        # cmd = "{} {}".format(self.bwa_sh_path, self.bwa_path)
        # cmd += " '@RG\\tID:1\\tLG:{}\\tLB:1\\tPL:illumina\\tSM:{}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq'".format(sample_name, sample_name)
        # cmd + " {} ".format(self.option("ref_fa").prop["path"])
        # cmd += " {} {} {} {}".format(self.option("fastq_l").prop["path"], self.option("fastq_r").prop["path"], sam_path, 8)  # -t 8；使用8个threads
        command = self.add_command("bwa_mapping_mem", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bwa mapping完成")
        else:
            self.set_error("bwa mapping失败", code="34501101")

    def run(self):
        super(BwaMemTool, self).run()
        self.run_bwa_mem()
        self.end()
