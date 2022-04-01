# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/4/28'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from Bio import SeqIO


class PilonAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(PilonAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "frags_bam", "type": "infile", "format": "align.bwa.bam", "required": True},
            {"name": "scf_seq", "type": "outfile", "format": "sequence.fasta"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        pass

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = "40G"


class PilonTool(Tool):
    def __init__(self, config):
        super(PilonTool, self).__init__(config)
        self.pilon = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/Pilon/pilon-1.22.jar "
        self.java = "program/sun_jdk1.8.0/bin/java "

    def run_pilon(self):
        """
        description
        :return:
        """
        cmd = self.java + " -d64 -Xmx15G -jar "  + self.pilon + " --genome %s --frags %s --fix snps,indels --output pilon_polished --vcf" % (self.option("ref_fa").prop["path"], self.option("frags_bam").prop["path"])
        command = self.add_command("pilon", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("pilon success")
        else:
            self.set_error("pilon error")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        # 将序列id按原数据名称进行修改
        seqrecords = SeqIO.parse(os.path.join(self.work_dir, "pilon_polished.fasta"), "fasta")
        ref_fas = SeqIO.to_dict(SeqIO.parse(self.option("ref_fa").prop["path"], "fasta"))
        newrecords = []
        for seqrecord in seqrecords:
            seq_id = seqrecord.id.replace("_pilon", "")
            seq_id = seq_id
            seqrecord.id = seq_id
            seqrecord.description = ref_fas[seq_id].description
            newrecords.append(seqrecord)
        new_fa = os.path.join(self.output_dir, "pilon_polished.fasta")
        SeqIO.write(newrecords, new_fa, "fasta")
        self.option("scf_seq").set_path(new_fa)

    def run(self):
        super(PilonTool, self).run()
        self.run_pilon()
        self.set_output()
        self.end()