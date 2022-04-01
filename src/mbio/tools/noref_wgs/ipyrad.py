# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20181217

import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class IpyradAgent(Agent):
    """
    ipyrad
    """
    def __init__(self, parent=None):
        super(IpyradAgent, self).__init__(parent)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "bsa.dir", "required": True},
            # fastq路径list.txt文件，第一列分析样本名，第二列当前样本名，第三列fastq_l路径，第四列fastq_r路径，第五列文库类型
            {"name": "analysis_step", "type": "string", "default": "123"},  # 分析步骤
            {"name": "enzyme_method", "type": "string", "required": True},  # 酶切方案,GBS/RAD
            {"name": "cutseq", "type": "string", "required": True},  # 酶切序列, 如：TCGA, TTAA
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("enzyme_method") not in ["RAD", "GBS"]:
            raise OptionError("酶切方案: %s只能是RAD/GBS" , variables=( self.option("enzyme_method")), code="35500905")
        if self.option("enzyme_method") == "GBS":
            if len(self.option("cutseq").split(", ")) != 2:
                raise OptionError("酶切方案为GBS的时候cutseq:%s 必须为2个" , variables=( self.option("enzyme_method")), code="35500906")

    def set_resource(self):
        self._cpu = 16
        self._memory = "250G"

    def end(self):
        super(IpyradAgent, self).end()


class IpyradTool(Tool):
    def __init__(self, config):
        super(IpyradTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/noRefWGS/miniconda2/bin")
        self.ipyrad = "bioinfo/noRefWGS/miniconda2/bin/ipyrad"

    def get_ipyrad_params(self):
        """
        得到ipyrad的params文件
        """
        if self.option("enzyme_method") == "GBS":
            self.method = "pairddrad"
            cutseq = self.option("cutseq") + ","
        else:
            self.method = "rad"
            cutseq = self.option("cutseq")
        self.params_txt = os.path.join(self.work_dir, "params-data.txt")
        with open(self.params_txt, "w") as w:
            w.write("------- ipyrad params file (v.0.7.28)-------------------------------------------\n")
            w.write("data                                   ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps\n")
            w.write(self.output_dir + "                     ## [1] [project_dir]: Project dir (made in curdir if not present)\n")
            w.write("                                       ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files\n")
            w.write("                                       ## [3] [barcodes_path]: Location of barcodes file\n")
            w.write(self.option("fastq_dir").prop["path"] + "/*     ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files\n")
            w.write("denovo                                 ## [5] [assembly_method]: Assembly method (denovo, reference, denovo+reference, denovo-reference)\n")
            w.write("                                       ## [6] [reference_sequence]: Location of reference sequence file\n")
            w.write(self.method + "                          ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.\n")
            w.write(cutseq + "               ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)\n")
            w.write("5                                      ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read\n")
            w.write("33                                     ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)\n")
            w.write("6                                      ## [11] [mindepth_statistical]: Min depth for statistical base calling\n")
            w.write("6                                      ## [12] [mindepth_majrule]: Min depth for majority-rule base calling\n")
            w.write("10000                                  ## [13] [maxdepth]: Max cluster depth within samples\n")
            w.write("0.85                                   ## [14] [clust_threshold]: Clustering threshold for de novo assembly\n")
            w.write("0                                      ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes\n")
            w.write("0                                      ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)\n")
            w.write("35                                     ## [17] [filter_min_trim_len]: Min length of reads after adapter trim\n")
            w.write("2                                      ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences\n")
            w.write("5, 5                                   ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)\n")
            w.write("8, 8                                   ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)\n")
            w.write("4                                      ## [21] [min_samples_locus]: Min # samples per locus for output\n")
            w.write("20, 20                                 ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)\n")
            w.write("8, 8                                   ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)\n")
            w.write("0.5                                    ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)\n")
            w.write("0, 0, 0, 0                             ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)\n")
            w.write("0, 0, 0, 0                             ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)\n")
            w.write("p, s, v                                ## [27] [output_formats]: Output formats (see docs)\n")
            w.write("                                       ## [28] [pop_assign_file]: Path to population assignment file\n")

    def run_ipyrad(self):
        """
        运行ipyrad
        """
        cmd = "{} -p {} -s {} -r --MPI -f -t 2 -c 8".format(self.ipyrad, self.params_txt, self.option("analysis_step"))
        command = self.add_command("ipyrad-" + self.option("analysis_step"), cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ipyrad-" + self.option("analysis_step") + "成功")
        else:
            self.set_error("ipyrad-" + self.option("analysis_step") + "失败，请检查", code="35500903")

    def run(self):
        super(IpyradTool, self).run()
        self.get_ipyrad_params()
        self.run_ipyrad()
        self.end()
