#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess
import shutil
import json
import glob
import os
import pandas as pd
import collections
from collections import OrderedDict
import re
import datetime
import unittest
import random


class VcfFilterSamtoolsAgent(Agent):
    """
    Annovar:用对处理vcf格式文件/注释突变信息
    version 1.0
    author: qindanhua
    last_modify: 2016.12.30
    """

    def __init__(self, parent):
        super(VcfFilterSamtoolsAgent, self).__init__(parent)
        options = [
            #{"name": "ref_genome", "type": "string"},  # 参考基因组类型
            {"name": "input_file", "type": "infile", "format": "gene_structure.vcf,gene_structure.vcf_dir"},  # 输入文件
            #{"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 输入文件,参考基因组文件
            #{"name": "combine_vcf", "type": "bool", "default": False},  # 输入文件

        ]
        self.add_option(options)
        self.step.add_steps('vcffiltersamtools')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.vcffiltersamtools.start()
        self.step.update()

    def step_end(self):
        self.step.vcffiltersamtools.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("input_file").is_set:
            raise OptionError("请输入VCF格式文件", code="33709602")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 11
        self._memory = '10G'

    def end(self):
        super(VcfFilterSamtoolsAgent, self).end()


class VcfFilterSamtoolsTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(VcfFilterSamtoolsTool, self).__init__(config)


    def vcf_filter(self):
        with open(self.vcf_path,"r") as rv,open(os.path.join(self.output_dir,"final.vcf"),"w") as fv:
            for line in rv.readlines():
                if line.startswith("#"):
                    fv.write(line)
                else:
                    line_info = line.strip().split()
                    qual = float(line_info[5])
                    va_dp=0.0
                    for info in line_info[9:]:
                        if info.strip().split(":")[0] != "./." and info.strip().split(":")[0] != "0/0":
                            va_dp += float(info.strip().split(":")[2])
                    qd=qual/va_dp
                    mapping_quality = re.search("MQ=[\d]*", line_info[7])
                    mapping_quality = mapping_quality.group().split("=")[-1]
                    if qd >= 2.0 and float(mapping_quality) >= 40.0 :
                        fv.write(line)

    def run(self):
        super(VcfFilterSamtoolsTool, self).run()
        self.vcf_path = self.option("input_file").prop["path"]
        self.vcf_filter()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "vcf_filter_samtools" + str(random.randint(1, 10000))+"-xxx",
            "type": "tool",
            "name": "ref_rna_v2.vcf_filter_samtools",
            "instant": False,
            "options": dict(
                #ref_dict=test_dir + "/" + "Mus_musculus.GRCm38.dna_rm.toplevel.clean.dict",
                #bam_list="/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/bamlist_new",
                #call_type="sentieon",
                #ref_fasta=test_dir+"/"+"Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                # scm="complete",
                # scd="correlation",
                # corr_method='pearson',
                #output=None,
                #ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                #input_file="/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna/output_vcf",
                #ref_genome="customer_mode",
                #combine_vcf=True,
                #des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/biomart/Mus_musculus.GRCm38.biomart_gene.txt",
                #des_type="type1"
                #des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Capsicum_annuum/NCBI/biomart/GCF_000710875.1_Pepper_Zunla_1_Ref_v1.0.biomart",
                #des_type="type3",
                #ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Capsicum_annuum/NCBI/gtf/GCF_000710875.1_Pepper_Zunla_1_Ref_v1.0.gtf",
                #ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Capsicum_annuum/NCBI/dna/GCF_000710875.1_Pepper_Zunla_1_Ref_v1.0_genomic.fna",
                #input_file="/mnt/ilustre/users/sanger-dev/workspace/20190617/Snp_tsg_33912_1820_7969/SnpRna/output_vcf"
                input_file="/mnt/ilustre/users/sanger-dev/workspace/20190620/Snp_tsg_33912_3063_9782/SamRna/BcftoolVcf/output/pop.variant.vcf"
                #input_file="/mnt/ilustre/users/sanger-dev/workspace/20190621/Snp_tsg_34423_8856_3754/CallSnpIndel/GvcfTypingV2/output/pop.variant.vcf"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()