# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from Bio import SeqIO
from biocluster.config import Config
from biocluster.module import Module


class WholeSnpModule(Module):
    def __init__(self, work_id):
        super(WholeSnpModule, self).__init__(work_id)
        options = [
            {"name": "ref_genome_custom", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "ref_genome", "type": "string"},  # 参考基因组类型
            # {"name": "bam_list", "type": "string"},  # 如果是gatk 样本+建库类型+bam文件，否则样本bam列表
            {"name": "in_bam", "type": "infile", "format": "align.bwa.bam_dir"},  # bam格式文件
            {"name": "bam_list", "type": "string"},  # 如果是gatk 样本+建库类型+bam文件，否则样本bam列表
            # {"name": "ref_dict", "type": "string"},  # 参考基因组的配置文件，里面可以解析出有多少条染色体
            {"name": "call_type", "type": "string", "default": "sentieon"},  # call snp的方式
            {"name": "des", "type": "string"},
            {"name": "des_type", "type": "string"},
            {"name": "ref_gtf", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "align_method", "type": "string", "default": "hisat"},
            {"name": "analysis_format", "type": "string", "default": "bam"},
            {"name": "algorithm", "type": "string", "default": "HaplotypeCaller"}
        ]
        self.add_option(options)
        self.config = Config()
        self.logger.info("in_bam:{}".format(self.option("in_bam")))

        global WORK_DIR
        WORK_DIR = self.work_dir

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(WholeSnpModule, self).run()
        if self.option("call_type").lower() == "sentieon":
            self.run_sentieon()
        elif self.option("call_type").lower() == "samtools":
            self.run_samtools()
        else:
            self.run_gatk()

    def run_sentieon(self):
        self.snp=self.add_module('medical_transcriptome.snp.call_snp_indel_sentieon')
        bam_list = self.work_dir + "/bamlist"
        with open(bam_list, "w") as w:
            for f in sorted(os.listdir(self.option("in_bam").prop["path"])):
                f_path = self.option("in_bam").prop["path"] + "/" + f
                w.write(f_path + "\n")
        opts = {
            'ref_fasta': self.option("ref_genome_custom"),
            'bam_list': bam_list,
            'ref_gtf': self.option("ref_gtf"),
            'des':  self.option("des"),
            'des_type': self.option('des_type'),
            'analysis_format': self.option('analysis_format'),
            'align_method': self.option('align_method'),
            "algorithm": self.option("algorithm")
          }
        self.snp.set_options(opts)
        self.snp.on('end', self.set_output, 'sentieon')
        self.snp.run()

    def run_samtools(self):
        self.snp = self.add_module('medical_transcriptome.snp.call_snp_indel_samtools')
        opts = {
            'ref_genome_custom': self.option("ref_genome_custom"),
            'in_bam': self.option("in_bam"),
            'ref_gtf': self.option("ref_gtf"),
            'des': self.option("des"),
            'des_type': self.option('des_type')
        }
        self.snp.set_options(opts)
        self.snp.on('end', self.set_output, 'samtools')
        self.snp.run()

    def run_gatk(self):
        self.snp = self.add_module('medical_transcriptome.snp.call_snp_indel_gatk')
        opts = {
            'ref_genome_custom': self.option("ref_genome_custom"),
            "ref_genome": "customer_mode",
            'in_bam': self.option("in_bam"),
            'ref_gtf': self.option("ref_gtf"),
            'des': self.option("des"),
            'des_type': self.option('des_type')
        }
        self.snp.set_options(opts)
        self.snp.on('end', self.set_output, 'gatk')
        self.snp.run()

    def set_output(self, event):
        self.logger.info("set output started!!!")
        obj = event["bind_object"]
        if event['data'] == 'sentieon':
            self.logger.info("llllllllllllooking for event data")
            self.logger.info("llllllllllllooking for event data{}".format(obj.output_dir))
            self.logger.info("llllllllllllooking for event data{}".format(os.path.join(obj.output_dir,"predeal")))
            for file in os.listdir(os.path.join(obj.output_dir,"predeal")):
                if file.endswith(".xls") or file.endswith("info"):
                    old = os.path.join(obj.output_dir,"predeal", file)
                    new = os.path.join(self.output_dir, file)
                    if os.path.exists(new):
                        os.remove(new)
                    os.link(old, new)
            self.logger.info("set output done")
            self.end()
        else:
            self.logger.info("llllllllllllooking for event data")
            for file in os.listdir(obj.output_dir):
                if file.endswith(".xls") or file.endswith("info"):
                    old = os.path.join(obj.output_dir, file)
                    new = os.path.join(self.output_dir, file)
                    if os.path.exists(new):
                        os.remove(new)
                    os.link(old, new)
            self.logger.info("set output done")
            self.end()

    def end(self):
        super(WholeSnpModule, self).end()




class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "medical_transcriptome.snp.whole_snp",
            "instant": False,
            "options": dict(
                # ref_dict=test_dir + "/" + "Mus_musculus.GRCm38.dna_rm.toplevel.clean.dict",
                # bam_list="/mnt/ilustre/users/sanger-dev/workspace/20190703/Single_bam_realign5598/BamRealign/output/bam.list",
                # bam_list="/mnt/ilustre/users/sanger-dev/workspace/20190722/Single_bam_realign5942/BamRealign/output/bam.list",
                call_type="samtools",
                in_bam="/mnt/ilustre/users/sanger-dev/workspace/20200825/Snp_medical_transcriptome_3735_2148/bam_folder/",
                ref_genome_custom="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/biomart/biomart.txt",
                des_type="type1",
                ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf",
                bam_list = "/mnt/ilustre/users/sanger-dev/workspace/20200825/Snp_medical_transcriptome_3735_2148/bamlist_new",
                # scm="complete",
                # scd="correlation",
                # corr_method='pearson',
                # output=None,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()