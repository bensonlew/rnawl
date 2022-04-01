# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20200609

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
import gevent.subprocess as subprocess
import pandas as pd
from biocluster.config import Config


class StarFusionModule(Module):
    """
    该module从bam文件开始，经过鉴定过滤替换等多个步骤，最终得到结果文件
    """
    def __init__(self, work_id):
        super(StarFusionModule, self).__init__(work_id)
        options = [
            {"name": "ref_genome_custom", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "min_junction_reads", "type": "int", "default": 1},
            # 最小junction_reads数  minimum fusion support = ( # junction_reads + # spanning_frags ) Default: 2   最小融合支持：junction+spanning_frags
            {"name": "min_sum_frags", "type": "int", "default": 2},
            # (minimum of junction reads required if breakpoint  lacks involvement of only reference junctions) 当breakpoint  不在reference junction时，junction的最小数值要求
            {"name": "min_novel_junction_support", "type": "int", "default": 3},
            # minimum number of rna-seq fragments required as fusion evidence if there are no junction reads (default: 5) 如果没有junction_reads只有spanning_frags需要多少的支持才进行考虑
            {"name": "min_spanning_frags_only", "type": "int", "default": 5},
            # minimum FFPM (fusion fragments per million rna-seq frags)  (default: 0.1)
            {"name": "min_FFPM", "type": "float", "default": 0.1},
            {"name": "ref_gtf", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "in_bam", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "min_FFPM", "type": "float", "default": 0.1},
            {"name": "circos", "type": "string", "default": None},
            {"name": "id_modify", "type": "string", "default": None},
            {"name": "sample_name", "type": "string", "default": None},
            {"name": "task_id", "type": "string"},
            {"name": "species", "type": "string", "default": None}
        ]
        self.add_option(options)
        self.ctat_genome_lib_build_dir = ""
        if self.option("species") == "Homo_sapiens":
            self.dfam_path = Config().SOFTWARE_DIR + "/database/DFAM_3.1/human/homo_sapiens_dfam.hmm"
        elif self.option("species") == "Mus_musculus":
            self.dfam_path = Config().SOFTWARE_DIR + "/database/DFAM_3.1/mouse/mus_musculus_dfam.hmm"
        else:
            self.dfam_path =  Config().SOFTWARE_DIR + "/database/DFAM_3.1/common/Dfam.hmm"
        self.pfam_path = Config().SOFTWARE_DIR + "/database/Annotation/other2019/pfam32/Pfam-A.hmm"
        self.make_lib = None
        self.star_fusion = self.add_tool("ref_rna_v3.gene_fusion.star_fusion")
        self.bam2fastq = self.add_tool("ref_rna_v3.gene_fusion.bam2fastq")
        self.result_filter = self.add_tool("ref_rna_v3.gene_fusion.fusion_result_filter")
        self.id_modify = self.add_tool("ref_rna_v3.gene_fusion.id_modify")
        self.fusioninspector = self.add_tool("ref_rna_v3.gene_fusion.fusioninspector")
        self.fusion_result = ""
        self.fusion_result_detail = ""
        self.final_result_path =""
        self.final_result_detail_path = ""


    def check_options(self):
        if not self.option("ref_genome_custom").is_set:
            raise OptionError("缺少ref_genome_custom参数")
        if not self.option("in_bam"):
            raise OptionError("请输入bam文件")
        if not self.option("ref_gtf"):
            raise OptionError("缺少ref_gtf参数")
        return True


    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()


    def check_lib(self):
        genome_path = self.option("ref_genome_custom").prop["path"]
        if os.path.exists(os.path.join(os.path.dirname(genome_path),"star_27")):
            self.logger.info("{}已建立对应的star2.7库并构建融合基因数据库,直接进行分析".format(self.option("species")))
            self.ctat_genome_lib_build_dir = os.path.join(os.path.dirname(genome_path),"star_27","ctat_genome_lib_build_dir")
            self.logger.info("数据库位置是{}".format(self.ctat_genome_lib_build_dir))
            self.run_bam2fastq()
        else:
            self.ctat_genome_lib_build_dir = os.path.join(os.path.dirname(genome_path), "star_27",
                                                          "ctat_genome_lib_build_dir")
            self.make_lib=self.add_tool("ref_rna_v3.gene_fusion.make_lib")
            self.make_lib.set_options({
                "ref_fasta":self.option("ref_genome_custom").prop["path"],
                "gtf_path" :self.option("ref_gtf").prop["path"],
                "dfam_path" : self.dfam_path ,
                "pfam_path": self.pfam_path
            })
            self.make_lib.on('end', self.run_bam2fastq)
            self.make_lib.run()

    def run_bam2fastq(self):
        self.logger.info("开始将bam文件转化成fastq文件")
        opts = {
            "in_bam": self.option("in_bam").prop["path"],
        }
        self.bam2fastq.set_options(opts)
        self.bam2fastq.on('end',self.run_star_fusion)
        self.bam2fastq.run()


    def run_star_fusion(self):
        self.logger.info("完成将bam文件转化成fastq文件")
        self.logger.info("开始运行star_fusion")
        opts = {
            "lib_path": self.ctat_genome_lib_build_dir,
            "fastq_l": self.bam2fastq.option("out_fq1"),
            "fastq_r": self.bam2fastq.option("out_fq2"),
            "min_junction_reads": self.option("min_junction_reads"),
            "min_sum_frags": self.option("min_sum_frags"),
            "min_novel_junction_support": self.option("min_novel_junction_support"),
            "min_spanning_frags_only": self.option("min_spanning_frags_only"),
            "min_FFPM": self.option("min_FFPM"),
        }
        self.star_fusion.set_options(opts)
        self.star_fusion.on("end",self.run_filter_scafford)
        self.star_fusion.run()

    def run_filter_scafford(self):
        self.logger.info("完成将star_fusion的运行")
        self.logger.info("开始对star_fusion结果根据组装水平进行过滤")
        fusion_result = os.path.join(self.star_fusion.output_dir,"star-fusion.fusion_predictions.abridged.tsv")
        fusion_result_detail = os.path.join(self.star_fusion.output_dir,"star-fusion.fusion_predictions.tsv")
        # assembly_level_file = os.path.join(os.path.dirname(self.option("ref_gtf").prop["path"]),"assembly_level.txt")
        assembly_level_file = os.path.join(os.path.dirname(self.option("ref_genome_custom").prop["path"]),"..","gtf" ,"assembly_level.txt")
        opts = {
            "fusion_result": fusion_result,
            "fusion_result_detail": fusion_result_detail
        }
        if self.option("circos"):
            opts.update({
                "circos":self.option("circos"),
                "assembly_level_file":assembly_level_file
            })
        self.result_filter.set_options(opts)
        self.result_filter.on("end",self.run_id_modify)
        self.result_filter.run()

    def run_id_modify(self):
        self.logger.info("完成对star_fusion结果根据组装水平进行过滤")
        self.logger.info("开始将检查结果是否需要用客户上传的gene_name进行替代")
        opts={
            "fusion_result":self.result_filter.option("fusion_result_out"),
            "fusion_result_detail":self.result_filter.option("fusion_result_detail_out"),
            "fastq_l": self.bam2fastq.option("out_fq1"),
            "fastq_r": self.bam2fastq.option("out_fq2")
        }
        if self.option("id_modify"):
            opts.update({
                "id_modify": self.option("id_modify"),
                "task_id":self.option("task_id")
            })
        self.id_modify.set_options(opts)
        self.id_modify.on("end", self.run_fusioninspector)
        self.id_modify.run()

    def run_fusioninspector(self):
        self.logger.info("完成将检查结果是否需要用客户上传的gene_name进行替代")
        self.logger.info("开始对分析结果进行检验和可视化")
        opts={
            "lib_path":self.ctat_genome_lib_build_dir,
            "fusion_result":self.id_modify.option("fusion_result_detail_out"),
            "left_fq":self.id_modify.option("evidence_fq1"),
            "right_fq": self.id_modify.option("evidence_fq2")
        }
        self.fusioninspector.set_options(opts)
        self.fusioninspector.on("end", self.final_adjust)
        self.fusioninspector.run()

    def final_adjust(self):
        self.logger.info("完成igv可视化操作")
        self.logger.info("开始最终校正,去除LeftGene/RightGene中的name^gene_id中的gene_name去掉")
        #将最终统计结果文件star-fusion.fusion_predictions.abridged.tsv中的LeftGene/RightGene中的name^gene_id中的gene_name去掉，用以最后导表
        final_result = pd.read_table(self.id_modify.option("fusion_result_out").prop["path"])
        final_result["LeftGene"] = final_result["LeftGene"].apply(lambda x: x.split("^")[-1])
        final_result["RightGene"] = final_result["RightGene"].apply(lambda x: x.split("^")[-1])
        if self.option("circos"):
            final_result["ciros"] = 1 if final_result.shape[0] >=2 else 0
        else:
            final_result["ciros"] = 0
        self.final_result_path = os.path.join(self.work_dir,"star-fusion.fusion_predictions.abridged.tsv")
        final_result.to_csv(self.final_result_path,index = False,sep = "\t")

        # 将最终统计结果文件fusion_result_detail_out中的LeftGene/RightGene中的name^gene_id中的gene_name去掉，用以最后导表
        final_result_detail = pd.read_table(self.id_modify.option("fusion_result_detail_out").prop["path"])
        final_result_detail["LeftGene"] = final_result["LeftGene"].apply(lambda x: x.split("^")[-1])
        final_result_detail["RightGene"] = final_result["RightGene"].apply(lambda x: x.split("^")[-1])
        if self.option("circos"):
            final_result_detail["ciros"] = 1 if final_result_detail.shape[0] >=2 else 0
        else:
            final_result_detail["ciros"] = 0
        self.final_result_detail_path = os.path.join(self.work_dir, "star-fusion.fusion_predictions.tsv")
        final_result_detail.to_csv(self.final_result_detail_path, index=False, sep="\t")
        self.set_output()

    def set_output(self):
        result_dir = os.path.join(self.output_dir,self.option("sample_name"))
        if os.path.exists(result_dir):
            os.system('rm -r %s' % result_dir)
        os.mkdir(result_dir)
        final_fusion_result = os.path.join(result_dir,"star-fusion.fusion_predictions.tsv")
        if os.path.isfile(final_fusion_result):
            os.remove(final_fusion_result)
        final_fusion_result_detail = os.path.join(result_dir,"star-fusion.fusion_predictions.abridged.tsv")
        if os.path.isfile(final_fusion_result_detail):
            os.remove(final_fusion_result_detail)
        os.link(self.final_result_path,final_fusion_result_detail)
        os.link(self.final_result_detail_path, final_fusion_result)
        fusioninspector_dir = os.path.join(result_dir,"fusion_inspector")
        if os.path.exists(fusioninspector_dir):
            os.system('rm -r %s' % fusioninspector_dir)
        os.mkdir(fusioninspector_dir)
        os.system('cp -r %s %s' % (self.fusioninspector.output_dir, fusioninspector_dir))
        self.end()





    def run(self):
        super(StarFusionModule, self).run()
        self.logger.info("开始检查是否已完成star2.7索引建立以及融合库")
        self.check_lib()


    def end(self):
        super(StarFusionModule, self).end()



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
            "id": "star_fusion_test" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "ref_rna_v3.star_fusion",
            "instant": False,
            "options": dict(
                ref_genome_custom="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Pisum_sativum/v1a_v1a/dna/Pisum_sativum_v1a.fa",
                ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Pisum_sativum/v1a_v1a/gtf/final.gtf",
                in_bam = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/data_pre/bam/bam/A1.bam",
                task_id = "ref_rna_v3.1_test",
                sample_name = "A1"

            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()