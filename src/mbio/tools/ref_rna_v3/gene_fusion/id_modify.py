# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modifiy:2020.06.08

import os
import re
import shutil
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from biocluster.config import Config



class IdModifyAgent(Agent):
    """
    这里对于客户上传了自定义gene_name/description的项目,将结果表和结果文件中的name进行替换

    """

    def __init__(self, parent):
        super(IdModifyAgent, self).__init__(parent)
        options = [
            # 原始fastq r1端文件(用来提取evidence_reads_1)
            {"name": "fastq_l", "type": "infile", "format": "ref_rna_v2.common"},
            # 原始fastq r2端文件(用来提取evidence_reads_2)
            {"name": "fastq_r", "type": "infile", "format": "ref_rna_v2.common"},
            #项目的task_id
            {"name": "task_id", "type": "string"},
            {"name": "id_modify", "type": "string", "default": None},
            # star_fusion的结果（简单)  star-fusion.fusion_predictions.abridged.tsv
            {"name":"fusion_result","type": "infile", "format": "ref_rna_v2.common"},
            # star_fusion的结果(具体）  star-fusion.fusion_predictions.tsv
            {"name": "fusion_result_detail", "type": "infile", "format": "ref_rna_v2.common"},
            #更名后的star_fusion结果(简单)
            {"name": "fusion_result_out", "type": "outfile", "format": "ref_rna_v2.common"},
            # 更名后的star_fusion的结果(具体）
            {"name": "fusion_result_detail_out", "type": "outfile", "format": "ref_rna_v2.common"},
            #根据更名后的表格生成的新的fq1(evidence fq1)
            {"name": "evidence_fq1", "type": "outfile", "format": "ref_rna_v2.common"},
            # 根据更名后的表格生成的新的fq1(evidence fq2)
            {"name": "evidence_fq2", "type": "outfile", "format": "ref_rna_v2.common"},
        ]
        self.add_option(options)
        self._memory_increase_step = 50
        self.step.add_steps('id_modify')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.id_modify.start()
        self.step.update()

    def step_end(self):
        self.step.id_modify.finish()
        self.step.update()

    def check_options(self):

        if not self.option("fusion_result").is_set:
            raise OptionError("请输入star_fusion的结果文件")
        if not self.option("fusion_result_detail").is_set:
            raise OptionError("请输入star_fusion的结果详情文件")
        # if not self.option("fastq_l").is_set:
        #     raise OptionError("请输入用于分析的left_fq文件!")
        # if not self.option("fastq_r").is_set:
        #     raise OptionError("请输入用于分析的right_fq文件!")




    def set_resource(self):
        self._cpu = 5
        self._memory = '15G'

    def end(self):
        super(IdModifyAgent, self).end()


class IdModifyTool(Tool):

    def __init__(self, config):
        super(IdModifyTool, self).__init__(config)
        # self.star_fusion_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/bin"
        # self.set_environ(PATH=self.star_fusion_path)
        # self.make_lib_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/lib/STAR-Fusion/ctat-genome-lib-builder/"
        # self.fusion_inspect_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/lib/STAR-Fusion/FusionInspector/"
        # self.get_fq_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/lib/STAR-Fusion/util/"
        # self.fusion_result_out = ""
        # self.fusion_result_detail_out = self.option("fusion_result_detail").prop["path"]

        self.star_fusion_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v3/gene_fusion/miniconda3/bin/"
        self.set_environ(PATH=self.star_fusion_path)
        self.make_lib_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v3/gene_fusion/miniconda3/lib/STAR-Fusion/ctat-genome-lib-builder/"
        self.fusion_inspect_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v3/gene_fusion/miniconda3/lib/STAR-Fusion/FusionInspector/"
        self.get_fq_path = self.config.SOFTWARE_DIR +"/bioinfo/ref_rna_v3/gene_fusion/miniconda3/lib/STAR-Fusion/util/"
        self.fusion_result_out = ""
        self.fusion_result_detail_out = self.option("fusion_result_detail").prop["path"]

    def prepare_annotation_file(self):
        """
        准备注释文件
        """
        self.logger.info("开始运行连接mongo库并获取gene,id/name文件")
        project_type = 'ref_rna_v2'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        task_id = self.option("task_id")
        annot_table = db['sg_annotation_query']
        annot_main = annot_table.find_one({"task_id": task_id})
        if not annot_main:
            self.set_error("Not Found in sg_annotation_query by query %s", variables=(task_id),
                                     )
        if "main_id" not in annot_main:
            annot_main_id = annot_main['_id']
        else:
            annot_main_id = annot_main['main_id']
        annot_detail = db['sg_annotation_query_detail']
        query_dict = dict(query_id=annot_main_id, )
        result_dict = dict(_id=0, gene_name=1, gene_id=1,gene_name_old=1,)
        result = annot_detail.find(query_dict, result_dict)
        result_pd = pd.DataFrame(list(result))
        result_pd.set_index("gene_id", inplace=True)
        result_pd = result_pd[~result_pd.index.duplicated(keep='first')]
        # result_pd = result_pd.loc[:, ["gene_name_old","gene_name"]]
        # result_pd.columns = ["gene_name_old","gene_name"]
        result_pd = result_pd.loc[:, ["gene_name"]]
        result_pd.columns = [ "gene_name"]
        gene_annot = os.path.join(self.work_dir, "seq_annot.xls")
        result_pd.to_csv(gene_annot, sep='\t', header=True, index=True)
        os.system(
            r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)

    def gene_name_replace(self):
        seq_detail=pd.read_table("seq_annot.xls")
        seq_detail.set_index("gene_id", inplace=True)
        seq_detail.fillna("-")
        seq_detail_index = seq_detail.to_dict(orient='index')
        def change_fusion_name(series, ndict):
            raw_name = series["#FusionName"]
            raw_left_name = raw_name.split("-")[0]
            raw_right_name = raw_name.split("-")[1]
            left_gene = series["LeftGene"].split("^")[-1]
            right_gene = series["RightGene"].split("^")[-1]
            new_left_gene = ndict[left_gene]["gene_name"] if ndict[left_gene]["gene_name"] != "-" else raw_left_name
            new_right_gene = ndict[right_gene]["gene_name"] if ndict[right_gene]["gene_name"] != "-" else raw_right_name
            return new_left_gene + "-" + new_right_gene
        furesult = pd.read_table(self.option("fusion_result").prop["path"])
        furesult["#FusionName"] = furesult.apply(change_fusion_name, args=(seq_detail_index,), axis=1)
        furesult_detail= pd.read_table(self.option("fusion_result_detail").prop["path"])
        furesult_detail["#FusionName"] = furesult_detail.apply(change_fusion_name, args=(seq_detail_index,), axis=1)
        self.fusion_result_out = os.path.join(self.output_dir,os.path.basename(self.option("fusion_result").prop["path"]))
        furesult.to_csv(self.fusion_result_out,index=False,sep="\t")
        self.fusion_result_detail_out = os.path.join(self.output_dir, os.path.basename(self.option("fusion_result_detail").prop["path"]))
        furesult_detail.to_csv(self.fusion_result_detail_out,index=False,sep="\t")


    def extract_evidence_reads(self):
        cmd = "{}get_FUSION_EVIDENCE_fastqs.pl ".format(self.get_fq_path)
        cmd += "--fusions {} ".format(self.fusion_result_detail_out)
        cmd += "--output_prefix star-fusion "
        cmd += "--left_fq {} ".format(self.option("fastq_l").prop["path"])
        cmd += "--right_fq {} ".format(self.option("fastq_r").prop["path"])
        self.logger.info("开始运行get_FUSION_EVIDENCE_fastqs.pl")
        command = self.add_command("get_fastq", cmd, ignore_error=True, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用开始运行get_FUSION_EVIDENCE_fastqs对新更换后的star_fusion结果生成evidence_fastq文件成功!")
        elif command.return_code in [1, -9]:  # 当返回码为1或-9，加内存重试
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用开始运行get_FUSION_EVIDENCE_fastqs对新更换后的star_fusion结果生成evidence_fastq文件失败!")

    def set_output(self):
        soure_fq1 = os.path.join(self.work_dir, "star-fusion.fusion_evidence_reads_1.fq")
        out_fq1 = os.path.join(self.output_dir, "star-fusion.fusion_evidence_reads_1.fq")
        if os.path.exists(out_fq1):
            os.remove(out_fq1)
        os.link(soure_fq1, out_fq1)
        soure_fq2 = os.path.join(self.work_dir, "star-fusion.fusion_evidence_reads_2.fq")
        out_fq2 = os.path.join(self.output_dir, "star-fusion.fusion_evidence_reads_2.fq")
        if os.path.exists(out_fq2):
            os.remove(out_fq2)
        os.link(soure_fq2, out_fq2)
        self.option("evidence_fq1").set_path(out_fq1)
        self.option("evidence_fq2").set_path(out_fq2)
        if self.option("id_modify"):
            self.option("fusion_result_out").set_path(self.fusion_result_out)
            self.option("fusion_result_detail_out").set_path(self.fusion_result_detail_out)
        else:
            self.option("fusion_result_out").set_path(self.option("fusion_result").prop["path"])
            self.option("fusion_result_detail_out").set_path(self.option("fusion_result_detail").prop["path"])
        self.logger.info("最终evidence_fq1是{}".format(out_fq1))
        self.logger.info("最终evidence_fq2是{}".format(out_fq2))
        self.logger.info("最终fusion_result_out是{}".format(self.option("fusion_result").prop["path"]))
        self.logger.info("最终fusion_result_detail_out是{}".format(self.option("fusion_result_detail").prop["path"]))




    def run(self):
        """
        运行
        """
        super(IdModifyTool, self).run()
        self.logger.info("运行开始")
        if self.option("id_modify"):
            self.prepare_annotation_file()
            self.gene_name_replace()
        self.extract_evidence_reads()
        self.set_output()
        self.logger.info("运行结束")
        self.end()





class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.gene_fusion.id_modify",
            "instant": False,
            "options": dict(
                # lib_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/data_pre/make_lib/ctat_genome_lib_build_dir",
                # fastq_l="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/data_pre/fastq_raw/fastq/A1.clean.1.fastq",
                # fastq_r="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/data_pre/fastq_raw/fastq/A1.clean.2.fastq",
                # fusion_result="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/test_raw_fastq/star-fusion.fusion_predictions.abridged.tsv",
                # fusion_result_detail = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/test_raw_fastq/star-fusion.fusion_predictions.tsv",
                fastq_l="/mnt/ilustre/users/sanger-dev/workspace/20200731/GeneFusion_tsg_38132_1089_4587/GeneFusion/StarFusion/Bam2fastq/output/Cav_1_KO_1.fq1",
                fastq_r="/mnt/ilustre/users/sanger-dev/workspace/20200731/GeneFusion_tsg_38132_1089_4587/GeneFusion/StarFusion/Bam2fastq/output/Cav_1_KO_1.fq2",
                fusion_result="/mnt/ilustre/users/sanger-dev/workspace/20200731/GeneFusion_tsg_38132_1089_4587/GeneFusion/StarFusion/FusionResultFilter/output/star-fusion.fusion_predictions.abridged.tsv",
                fusion_result_detail="/mnt/ilustre/users/sanger-dev/workspace/20200731/GeneFusion_tsg_38132_1089_4587/GeneFusion/StarFusion/FusionResultFilter/output/star-fusion.fusion_predictions.tsv",
                task_id = "tsg_38132",
                id_modify = "1"

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
