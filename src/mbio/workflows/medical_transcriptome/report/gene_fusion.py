# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import glob
import pandas as pd
import os
import re
from mbio.packages.ref_rna_v3.gene_fusion.extract_pos_bygtf import extract_all_gene_pos,extract_all_chr_length,gtf_check
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.config import Config
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart import Chart


class GeneFusionWorkflow(Workflow):
    """
    基因融合分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GeneFusionWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'sample_list_str', 'type': 'string'},
            dict(name='bamlist', type='string'),
            dict(name="ref_gtf", type="infile", format="ref_rna_v2.common"),
            dict(name="ref_genome_custom", type="infile", format="ref_rna_v2.common"),
            dict(name="gene_fusion_main_id", type="string"),
            dict(name="min_junction_reads", type="int", default="1"),
            dict(name="min_sum_frags", type="int", default="2"),
            dict(name="min_novel_junction_support", type="int", default="3"),
            dict(name="min_spanning_frags_only", type="int", default="5"),
            dict(name="min_FFPM", type="float", default="0.1"),
            dict(name="species", type="string", default=None),
            dict(name="task_id", type="string"),
            dict(name="group_dict", type="infile", format="ref_rna_v2.common"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.gene_fusion = self.add_module("medical_transcriptome.gene_fusion")
        self.circos = None
        # database = Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]
        # collection = database['sg_annotation_query']
        # anno_info =  collection.find_one({'task_id': self.option('task_id')})
        self.assemble_level_file = ""
        self.chr_length_path = os.path.join(os.path.dirname(self.option("ref_genome_custom").prop["path"]),"star_27","ctat_genome_lib_build_dir","ref_genome.fa.star.idx","chrNameLength.txt")
        # if "flag" in anno_info:
        #     self.id_modify = "yes"
        # else:
        #     self.id_modify = None
        self.gtf = os.path.join(self.work_dir,"final.gtf")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/03 Gene_structure_analysis/03 GeneFusion')
        self.inter_dirs = []


    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(GeneFusionWorkflow, self).send_log(data)

    def run(self):
        self.gene_fusion.on("end", self.set_db)
        self.logger.info("在基因融合分析之前")
        self.get_run_log()
        self.run_star_fusion()
        self.logger.info("在基因融合分析之后")
        super(GeneFusionWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_gene_fusion", main_id=self.option('gene_fusion_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.api_fusion = self.api.api("medical_transcriptome.gene_fusion")
        # add result info
        s3_output = os.path.join(self._sheet.output,"star_fusion")
        if self.circos:
            extract_all_gene_pos(self.option("ref_gtf").prop["path"],os.path.join(self.work_dir,"gene_pos"))
            extract_all_chr_length(self.assemble_level_file,self.chr_length_path,os.path.join(self.work_dir,"chr_length"))
            pos_file = os.path.join(self.work_dir,"gene_pos")
            chr_length_file = os.path.join(self.work_dir,"chr_length")
            fusion_detail_anno = self.gene_fusion.output_dir
            self.api_fusion.add_fusion_main(fusion_result=fusion_detail_anno,pos_file=pos_file,chr_length = chr_length_file,circos= self.circos,main_id=self.option('gene_fusion_main_id'),s3_output = s3_output, group=self.option('sample_list_str'))
        else:
            fusion_detail_anno = self.gene_fusion.output_dir
            self.api_fusion.add_fusion_main(fusion_result=fusion_detail_anno, main_id=self.option('gene_fusion_main_id'),s3_output = s3_output, group=self.option('sample_list_str'))
        self.end()
        #os.link(os.path.join(snp_anno, "data_anno_pre.xls"), os.path.join(self.output_dir, "snp_anno.xls"))

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        sample_names = os.listdir(os.path.join(self.gene_fusion.output_dir,"star_fusion"))
        genefusion = {"result_dir" :self.gene_fusion.output_dir+"/star_fusion/{sample_name}/star-fusion.fusion_predictions.abridged.tsv",
                        "chr_length_path":os.path.join(self.work_dir,"chr_length"),
                        "gene_pos":os.path.join(self.work_dir,"gene_pos")}
        fusion_files = [genefusion["result_dir"].format(sample_name=sample) for sample in sample_names]
        chart.chart_gene_fusion_ciros(sample_names, fusion_files, genefusion)
        chart.to_pdf()
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.gene_fusion.output_dir + "/" + os.path.basename(p))


    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.gene_fusion.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.gene_fusion.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.gene_fusion.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.gene_fusion.output_dir)
        self.inter_dirs = [
            ["03 Gene_structure_analysis", "", "基因结构分析数据挖掘结果目录", 0],
            ["03 Gene_structure_analysis/03 GeneFusion", "", "基因融合结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            [r'.', '', '基因融合分析文件', 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        result_dir.add_regexp_rules([
            [r'fusion_stat.txt', 'txt', '基因融合信息统计表', 0],
            [r'.*\.pdf', 'pdf', '基因融合circos图', 0],
            [r'star_fusion', 'dir', '基因融合信息详情', 0],
            [r'star_fusion/.*','dir', '单个样本融合信息详情', 0],
            [r'star_fusion/.*/star-fusion.fusion_predictions.abridged.tsv', 'tsv', '基因融合信息详情缩减表', 0],
            [r'star_fusion/.*/star-fusion.fusion_predictions.tsv', 'tsv', '基因融合信息详情detail表', 0],
            [r'star_fusion/.*/fusion_inspector', 'dir', 'fusion_inspector文件', 0],
            [r'star_fusion/.*/fusion_inspector/output', 'dir', 'fusion可视化所用文件', 0],
            [r'star_fusion/.*/fusion_inspector/output/finspector.fa', '', 'fusion序列文件', 0],
            [r'star_fusion/.*/fusion_inspector/output/finspector.bed', '', 'fusion注释文件', 0],
            [r'star_fusion/.*/fusion_inspector/output/finspector.junction_reads.bam', '', 'fusion bam文件', 0],
            [r'star_fusion/.*/fusion_inspector/output/finspector.spanning_reads.bam', '', 'fusion bam文件', 0],
        ])

        super(GeneFusionWorkflow, self).end()


    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find', code="15600505")
        return to_path


    def run_star_fusion(self):
        if os.path.exists(self.option("ref_genome_custom").prop['path']):
            ref_genome_custom1 = self.option("ref_genome_custom").prop['path']
        elif re.match(r'/mnt/ilustre/', self.option("ref_genome_custom").prop['path']):
            ref_genome_custom1 = self.option("ref_genome_custom").prop['path'].replace('ilustre','lustre')
        else:
            self.set_error("ref genome file not exist", code="15600506")
        if os.path.exists(self.option("ref_gtf").prop['path']):
            ref_gtf1 = self.option("ref_gtf").prop['path']
        elif re.match(r'/mnt/ilustre/', self.option("ref_gtf").prop['path']):
            ref_gtf1 = self.option("ref_gtf").prop['path'].replace('ilustre','lustre')
        else:
            self.set_error("ref gtf file not exist", code="15600507")
        if os.path.exists(os.path.join(os.path.dirname(ref_gtf1),"assembly_level.txt")):
            self.circos = "yes"
            self.assemble_level_file = os.path.join(os.path.dirname(ref_gtf1),"assembly_level.txt")

        if not os.path.exists(self.work_dir + "/bam_folder/"):
            os.mkdir(self.work_dir + "/bam_folder/")
        with open(self.option("bamlist")) as f, open(self.work_dir + '/bamlist_new', 'w') as bam_w:
            for line in f:
                bam_name = os.path.basename(line.strip())
                if re.match(r'^\w+://\S+/.+$', line.strip()) or re.match(r'/mnt/ilustre', line.strip()):
                    transfer = MultiFileTransfer()
                    transfer.add_download(line.strip(), self.work_dir + "/bam_folder/")
                    transfer.perform()
                    #path = self.download_s3_file(line.strip(), bam_name)
                    bam_w.write(self.work_dir + "/bam_folder/" + bam_name + '\n')
                    # if os.path.exists(self.work_dir + "/bam_folder/" + bam_name):
                    #     os.remove(self.work_dir + "/bam_folder/" + bam_name)
                    # os.link(path, self.work_dir + "/bam_folder/" + bam_name)
                else:
                    if os.path.exists(self.work_dir + "/bam_folder/" + bam_name):
                        os.remove(self.work_dir + "/bam_folder/" + bam_name)
                    os.link(line.strip(), self.work_dir + "/bam_folder/" + bam_name)
                    bam_w.write(self.work_dir + "/bam_folder/" + bam_name + '\n')

        gtf_check(ref_gtf1,self.gtf)
        opts = {
                "ref_genome_custom": ref_genome_custom1,
                "ref_gtf": self.gtf,
                "min_junction_reads":self.option("min_junction_reads"),
                "min_sum_frags": self.option("min_sum_frags"),
                "min_novel_junction_support": self.option("min_novel_junction_support"),
                "min_spanning_frags_only": self.option("min_spanning_frags_only"),
                "min_FFPM": self.option("min_FFPM"),
                "bam_list" :self.work_dir + '/bamlist_new',
                "species": self.option("species"),
                "task_id": self.option("task_id")
        }
        if self.circos:
            opts.update({"circos":"yes"})
        # if self.id_modify:
        #     opts.update({"id_modify": "yes"})
        self.gene_fusion.set_options(opts)
        self.logger.info("star_fusion运行开始")
        self.gene_fusion.run()
        self.logger.info("star_fusion运行结束")

