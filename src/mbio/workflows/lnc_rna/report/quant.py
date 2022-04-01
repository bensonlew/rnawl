# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import pandas as pd
import glob
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.config import Config
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder

class QuantWorkflow(Workflow):
    """
    quant
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(QuantWorkflow, self).__init__(wsheet_object)
        options = [
            #dict(name="transcriptome", type="infile", format="lnc_rna.fasta"),
            #dict(name="fastq", type="infile", format="lnc_rna.fastq_list"),
            dict(name="method", type="string", default="rsem"),
            dict(name="libtype", type="string", default=None),
            #dict(name="t2g", type="infile", format="lnc_rna.common"),
            dict(name="read_len", type="int", default=149),
            dict(name="transcript_main_id", type="string"),
            dict(name="gene_main_id", type="string"),
            dict(name="task_type", type="string"),
            dict(name="raw_task_id", type="string"),
            dict(name="submit_location", type="string"),
            dict(name="exp_type", type="string"),
            dict(name="task_id", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.tool = self.add_tool("ref_rna_v2.quant")
        self.tool = self.add_module("lnc_rna.quant")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/01 Express/01 Exp_Analysis')
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
        super(QuantWorkflow, self).send_log(data)

    def run(self):
        self.get_run_log()
        #db = Config().get_mongo_client(mtype="lnc_rna")[Config().get_mongo_dbname("lnc_rna")]
        #col = db["sg_task"]
        #result = col.find_one({"task_id": self.option('raw_task_id')})
        #basepath = result["fastq"].split("/fq_list.txt")[0]
        #with open (self.option("fastq").prop["path"], "r") as f, open(self.work_dir + "/fq_list.txt", "w") as w:
         #   for line in f:
         #       items = line.strip().split("\t")
          #      if (r'workspace', line):
           #         if len(items) == 3:
            #            w.write(items[0] + "\t" + basepath + "/" + os.path.basename(items[1]) + "\t" + basepath + "/" + os.path.basename(items[2]) + "\n")
             #       else:
                #        w.write(items[0] + "\t" + basepath + "/" + os.path.basename(items[1]) + "\n")
              #  else:
               #     w.writeline(line)
        #self.option("fastq", self.work_dir + "/fq_list.txt")
        if self.option('exp_type').lower() == 'fpkm' and self.option('method').lower() == 'rsem':
            self.start_listener()
            self.fire("start")
            self.dump_fpkm(self.option('raw_task_id'),
                           self.option('transcript_main_id'),
                           self.option('gene_main_id'))
            # self.fire("end")
            # super(QuantWorkflow, self).end()
            self.end()
        elif self.option('exp_type').lower() == 'tpm' and self.option('method').lower() == 'rsem':
            self.start_listener()
            self.fire("start")
            self.dump_tpm(self.option('raw_task_id'),
                           self.option('transcript_main_id'),
                           self.option('gene_main_id'))
            # self.fire("end")
            # super(QuantWorkflow, self).end()
            self.end()
        else:
            self.tool.on("end", self.set_db)
            self.run_tool()
            super(QuantWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_exp", main_id=self.option('gene_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        all_exp = self.api.api("lnc_rna.all_exp")
        # add transcript info
        # exp_matrix = self.work_dir + '/transcript.tpm.matrix'
        exp_matrix = self.tool.output_dir + '/transcript.tpm.matrix'
        all_exp.add_exp(exp_matrix,
                        quant_method=self.option('method'),
                        exp_level='T',
                        main_id=self.option('transcript_main_id'),
                        exp_type='tpm',
                        add_distribution=False,
                        task_id=self.option("raw_task_id"),)
        # add gene info
        exp_matrix = self.tool.output_dir + '/gene.tpm.matrix'
        all_exp.add_exp(exp_matrix,
                        quant_method=self.option('method'),
                        exp_level='G',
                        main_id=self.option('gene_main_id'),
                        exp_type='tpm',
                        add_distribution=False,
                        task_id=self.option("raw_task_id"),)
        # paste annot
        self.paste_annotation()
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)

        self.inter_dirs = [
            ["01 Express", "", "表达量分析结果目录", 0],
            ["01 Express/01 Exp_Analysis", "", "表达量分析", 0]
        ]

        result_dir.add_relpath_rules([
            [".", "", "表达量分析文件", 0, "211097"],
            ["./gene.fpkm.matrix.annot.xls", "xls", "基因fpkm表达定量注释结果表", 0],
            ["./transcript.fpkm.matrix.annot.xls", "xls", "转录本fpkm表达定量注释结果表", 0],
            ["./gene.tpm.matrix.annot.xls", "xls", "基因tpm表达定量注释结果表", 0],
            ["./transcript.tpm.matrix.annot.xls", "xls", "转录本tpm表达定量注释结果表", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        super(QuantWorkflow, self).end()

    def run_tool(self):
        options = {
            "method": self.option('method'),
            "libtype": self.option('libtype'),
            "read_len": self.option('read_len'),
            "t2g": self.option('t2g').prop['path'],
            "transcriptome": self.option('transcriptome'),
            "fastq": self.option('fastq'),
        }
        self.tool.set_options(options)
        self.tool.run()

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        self.logger.info("path: {}".format(path))
        self.logger.info("path: {}".format(to_path))
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find')
        return to_path


    def dump_fpkm(self, task_id, transcript_main_id, gene_main_id):
        all_exp = self.api.api("lnc_rna.all_exp")
        collection = all_exp.db['sg_exp']
        results = collection.find({"task_id": task_id, "method": "RSEM"})
        transcript_fpkm_file = None
        gene_fpkm_file = None
        transcript_anno_file = None     # anno files added on 20201209
        gene_anno_file = None
        for each in results:
            if each['exp_level'] == 'T' and 'count_file' in each:
                transcript_fpkm_file = each['count_file'].replace('.count.matrix', '.fpkm.matrix')
                transcript_anno_file = each['count_file'].replace('.count.matrix', '.fpkm.matrix.annot')
                transcript_fpkm_file = self.download_s3_file(transcript_fpkm_file, os.path.basename(transcript_fpkm_file))
                transcript_anno_file = self.download_s3_file(transcript_anno_file,
                                                             os.path.join(self.tool.output_dir, os.path.basename(transcript_anno_file)))
            if each['exp_level'] == 'G' and 'count_file' in each:
                gene_fpkm_file = each['count_file'].replace('.count.matrix', '.fpkm.matrix')
                gene_anno_file = each['count_file'].replace('.count.matrix', '.fpkm.matrix.annot')
                gene_fpkm_file = self.download_s3_file(gene_fpkm_file, os.path.basename(gene_fpkm_file))
                gene_anno_file = self.download_s3_file(gene_anno_file,
                                                       os.path.join(self.tool.output_dir, os.path.basename(gene_anno_file)))
            if transcript_fpkm_file and gene_fpkm_file and transcript_anno_file and gene_anno_file:
                break
        if transcript_fpkm_file is not None:
            all_exp.add_exp(transcript_fpkm_file, quant_method="RSEM", exp_level='T',
                            main_id=transcript_main_id, exp_type='fpkm', add_distribution=False, task_id=self.option('raw_task_id'))
        else:
            self.set_error('transcript_fpkm/count_file is not found', code = "13702101")
        if gene_fpkm_file is not None:
            all_exp.add_exp(gene_fpkm_file, quant_method="RSEM", exp_level='G',
                            main_id=gene_main_id, exp_type='fpkm', add_distribution=False, task_id=self.option('raw_task_id'))
        else:
            self.set_error('gene_fpkm/count_file is not found', code = "13702102")

    def dump_tpm(self, task_id, transcript_main_id, gene_main_id):
        all_exp = self.api.api("lnc_rna.all_exp")
        collection = all_exp.db['sg_exp']
        results = collection.find({"task_id": task_id, "method": "RSEM"})
        transcript_fpkm_file = None
        gene_fpkm_file = None
        transcript_anno_file = None  # anno files added on 20201209
        gene_anno_file = None
        for each in results:
            if each['exp_level'] == 'T' and 'count_file' in each:
                transcript_fpkm_file = each['count_file'].replace('.count.matrix', '.tpm.matrix')
                transcript_anno_file = each['count_file'].replace('.count.matrix', '.tpm.matrix.annot')
                transcript_fpkm_file = self.download_s3_file(transcript_fpkm_file, os.path.basename(transcript_fpkm_file))
                transcript_anno_file = self.download_s3_file(transcript_anno_file,
                                                             os.path.join(self.tool.output_dir,os.path.basename(transcript_anno_file)))
            if each['exp_level'] == 'G' and 'count_file' in each:
                gene_fpkm_file = each['count_file'].replace('.count.matrix', '.tpm.matrix')
                gene_anno_file = each['count_file'].replace('.count.matrix', '.tpm.matrix.annot')
                gene_fpkm_file = self.download_s3_file(gene_fpkm_file, os.path.basename(gene_fpkm_file))
                gene_anno_file = self.download_s3_file(gene_anno_file,
                                                       os.path.join(self.tool.output_dir,os.path.basename(gene_anno_file)))
            if transcript_fpkm_file and gene_fpkm_file and transcript_anno_file and gene_anno_file:
                break
        if transcript_fpkm_file is not None:
            all_exp.add_exp(transcript_fpkm_file, quant_method="RSEM", exp_level='T',
                            main_id=transcript_main_id, exp_type='TPM', add_distribution=False, task_id=self.option('raw_task_id'))
        else:
            self.set_error('transcript_tpm/count_file is not found', code = "13702103")
        if gene_fpkm_file is not None:
            all_exp.add_exp(gene_fpkm_file, quant_method="RSEM", exp_level='G',
                            main_id=gene_main_id, exp_type='TPM', add_distribution=False, task_id=self.option('raw_task_id'))
        else:
            self.set_error('gene_tpm/count_file is not found', code = "13702104")

    def paste_annotation(self):
        all_exp = self.api.api("lnc_rna.all_exp")
        conn = all_exp.db['sg_annotation_stat']
        task_id = self.option('task_id')
        try:
            find_result = conn.find_one({"task_id": task_id, "type": "latest", "status": "end"})
        except:
            self.set_error('annotation result file cannot found', code = "13702105")
        else:
            if find_result is None:
                find_result = conn.find_one({"task_id": task_id, "type": "origin"})
        annot = find_result['result_dir'] + '/allannot_class/all_annot.xls'
        annot = self.download_s3_file(annot, "all_annot.xls")
        all_annot = pd.read_table(annot, header=0, index_col=0)
        # for gene
        annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(columns=['transcript_id', 'is_gene'])
        exp_pd = pd.read_table(self.tool.output_dir + '/gene.tpm.matrix', header=0, sep='\t', index_col=0)
        exp_pd = exp_pd.join(annot_pd, how='left')
        out_exp = self.tool.output_dir + '/gene.tpm.matrix' + '.annot.xls'
        exp_pd.to_csv(out_exp, sep='\t', header=True, index=True)
        # for transcript
        annot_pd = all_annot.reset_index().drop(columns=['is_gene']).set_index('transcript_id')
        exp_pd = pd.read_table(self.tool.output_dir + '/transcript.tpm.matrix', header=0, sep='\t', index_col=0)
        exp_pd = exp_pd.join(annot_pd, how='left')
        out_exp = self.tool.output_dir + '/transcript.tpm.matrix' + '.annot.xls'
        exp_pd.to_csv(out_exp, sep='\t', header=True, index=True)
