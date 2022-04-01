# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import pandas as pd
import glob
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.config import Config


class QuantWorkflow(Workflow):
    """
    quant
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(QuantWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name="biotype", type="infile", format="prok_rna.common"),
            dict(name="transcriptome", type="infile", format="prok_rna.fasta"),
            dict(name="fastq", type="infile", format="prok_rna.fastq_list"),
            dict(name="method", type="string", default="rsem"),
            dict(name="libtype", type="string", default=None),
            dict(name="read_len", type="int", default=149),
            dict(name="task_type", type="string"),
            dict(name="raw_task_id", type="string"),
            dict(name="submit_location", type="string"),
            dict(name="exp_type", type="string"),
            dict(name="id2name", type="infile", format="prok_rna.common"),  # for convert gene_id to gene_name
            # dict(name="exp_level", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        if options[3]:
            options.append(dict(name="main_id_ms", type="string", default=None))
            options.append(dict(name="main_id_m", type="string"))
            options.append(dict(name="main_id_s", type="string", default=None))
        else:
            options.append(dict(name="main_id_m", type="string"))
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_module("prok_rna.quant")
        self.extract_mrna = self.add_tool("prok_rna.extract_mrna")

    def run(self):
        db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
        col = db["sg_task"]
        result = col.find_one({"task_id": self.option('raw_task_id')})
        if "version" in result:
            self.version = col.find_one({"task_id": self.option('raw_task_id')})["version"]
        else:
            self.version = None
        if self.version and self.version in ["v1.1", "v3", "v3.1"]:
            result_dir = db["sg_exp"].find_one({"task_id": self.option('raw_task_id')})["result_dir"]
            if self.option('exp_type').lower() == 'fpkm' and self.option('method').lower() == 'rsem':
                self.start_listener()
                self.fire("start")
                self.dump_fpkm(self.option('raw_task_id'), result_dir)
                # self.fire("end")
                self.end()
            elif self.option('exp_type').lower() == 'tpm' and self.option('method').lower() == 'rsem':
                self.start_listener()
                self.fire("start")
                self.dump_tpm(self.option('raw_task_id'), result_dir)
                # self.fire("end")
                self.end()
        else:
            self.extract_mrna.on("end", self.run_tool)
            self.tool.on("end", self.set_db)
            self.run_extract_mrna()
            super(QuantWorkflow, self).run()

    def dump_fpkm(self, task_id, result_dir):
        transcript_fpkm_file = result_dir + "/transcript.fpkm.matrix"
        transcript_count_file = result_dir + "/transcript.count.matrix"
        transcript_count_file_path = result_dir + "/transcript.count.matrix"
        transcript_fpkm_anno_file = result_dir + "/transcript.fpkm.matrix.annot.xls"
        transcript_fpkm_file = self.download_s3_file(transcript_fpkm_file, os.path.basename(transcript_fpkm_file))
        transcript_count_file = self.download_s3_file(transcript_count_file, os.path.basename(transcript_count_file))
        transcript_fpkm_anno_file = self.download_s3_file(transcript_fpkm_anno_file, os.path.basename(transcript_fpkm_anno_file))
        all_exp = self.api.api("prok_rna.all_exp")
        all_exp.add_exp(transcript_fpkm_file, quant_method=self.option('method'), main_id_m=self.option('main_id_m'),
                                    main_id_ms=self.option('main_id_ms'), main_id_s=self.option('main_id_s'),
                                    add_distribution=False, exp_type='fpkm', lib_type=self.option('libtype'),
                                    task_id=self.option('raw_task_id'), count_file_path=transcript_count_file_path)
        ## set output
        if transcript_fpkm_anno_file is not None:
            os.link(transcript_fpkm_anno_file, os.path.join(self.output_dir, os.path.basename(transcript_fpkm_anno_file)))
        if transcript_fpkm_file is not None:
            os.link(transcript_fpkm_file, os.path.join(self.output_dir, os.path.basename(transcript_fpkm_file)) + ".xls")
        if transcript_count_file is not None:
            os.link(transcript_count_file, os.path.join(self.output_dir, os.path.basename(transcript_count_file)) + ".xls")

    def dump_tpm(self, task_id, result_dir):
        transcript_tpm_file = result_dir + "/transcript.tpm.matrix"
        transcript_count_file = result_dir + "/transcript.count.matrix"
        transcript_tpm_anno_file = result_dir + "/transcript.tpm.matrix.annot.xls"
        transcript_tpm_file = self.download_s3_file(transcript_tpm_file, os.path.basename(transcript_tpm_file))
        transcript_count_file = self.download_s3_file(transcript_count_file, os.path.basename(transcript_count_file))
        transcript_tpm_anno_file = self.download_s3_file(transcript_tpm_anno_file, os.path.basename(transcript_tpm_anno_file))
        all_exp = self.api.api("prok_rna.all_exp")
        all_exp.add_exp(transcript_tpm_file, quant_method=self.option('method'), main_id_m=self.option('main_id_m'),
                                    main_id_ms=self.option('main_id_ms'), main_id_s=self.option('main_id_s'),
                                    add_distribution=False, exp_type='tpm', lib_type=self.option('libtype'),
                                    task_id=self.option('raw_task_id'), )
        ## set output
        if transcript_tpm_anno_file is not None:
            os.link(transcript_tpm_anno_file, os.path.join(self.output_dir, os.path.basename(transcript_tpm_anno_file)))
        if transcript_tpm_file is not None:
            os.link(transcript_tpm_file, os.path.join(self.output_dir, os.path.basename(transcript_tpm_file)) + ".xls")
        if transcript_count_file is not None:
            os.link(transcript_count_file, os.path.join(self.output_dir, os.path.basename(transcript_count_file)) + ".xls")

    def run_extract_mrna(self):
        self.logger.info("开始提取mRNA+sRNA序列")
        opts = {
            "fasta" : self.option("transcriptome").prop['path'],
            "biotype" : self.option("biotype"),
            "rna_type" : "mRNA+sRNA",
        }
        self.extract_mrna.set_options(opts)
        self.extract_mrna.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        self.get_anno_file()
        all_exp = self.api.api("prok_rna.all_exp")
        if self.option('method').lower() == 'rsem':
            exp_matrix_tpm = self.tool.output_dir + '/transcript.tpm.matrix'
            exp_matrix_fpkm = self.tool.output_dir + '/transcript.fpkm.matrix'
            if self.option('libtype'):
                if self.option("exp_type").lower() == "tpm":
                    all_exp.add_exp(exp_matrix_tpm, quant_method=self.option('method'), main_id_m=self.option('main_id_m'),
                                    main_id_ms=self.option('main_id_ms'), main_id_s=self.option('main_id_s'),
                                    add_distribution=False, exp_type='tpm', lib_type=self.option('libtype'),
                                    task_id=self.option('raw_task_id'), )
                else:
                    all_exp.add_exp(exp_matrix_fpkm, quant_method=self.option('method'), main_id_m=self.option('main_id_m'),
                                    main_id_ms=self.option('main_id_ms'), main_id_s=self.option('main_id_s'),
                                    add_distribution=False, exp_type='fpkm', lib_type=self.option('libtype'),
                                    task_id=self.option('raw_task_id'), )
                self.paste_annotation(exp_matrix_tpm)
                self.paste_annotation(exp_matrix_fpkm)
                self.end()
            else:
                if self.option("exp_type").lower() == "tpm":
                    all_exp.add_exp(exp_matrix_tpm, quant_method=self.option('method'), main_id_m=self.option('main_id_m'),
                                add_distribution=False, exp_type='tpm', task_id=self.option('raw_task_id'))
                else:
                    all_exp.add_exp(exp_matrix_fpkm, quant_method=self.option('method'), main_id_m=self.option('main_id_m'),
                                add_distribution=False, exp_type='fpkm', task_id=self.option('raw_task_id'))
                self.paste_annotation(exp_matrix_tpm)
                self.paste_annotation(exp_matrix_fpkm)
                self.end()
        else:
            exp_matrix = self.tool.output_dir + '/transcript.tpm.matrix'
            if self.option('libtype'):
                all_exp.add_exp(exp_matrix, quant_method=self.option('method'), main_id_m=self.option('main_id_m'),
                                main_id_ms=self.option('main_id_ms'), main_id_s=self.option('main_id_s'),
                                add_distribution=False, exp_type='tpm', lib_type=self.option('libtype'),
                                task_id=self.option('raw_task_id'), )
            else:
                all_exp.add_exp(exp_matrix, quant_method=self.option('method'), main_id_m=self.option('main_id_m'),
                                add_distribution=False, exp_type='tpm', task_id=self.option('raw_task_id'))
            self.paste_annotation(exp_matrix)
            self.end()

    def end(self):
        if self.version in ["v1.1", "v3", "v3.1"]:
            result_dir = self.add_upload_dir(self.output_dir)
        else:
            result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_regexp_rules([
            [".", "", "表达定量分析结果目录"],
            [r"transcript.count.matrix.xls", "", "transcript count表达定量结果表"],
            [r"transcript.tpm.matrix.xls", "", "transcript tpm表达定量结果表"],
            [r"transcript.tpm.matrix.annot.xls", "", "transcript tpm表达定量注释结果表"],
            [r"transcript.fpkm.matrix.xls", "", "transcript fpkm表达定量结果表"],
            [r"transcript.fpkm.matrix.annot.xls", "", "transcript fpkm表达定量注释结果表"],
        ])
        super(QuantWorkflow, self).end()

    def run_tool(self):
        options = {
            "method": self.option('method'),
            "libtype": self.option('libtype'),
            "read_len": self.option('read_len'),
            "transcriptome": self.extract_mrna.option("out_fasta"),
            "fastq": self.option('fastq'),
            "id2name": self.option('id2name').prop['path'],
            "biotype": self.option("biotype"),
            "task_id": self.option("raw_task_id")
        }
        self.tool.set_options(options)
        self.tool.run()

    def paste_annotation(self, matrix):
        annot = os.path.join(self.work_dir, "annot.xls")
        all_annot = pd.read_table(annot, header=0, index_col=0)
        all_annot.rename(columns={'gene_id': 'seq_id'}, inplace=True)
        exp_pd = pd.read_table(matrix, header=0, sep='\t', index_col=0)
        exp_pd.drop(['gene_name'], axis=1, inplace=True)
        exp_pd = exp_pd.join(all_annot, how='left')
        out_exp = matrix + '.annot.xls'
        exp_pd.to_csv(out_exp, sep='\t', header=True, index=True)

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
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path

    def get_anno_file(self):
        all_exp = self.api.api("prok_rna.all_exp")
        conn = all_exp.db['sg_annotation_stat']
        task_id = self.option('raw_task_id')
        find_result = conn.find_one({"task_id": task_id, "type": "origin", "status": "end"})
        annot = os.path.join(find_result['result_dir'], 'summary/all_anno_detail.xls')
        if os.path.exists(annot):
            os.system('cp {} {}'.format(annot, os.path.join(self.work_dir, "annot.xls")))
        else:
            annot = self.download_s3_file(annot, "annot.xls")
