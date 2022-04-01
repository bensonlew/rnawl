# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190128

import os
import time
import json
import shutil
import gevent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from upload_s5cmd import UploadS5cmd


class SampleQcWorkflow(Workflow):
    """
    样本质控，用于动植物基因组、微生物基因组、宏基因组、常规RNA、miRNA、lncRNA、原核转录组
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SampleQcWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "project_params", "type": "infile", "format": "datasplit.library_params"},  # 参数文件
            {"name": "run_type", "type": "string", "default": "auto"},  # 是否进行自动拆分,自动进行下面的样本拆分
            {"name": "split_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "is_upload_cleandata","type":"string"}, #增加释放cleandata的部分 from Qinwen20211019
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.start_times, self.end_times = 0, 0
        self.all_modules = []
        # self.md5sum = self.add_tool("datasplit_v2.md5sum")

    def check_options(self):
        if not self.option("project_params").is_set:
            raise OptionError("请设置参数文件")

    def get_params(self):
        """
        解析输入json文件，得到其中的参数信息
        """
        f = open(self.option("project_params").prop["path"], "rb")
        try:
            self.json_dict = json.loads(f.read())
        except:
            raise OptionError("json格式不正确")
        return self.json_dict

    def run_fastp_qc(self):
        """
        将属于meta_genomic(宏基因组)/rna(常规RNA)/dna(动植物基因组)项目
        的样本都用fastp进行质控
        """
        for project_type in self.json_dict.keys():
            if project_type in ["meta_genomic", "dna", "rna"]:
                for i in range(len(self.json_dict[project_type])):
                    qc_params = self.json_dict[project_type][i]
                    types = qc_params.keys()
                    options = {"sample_path": qc_params["sample_path"]}
                    if "qualified_quality_phred" in types:
                        options["qualified_quality_phred"] = qc_params["qualified_quality_phred"]
                    if "length_required" in types:
                        options["length_required"] = qc_params["length_required"]
                    if "cut_by_quality5" in types:
                        options["cut_by_quality5"] = qc_params["cut_by_quality5"]
                    if "cut_by_quality3" in types:
                        options["cut_by_quality3"] = qc_params["cut_by_quality3"]
                    if "cut_mean_quality" in types:
                        options["cut_mean_quality"] = qc_params["cut_mean_quality"]
                    if "n_base_limit" in types:
                        options["n_base_limit"] = qc_params["n_base_limit"]
                    if "cut_window_size" in types:
                        options["cut_window_size"] = qc_params["cut_window_size"]
                    if "adapter_sequence" in types:
                        options["adapter_sequence"] = qc_params["adapter_sequence"]
                    if "adapter_sequence_r2" in types:
                        options["adapter_sequence_r2"] = qc_params["adapter_sequence_r2"]
                    self.fastp_qc = self.add_module("datasplit_v2.fastp")
                    # self.start_times += 1
                    self.fastp_qc.set_options(options)
                    self.fastp_qc.on("end", self.set_output, project_type)
                    # self.fastp_qc.run()
                    self.all_modules.append(self.fastp_qc)
                    self.start_times += 1

    def run_rna_rfam(self):
        """
        prokaryotic_rna(原核mRNA)/lncrna(LncRNA)项目质控加rfam统计核糖体占比
        """
        for project_type in self.json_dict.keys():
            if project_type in ["prokaryotic_rna", "lncrna"]:
                for i in range(len(self.json_dict[project_type])):
                    qc_params = self.json_dict[project_type][i]
                    types = qc_params.keys()
                    options = {"sample_path": qc_params["sample_path"], "project_type": project_type}
                    if "qualified_quality_phred" in types:
                        options["qualified_quality_phred"] = qc_params["qualified_quality_phred"]
                    if "length_required" in types:
                        options["length_required"] = qc_params["length_required"]
                    if "cut_by_quality5" in types:
                        options["cut_by_quality5"] = qc_params["cut_by_quality5"]
                    if "cut_by_quality3" in types:
                        options["cut_by_quality3"] = qc_params["cut_by_quality3"]
                    if "cut_mean_quality" in types:
                        options["cut_mean_quality"] = qc_params["cut_mean_quality"]
                    if "n_base_limit" in types:
                        options["n_base_limit"] = qc_params["n_base_limit"]
                    if "adapter_sequence" in types:
                        options["adapter_sequence"] = qc_params["adapter_sequence"]
                    if "adapter_sequence_r2" in types:
                        options["adapter_sequence_r2"] = qc_params["adapter_sequence_r2"]
                    self.rna_qc_rfam = self.add_module("datasplit_v2.rna_qc_rfam")
                    # self.start_times += 1
                    self.rna_qc_rfam.set_options(options)
                    self.rna_qc_rfam.on("end", self.set_output, project_type)
                    # self.rna_qc_rfam.run()
                    self.all_modules.append(self.rna_qc_rfam)
                    self.start_times += 1

    def run_microbial_genome_qc(self):
        """
        microbial_genome(微生物基因组质控)
        """
        microbialgenome_info = self.json_dict["microbial_genome"]
        for i in range(len(microbialgenome_info)):
            microbialgenome_step = microbialgenome_info[i]
            options = {"sample_info": microbialgenome_step["sample_path"]}
            types = microbialgenome_step.keys()
            if "readl" in types:
                options["readl"] = microbialgenome_step["readl"]
            if "qualified_quality_phred" in types:
                options["qualified_quality_phred"] = microbialgenome_step["qualified_quality_phred"]
            if "length_required" in types:
                options["length_required"] = microbialgenome_step["length_required"]
            if "cut_by_quality5" in types:
                options["cut_by_quality5"] = microbialgenome_step["cut_by_quality5"]
            if "cut_by_quality3" in types:
                options["cut_by_quality3"] = microbialgenome_step["cut_by_quality3"]
            if "cut_mean_quality" in types:
                options["cut_mean_quality"] = microbialgenome_step["cut_mean_quality"]
            if "n_base_limit" in types:
                options["n_base_limit"] = microbialgenome_step["n_base_limit"]
            if "adapter_sequence" in types:
                options["adapter_sequence"] = microbialgenome_step["adapter_sequence"]
            if "adapter_sequence_r2" in types:
                options["adapter_sequence_r2"] = microbialgenome_step["adapter_sequence_r2"]
            self.microbialgenome_fastp = self.add_module("datasplit_v2.microbial_genome_qc")
            # self.start_times += 1
            self.microbialgenome_fastp.set_options(options)
            self.microbialgenome_fastp.on("end", self.set_output, "microbial_genome")
            # self.microbialgenome_fastp.run()
            self.all_modules.append(self.microbialgenome_fastp)
            self.start_times += 1

    def run_mirna_qc(self):
        """
        miRNA质控
        """
        mirna_info = self.json_dict["mirna"]
        for i in range(len(mirna_info)):
            mirna_step = mirna_info[i]
            options = {"list_file": mirna_step["sample_path"]}
            types = mirna_step.keys()
            if "length" in types:
                options["length"] = mirna_step["length"]
            if "adapter" in types:
                options["adapter"] = mirna_step["adapter"]
            if "phred_score" in types:
                options["phred_score"] = mirna_step["phred_score"]
            if "minlen" in types:
                options["minlen"] = mirna_step["minlen"]
            if "max_length" in types:
                options["max_length"] = mirna_step["max_length"]
            if "cut_left" in types:
                options["cut_left"] = mirna_step["cut_left"]
            self.mirna = self.add_module("datasplit_v2.mirna_qc")
            # self.start_times += 1
            self.mirna.set_options(options)
            self.mirna.on("end", self.set_output, "mirna")
            # self.mirna.run()
            self.all_modules.append(self.mirna)
            self.start_times += 1

    def run_md5sum(self):
        self.logger.info("开始进行md5校验")
        self.md_tool = []
        if self.option("is_upload_cleandata"):#增加了新的拆分参数判定用于释放cleandata from Qinwen20211019
            for f in os.listdir(self.output_dir):
                if f.split("-")[0] in ["meta_genomic", "microbial_genome", "rna", "dna", "prokaryotic_rna", "lncrna","mirna"]:
                    options = {"fastq_dir": os.path.join(self.output_dir, f + "/fastq")}
                self.md5sum = self.add_tool("datasplit_v2.md5sum")
                self.md5sum.set_options(options)
                self.md_tool.append(self.md5sum)
        else:
            for f in os.listdir(self.output_dir):
                if f in ["meta_genomic_clean", "microbial_genome_clean", "rna_clean", "dna_clean", "prokaryotic_rna_clean", "lncrna_clean"]:
                    options = {"fastq_dir": os.path.join(self.output_dir, f + "/fastq")}
                if f == "mirna_clean":
                    options = {"fastq_dir": os.path.join(self.output_dir, f + "/fasta")}
                self.md5sum = self.add_tool("datasplit_v2.md5sum")
                self.md5sum.set_options(options)
                self.md_tool.append(self.md5sum)
                if f == "mirna_clean":
                    options = {"fastq_dir": os.path.join(self.output_dir, f + "/fastq")}#mirna的fq.gz上传md5值
                    self.md5sum = self.add_tool("datasplit_v2.md5sum")
                    self.md5sum.set_options(options)
                    self.md_tool.append(self.md5sum)
        if len(self.md_tool) > 1:
            # self.on_rely(self.md_tool, self.set_db)
            self.on_rely(self.md_tool, self.set_db)
            for md_tool in self.md_tool:
                md_tool.run()
        elif len(self.md_tool) == 1:
            # self.md_tool[0].on("end", self.set_db)
            self.md_tool[0].on("end", self.set_db)
            self.md_tool[0].run()
        else:
            # self.set_db()
            self.end()

    def link_dir(self, old_dir, new_dir):
        if os.path.exists(new_dir):
            shutil.rmtree(new_dir)
        os.mkdir(new_dir)
        for f in os.listdir(old_dir):
            old = os.path.join(old_dir, f)
            new = os.path.join(new_dir, f)
            if os.path.isdir(old):
                self.link_dir(old, new)
            else:
                if os.path.exists(new):
                    raise Exception("{}文件重复，请检查".format(new))
                os.link(old, new)

    def set_output(self, event):
        obj = event["bind_object"]
        if self.option("is_upload_cleandata"):#释放cleandata可能会有相同类型的产品不同的参数，加了endtime去分别from Qinwen20211019
            self.end_times += 1
            self.link_dir(obj.output_dir, os.path.join(self.output_dir, event["data"]+"-"+str(self.end_times) + "_clean"))
        else:
            self.link_dir(obj.output_dir, os.path.join(self.output_dir, event["data"] + "_clean"))
            self.end_times += 1
        if self.start_times == self.end_times:
            time.sleep(40)# 将timesleep转移到这里
            # self.set_db()
            self.run_md5sum()
        #     self.end()

    def set_db(self):
        # time.sleep(40)
        self.logger.info("开始进行样本信息导表")
        if self.option("is_upload_cleandata") and self.option("split_id"):#释放cleandata的导表，导入到sg_qc from Qinwen20211019
            datasplit_api = self.api.api("datasplit.datasplit_qc")
            qc_id = self.option("split_id")
            for pro in os.listdir(self.output_dir):
                fq_dir = os.path.join(self.output_dir, pro+ "/" + "fastq")
                fq_s3_upload_dir = os.path.join(self._sheet.output, pro + "/" + "fastq")
                datasplit_api.update_sg_qc_specimen(qc_id,fq_dir,fq_s3_upload_dir)
        elif self.option("split_id"):
            datasplit_api = self.api.api("datasplit.datasplit_new")
            split_id = self.option("split_id")
            for pro in os.listdir(self.output_dir):
                product_type = pro.split("_clean")[0]
                pro_dir = os.path.join(self.output_dir, pro)
                if pro == "mirna_clean":
                    fasta_dir = os.path.join(pro_dir, "fasta")
                    fastq_dir = os.path.join(pro_dir, "fastq")
                    # s3_upload_dir = os.path.join(self._sheet.output, pro + "/" + "fasta")
                    fa_s3_upload_dir = os.path.join(self._sheet.output)
                    # fq_s3_upload_dir = os.path.join(self._sheet.output, pro + "/" + "fastq")
                    fq_s3_upload_dir = os.path.join(self._sheet.output)
                    fastq_stat = os.path.join(pro_dir, "fastq_stat.xls")
                    sample_qc_stat = os.path.join(pro_dir, "Sample_QC_stat.xls")
                    fasta_length_dir = os.path.join(pro_dir, "fasta_length")
                    datasplit_api.update_sample_path(split_id, fasta_dir, fa_s3_upload_dir, product_type)
                    datasplit_api.update_sample_path(split_id, fastq_dir, fq_s3_upload_dir, product_type)
                    datasplit_api.add_sg_split_qc_mirna(split_id, fastq_stat, sample_qc_stat)
                    datasplit_api.add_mirna_sg_bar(split_id, fasta_length_dir)
                else:
                    for f in os.listdir(pro_dir):
                        if f == "json_dir":
                            json_dir = os.path.join(pro_dir, f)
                            datasplit_api.add_fastp_json_stat(split_id, json_dir, product_type)
                        # elif f == "rfam_dir":
                        #     rfam_dir = os.path.join(pro_dir, f)
                        #     datasplit_api.add_rfam_stat(split_id, rfam_dir, product_type)
                        elif f == "fastq":
                            fastq_dir = os.path.join(pro_dir, f)
                            s3_upload_dir = os.path.join(self._sheet.output, pro + "/" + f)
                            datasplit_api.update_sample_path(split_id, fastq_dir, s3_upload_dir, product_type)
            datasplit_api.delete_sg_split_clean_qc(split_id)
            datasplit_api.delete_sg_split_library_qc(split_id)
        self.logger.info("质控导表完成")
        # self.run_md5sum()
        self.end()

    def run(self):
        self.get_params()
        kk = self.json_dict.keys()
        if kk == ["meta"] or not kk:
            self.start_listener()
            self.fire("start")
            gevent.spawn_later(5, self.end)
        self.run_fastp_qc()
        self.run_rna_rfam()
        if "microbial_genome" in self.json_dict.keys():
            self.run_microbial_genome_qc()
        if "mirna" in self.json_dict.keys():
            self.run_mirna_qc()
        self.logger.info(len(self.all_modules))
        if len(self.all_modules) > 1:
            # self.on_rely(self.all_modules, self.run_md5sum)
            for module in self.all_modules:
                module.run()
        elif len(self.all_modules) == 1:
            # self.all_modules[0].on("end", self.run_md5sum)
            self.all_modules[0].run()
        super(SampleQcWorkflow, self).run()

    def end(self):
        # self.add_upload_dir(self.output_dir)
        # u = UploadS5cmd()
        # u.upload(self.output_dir, self._sheet.output)
        if self.option("is_upload_cleandata"):#释放clean_data需要上传全部cleandata
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".","","结果输出目录"]
            ])
            result_dir.add_regexp_rules([
                ["","",""]
            ])
        else:
            mirna_clean = os.path.join(self.output_dir, "mirna_clean/fastq")  # 修改时注意对应的s3路径的更新
            if os.path.exists(mirna_clean):
                self.add_upload_dir(mirna_clean)
            mirna_clean = os.path.join(self.output_dir, "mirna_clean/fasta")
            if os.path.exists(mirna_clean):
                self.add_upload_dir(mirna_clean)
        super(SampleQcWorkflow, self).end()
