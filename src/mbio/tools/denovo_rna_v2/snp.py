# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
import re
from biocluster.core.exceptions import OptionError
import pandas as pd
import glob

class SnpAgent(Agent):
    """
    对无参转录组序列进行snp的预测。
    检查tool的详细的错误可以去tool的结果里面一个err的文件查看
    """
    def __init__(self, parent):
        super(SnpAgent, self).__init__(parent)
        options = [
            {"name": "trinity_fa", "type": "infile", "format":"denovo_rna_v2.trinity_fasta"},  # 这是fa文件的格式检查
            {"name": "bamlist", "type": "infile", "format": "denovo_rna_v2.bamlist"},
            {"name": "call_vcf", "type": "outfile", "format": "denovo_rna_v2.vcf"},
            {"name": "bed_file", "type": "string"},  # 切割后的转录本文件
            #{"name": "isoform_unigene", "type": "string"},#实际为一个文件，只是不检查

        ]
        self.add_option(options)
        self.step.add_steps("snp_predict")
        self.on("start", self.step_start)
        self.on("end", self.step_end)


    def step_start(self):
        self.step.snp_predict.start()
        self.step.update()

    def step_end(self):
        self.step.snp_predict.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数

        :return:
        """
        if not self.option("trinity_fa").is_set:
            raise OptionError("必须设置输入文件:组装完成的unigene.fa文件", code = "32005501")

        if not self.option("bamlist"):
            raise OptionError("必须设置输入文件:rsem比对得到的bamlist文件", code = "32005502")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_regexp_rules([
            [".", "", "结果输出目录"],
            [r"call.vcf$", " ", "bcftools call得到的vcf长度分布信息统计文件"],
        ])
        super(SnpAgent, self).end()


class SnpTool(Tool):
    def __init__(self, config):
        super(SnpTool, self).__init__(config)
        # self.samtools_path = "/bioinfo/align/samtools-1.6/samtools-1.6/samtools"
        #self.bcftools_path = "/bioinfo/seq/bcftools-1.4/bin/bcftools"
        # self.bcftools_path = "/bioinfo/align/samtools-1.6/bcftools-1.6/bcftools"
        self.samtools_path = "bioinfo/rna/miniconda2/bin/samtools"
        # self.bcftools_path = "/bioinfo/seq/bcftools-1.4/bin/bcftools"
        self.bcftools_path = "bioinfo/rna/miniconda2/bin/bcftools"

    def mpileup(self):
        self.logger.info(self.option("trinity_fa").prop["path"])
        self.logger.info(vars(self.option("bamlist")))
        cmd2 = "{} mpileup -t DP,AD -u -f {} -q 10 -Q 20 -C 50 -b {} -go {}.bcf " \
               "-l {}" \
        .format(self.samtools_path, self.option("trinity_fa").prop["path"], self.option("bamlist").prop["path"], os.path.join(self.work_dir, os.path.basename(self.option('bed_file'))), self.option("bed_file"))
        # cmd2 = "{} mpileup -t DP,AD -u -f {} -q 10 -Q 20 -C 50 -b {} -go {}".format(self.samtools_path, self.option("trinity_fa").prop["path"], self.option("bamlist").prop["path"], self.work_dir + "/bcf")
        self.logger.info("运行mpileup")
        self.logger.info(cmd2)
        cmd2_obj = self.add_command("cmd2", cmd2).run()
        self.wait(cmd2_obj)
        if cmd2_obj.return_code == 0:
            self.logger.info("运行mpileup完成")
        else:
            cmd2_obj.rerun()
            self.wait(cmd2_obj)
            if cmd2_obj.return_code == 0:
                self.logger.info("运行mpileup完成")
            else:
                self.set_error("运行mpileup出错", code = "32005503")

    def download_s3_file(self, path, to_path):
        """
        对象存储下载单个文件：
        返回结果路径
        path: s3 文件路径
        to_path: 下载的s3 文件目的路径 如果前面包含路径会自动转为 base_name
        """
        if os.path.exists(to_path):
            pass
        elif os.path.exists(path):
            to_path = path
        elif path.startswith("s3"):
            self.download_from_s3(path, to_path)
            m = re.match(r"^([\w\-]+)://([\w\-]+)/(.*)$", path)
            if not m:
                self.set_error("下载路径%s格式不正确!" , variables=( path), code="32005505")
            to_path = self.work_dir +  "/" + os.path.basename(to_path)
        else:
            path = "s3://" + path
            self.download_from_s3(path, to_path)
            m = re.match(r"^([\w\-]+)://([\w\-]+)/(.*)$", path)
            if not m:
                self.set_error("下载路径%s格式不正确!" , variables=( path), code="32005506")
            to_path = self.work_dir +  "/" + os.path.basename(to_path)
        return to_path

    def download_bam(self):
        with open(self.option("bamlist").prop['path'], 'r') as bams_in, open("bamlist_download", 'w') as bams_out:
            for bam in bams_in.readlines():
                bam_new = self.download_s3_file(bam.strip(), os.path.basename(bam.strip()))
                bams_out.write("{}\n".format(bam_new))

        self.option("bamlist", os.path.join(self.work_dir, "bamlist_download"))

    def bcfcall(self):
        # os.mkdir(os.path.join(self.output_dir, os.path.basename(self.option("bed_file"))))
        cmd3 = "{} call -mv --format-fields GQ,GP --output-type v {} -o {}.call.vcf".format(self.bcftools_path, os.path.join(self.work_dir, os.path.basename(self.option('bed_file'))) + ".bcf", os.path.join(self.output_dir, os.path.basename(self.option("bed_file"))))
        self.logger.info("运行bcfcall")
        self.logger.info(cmd3)
        cmd3_obj = self.add_command("cmd3", cmd3).run()
        self.wait(cmd3_obj)
        if cmd3_obj.return_code == 0:
            self.logger.info("运行bcfcall完成")
        else:
            cmd3_obj.rerun()
            self.wait(cmd3_obj)
            if cmd3_obj.return_code == 0:
                self.logger.info("运行bcfcall完成")
            else:
                self.set_error("运行bcfcall出错", code = "32005504")
                self.set_error("运行bcfcall出错", code="32005507")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("把所有文件移动到这个tool的self.output_dir")

    def t2u(self, call_vcf, isoform_unigene):
        with open(call_vcf, "r") as filter, open(self.work_dir + "/new_call.vcf", "w") as ncf:
            for line in filter:
                # if len(line.split("#")) > 2:
                if line.startswith("##"):
                    continue
                else:
                    ncf.write(line)
        new_snp_rewrite = pd.read_table(self.work_dir + '/new_call.vcf', header=0, sep = '\t')
        new_snp_rewrite.rename(columns={new_snp_rewrite.columns[0]: "CHROM"},inplace=True)
        df_2 = pd.read_table(isoform_unigene, header=None, sep = '\t')
        df_2.columns = ["CHROM", "unigene", "yes_no"]
        merge_df12 = pd.merge(df_2, new_snp_rewrite, on="CHROM")
        merge_df12_new = merge_df12[merge_df12.yes_no.isin(['yes'])]
        merge_df12_new.drop(['CHROM', 'yes_no'], axis=1,inplace=True)
        merge_df12_new.rename(columns={merge_df12_new.columns[0]: "CHROM"},inplace=True)
        merge_df12_new.to_csv(self.work_dir + "/call_vcf_input",sep = "\t", index=False)
        merge_df12_new.drop(['FILTER'], axis=1).to_csv(self.work_dir + "/call_vcf_new",sep = "\t", index=False)# 这个是作为snpfinal这个tool的输入文件，没有删除filter

    def t2u_new(self, Trinity_fasta_t2g2u, call_vcf):
        def _add111(matched):
            intStr = matched.group("number")
            intValue = int(intStr)
            addedValue = intValue + 111
            addedValueStr = str(addedValue)
            return addedValueStr
        with open(Trinity_fasta_t2g2u, "r") as t2u:
            t2g_dcit = dict()
            for line_new in t2u:
                line_new = line_new.strip().split("\t")
                if line_new[2] == "yes":
                    t2g_dcit[line_new[0]] = line_new[1]
        # os.mkdir(open(os.path.join(self.output_dir, os.path.basename(self.option("bed_file")))))
        with open(call_vcf, "r") as call, open(os.path.join(self.output_dir,os.path.basename(self.option("bed_file"))) + '.origin', "w") as f_or:
            pattern = re.compile(r'>*ID=(.*),>*')
            for line in call:
                if not line.startswith("##contig=") and "##" in line:
                    f_or.write(line)
                elif line.startswith("##contig=") and pattern.search(line).group(1) in t2g_dcit.keys():
                    regex = pattern.search(line).group(1)
                    replacedStr = re.sub(regex, t2g_dcit[regex], line)
                    f_or.write(replacedStr + "\n")
                elif "#CHROM	POS	ID" in line:
                    f_or.write(line)
                elif "#" not in line:
                    line = line.strip().split("\t")
                    if line[0] in t2g_dcit.keys():
                        line[0] = t2g_dcit[line[0]]
                        seq = line[1:]
                        f_or.write(line[0] + "\t" + "\t".join(seq) + "\n")
                else:
                    pass

    def run(self):
        super(SnpTool, self).run()
        self.download_bam()
        self.mpileup()
        self.bcfcall()
        # self.t2u_new(self.option("isoform_unigene"), os.path.join(self.output_dir, os.path.basename(self.option("bed_file"))) + ".call.vcf")
        self.set_output()
        self.end()
