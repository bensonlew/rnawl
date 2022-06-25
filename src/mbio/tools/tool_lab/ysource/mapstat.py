# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import datetime
import random
import re
# import xlrd
import os
import sys
import shutil

class MapstatAgent(Agent):
    def __init__(self, parent):
        super(MapstatAgent, self).__init__(parent)
        options = [
            {"name":"bam_file","type":"infile","format":"align.bwa.bam"},
            {"name":"bam_prermdup","type":"infile","format":"align.bwa.bam"},
            {"name":"vcf_yfull","type":"infile","format":"sequence.vcf"},
            {"name":"clean_statistic","type":"string"},
            {"name":"sample_name","type":"string"},
            # {"name":"posdb_bed","type":"string"},
            # {"name":"EX_autosnp_bed","type":"string"},
            # {"name":"chrY_bed","type":"string"},
            # {"name":"ref","type":"string"}
        ]
        self.add_option(options)

    def check_option(self):
        '''
        参数检查
        '''
        if not self.option("bam_file"):
            raise OptionError("没有找到bam文件，请检查上一步")

    def set_resource(self):
        self._cpu = 2 
        self._memory = '20G'

    def end(self):
        super(MapstatAgent,self).end()

class MapstatTool(Tool):
    def __init__(self, config):
        super(MapstatTool, self).__init__(config)
        self._version = "v1.0"
        self.software_dir = self.config.SOFTWARE_DIR
        self.samtools = os.path.join(self.software_dir, 'miniconda2/bin/samtools')
        self.bam = self.option('bam_file').prop["path"]
        self.gatk_path = self.software_dir + "/bioinfo/gene-structure/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"
        self.java = os.path.join(self.software_dir, "/program/sun_jdk1.8.0/bin/java")
        self.sample_name = self.option("sample_name")
        # self.posdb_bed = self.option('posdb_bed')
        # self.EX_autosnp_bed = self.option('EX_autosnp_bed')
        # self.chrY_bed = self.option('chrY_bed')
        # self.ref = self.option('ref')
        self.posdb_bed = os.path.join(self.software_dir,"database/Tool_lab/ysource/posdb_update_20200801.bed")
        self.ref = os.path.join(self.software_dir, 'database/Tool_lab/ysource/reference/hg38modMt.fa')
        self.EX_autosnp_bed = os.path.join(self.software_dir,"database/Tool_lab/ysource/EX_autosnp2_2.bed")
        self.chrY_bed = os.path.join(self.software_dir,"database/Tool_lab/ysource/Y_chr_hg38.bed")
        self.bam_prermdup = self.option("bam_prermdup").prop["path"]

    def run(self):
        super(MapstatTool, self).run()
        self.run_gatk_doc("exon_result",self.posdb_bed)
        self.run_gatk_doc("EX_target_result",self.EX_autosnp_bed)
        self.run_gatk_doc("YFULL_target_result",self.chrY_bed)
        self.run_samstat(self.bam_prermdup,"bam_statistic_pre_dedup")
        self.run_samstat(self.bam,"bam_statistic")
        self.run_samdepth()
        self.run_sample_MAPOT()
        self.run_sample_dep()
        self.set_output()
        self.end()

        
    def run_gatk_doc(self,o_name,bed):
        '''
        GATK 3.8
        '''
        cmd = "{} -Xmx6g -jar {} -T DepthOfCoverage -R {} ".format(self.java,self.gatk_path,self.ref)
        in_bam = self.option("bam_file").prop["path"]
        cmd += "-o {} -I {} -L {} ".format(o_name, in_bam, bed)
        cmd += "--omitDepthOutputAtEachBase --omitIntervalStatistics -ct 1 -ct 5 -ct 10 -ct 20"
        command = self.add_command("gatk_depthofcoverage_{}".format(o_name.lower()), cmd, ignore_error = True)
        command.run()
        self.wait(command)
        n = 0
        if command.return_code == 0:
            self.logger.info("运行GATK DepthOfCoverage完成")
        # elif n >1:
        #     self.logger.info("重复运行GATKDOC失败,跳过此步骤")
        else:
            while 1:
                if n > 5:
                    self.logger.info("重复运行GATKDOC失败,跳过此步骤")
                    break
                n += 1
                self.logger.warn("运行GATKDOC失败,尝试重新运行")
                command.rerun()
                self.wait()
                if command.return_code is 0:
                    self.logger.info("运行GATK DepthOfCoverage完成")
                    break


    def run_samstat(self,bam_file, out_name):
        # bam_file = self.option("bam_file")
        # output = self.option("output_name")
        cmd = "{} stats {} > {}".format(self.samtools,bam_file,out_name)
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.software_dir + '/bioinfo/tool_lab/Yoogene/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "samstat_{}.sh".format(now_time)
        self.logger.info(cmd)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash'+"\n")
            w.write(cmd)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path))
        shell = "/bioinfo/tool_lab/Yoogene/script_temp/{}".format(
            os.path.basename(file_path))
        self.logger.info("开始生成fastq压缩文件")
        self.logger.info(shell)
        command1 =self.add_command("samstat_{}".format(out_name.lower()), shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行samtools完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        os.system('rm {}'.format(file_path))
        # command = self.add_command("samstat", cmd, ignore_error = True)
        # command.run()
        # self.wait(command)
        # if command.return_code == 0:
        #     self.logger.info("运行samtools stat完成")
        # else:
        #     self.logger.info("运行samtools stat失败")


    def run_samdepth(self):
        bam_file = self.option("bam_file").prop['path']
        # output = self.option("output_name")
        bed = self.posdb_bed
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.software_dir + '/bioinfo/tool_lab/Yoogene/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "samdepth_{}.sh".format(now_time)
        cmd = "{} depth -a -b {} {} > bam_depth_statistic".format(self.samtools,bed, bam_file)
        # command = self.add_command("samdepth", cmd, ignore_error = True)
        # command.run()
        # self.wait(command)
        # if command.return_code == 0:
        #     self.logger.info("运行samtools depth 完成")
        # else:
        #     self.logger.info("运行samtools depth失败")
        self.logger.info(cmd)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash'+"\n")
            w.write(cmd)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path))
        shell = "/bioinfo/tool_lab/Yoogene/script_temp/{}".format(
            os.path.basename(file_path))
        self.logger.info("开始生成fastq压缩文件")
        self.logger.info(shell)
        command1 =self.add_command("samsdepth", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行samtools完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        os.system('rm {}'.format(file_path))

    def run_sample_dep(self):
        self.vcf_yfull = self.option("vcf_yfull").prop['path']
        depget= 0
        depget_Y = 0
        chrY = 0
        CHR = 0 
        y_1 = 0
        y_2 = 0
        y_5 = 0
        y_7 = 0
        y_10 = 0
        y_15 = 0
        y_20 = 0
        chrydep1num = 0
        with open(self.vcf_yfull,"r") as vy:
            while 1:
                line = vy.readline()
                if line[0:2] == '##':
                    pass
                elif line[0] == '#':
                    break

            while 1 :
                line = vy.readline()
                if not line:
                    break
                fields = line.rstrip().split('\t')
                if fields[7].split(";")[0] == 'INDEL':
                    continue
                CHR += 1
                info = fields[7].split(";")
                sample_vcf = {}
                for it in info:
                    tt = it.split("=")
                    sample_vcf[tt[0]] = tt[1]
                nn = "DP"
                if fields[0] == "chrY" and int(sample_vcf[nn]) == 1:
                    chrydep1num += 1
                if fields[0] == "chrY" and int(sample_vcf[nn]) >= 1:
                    y_1 += 1
                if fields[0] == "chrY" and int(sample_vcf[nn]) >= 2:
                    y_2 += 1 
                if fields[0] == "chrY" and int(sample_vcf[nn]) >= 5:
                    y_5 += 1
                if fields[0] == "chrY" and int(sample_vcf[nn]) >= 7:
                    y_7 += 1
                if fields[0] == "chrY" and int(sample_vcf[nn]) >= 10:
                    y_10 += 1
                if fields[0] == "chrY" and int(sample_vcf[nn]) >= 15:
                    y_15 += 1
                if fields[0] == "chrY" and int(sample_vcf[nn]) >= 20:
                    y_20 += 1
                if fields[0] == "chrY":
                    chrY += 1
        with open("{}_sample_dep.xls".format(self.sample_name),"w") as sd:
            sd.write("sample\tpos\t1X\t2X\t5X\t7X\t10X\t15X\t20X\n")
            sd.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                self.sample_name,chrY,y_1,y_2,y_5,y_7,y_10,y_15,y_20
            ))


    def run_sample_MAPOT(self):
        map_base = ""
        dedup_map_base = ""
        dedup_ratio = 0
        if os.path.exists("bam_statistic") and os.path.exists("bam_statistic_pre_dedup") and os.path.exists("YFULL_target_result.sample_summary"):
            with open("bam_statistic_pre_dedup",'r') as pre_stat:
                # line = pre_stat.readline()
                while 1 :
                    line = pre_stat.readline()
                    if re.match(r'SN\tbases mapped:\t(\d+)',line):
                        map_base = re.findall(r"\d+",line)[0]
                        break
                    if not line:
                        break
            with open("bam_statistic","r") as dedup_stat:
                while 1 :
                    line = dedup_stat.readline()
                    if re.match(r'SN\tbases mapped:\t(\d+)',line):
                        dedup_map_base = re.findall(r"\d+",line)[0]
                        break
                    if not line:
                        break
            if float(map_base) == 0:
                dedup_ratio = "Null"
            else:
                dedup_ratio = float(dedup_map_base)/float(map_base)*100

            bed_base = 0
            OT = 0
            with open("YFULL_target_result.sample_summary","r") as Y_ss:
                Y_ss.readline()
                line = Y_ss.readline()
                fields = line.rstrip().split("\t")
                bed_base = float(fields[1])
                OT = bed_base/float(map_base)*100
            # clean_stat = xlrd.open_workbook(self.option("clean_statistic"))
            # table = clean_stat.sheet_by_index(0)
            # nrows = table.nrows
            # clean_stat = 0
            # for i in range(nrows):
            #     line = table.row_values(i)
            #     if line[0] == "CLEAN":
            #         clean_base = line[2]
            # mapping_ratio = float(map_base)/float(clean_base)*100
            clean_base = 0
            with open(self.option("clean_statistic"),"r") as cs:
                while 1:
                    line = cs.readline()
                    if not line:
                        self.logger.info("clean_statistic no content")
                        break
                    fields = line.rstrip().split('\t')
                    if fields[0] == "CLEAN":
                        clean_base = float(fields[2])
            if float(clean_base) == 0:
                mapping_ratio = "Null"
            else:
                mapping_ratio = float(map_base)/float(clean_base)*100
            with open("{}_sample_MAP_OT.xls".format(self.sample_name),"w") as smo:
                smo_text =  "sample\tclean_bp\tmapped_bp\tontarget_bp\tmapping_ratio\tOT\tdedup_ratio\n"
                smo_text = smo_text + "{}\t{}\t{}\t{}\t{}%\t{}%\t{}%\n".format(
                    self.sample_name,clean_base,map_base,bed_base,mapping_ratio,OT,dedup_ratio
                )
                smo.write(smo_text)
            # with open("clean_statistic_reform.xls","r") as csr:
    
    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if len(self.output_dir) > 0:
                shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.work_dir,"{}_sample_MAP_OT.xls".format(self.sample_name)),os.path.join(self.output_dir,"{}_sample_MAP_OT.xls".format(self.sample_name)))
            os.link(os.path.join(self.work_dir,"{}_sample_dep.xls".format(self.sample_name)),os.path.join(self.output_dir,"{}_sample_dep.xls".format(self.sample_name)))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))

        

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "mapstat_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.mapstat",
            "options": dict(
                bam_file = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/2_bwa/align_final_sort.bam",
                bam_prermdup = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/QC_and_mapping_temp/aln-pe.sort.bam",
                clean_statistic = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/1_QC/clean_statistic_reform.xls",
                vcf_yfull = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/Y_full_vcf/align_final_sort.vcf",
                sample_name = "YG202003761",
            #     posdb_bed = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/file/posdb_update_20200801.bed",
            #     EX_autosnp_bed = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/file/EX_autosnp2_2.bed",
            #     chrY_bed = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/file/Y_chr_hg38.bed",
            #     ref = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/file/GRCh38/Homo_sapiens_assembly38.fasta"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()