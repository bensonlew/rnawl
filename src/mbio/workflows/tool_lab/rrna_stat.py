# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from mbio.files.sequence.file_sample import FileSampleFile
from biocluster.core.exceptions import OptionError
from mbio.packages.tool_lab.file_compress.file_compress import FileCompress
import gevent
from mbio.packages.rna.annot_config import AnnotConfig
import glob
import shutil


class RrnaStatWorkflow(Workflow):
    """
    功能: 统计fastq序列中核糖体rna的含量
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RrnaStatWorkflow, self).__init__(wsheet_object)
        options = [
            # {"name": "fastq_l", "type": "infile", "format": 'ref_rna_v2.common'},  # 原始输入文件,fastq R1端文件
            # {"name": "fastq_r", "type": "infile", "format": 'ref_rna_v2.common'},  # 原始输入文件,fastq R2端文件
            # {"name": "fastq_s", "type": "infile", "format": 'ref_rna_v2.common'},  # 原始输入文件,fastq 单端文件
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'ref_rna_v2.common_dir'},#
            {"name": "seq_type", "type": "string", "default": "fastq"},  # 评估软件 ["bowtie2","blast"]
            {"name": "assess_soft", "type": "string", "default": "bowtie2"},  # 评估软件 ["bowtie2","blast"]
            {"name": "assess_method", "type": "string", "default": "both"},  # 评估方式 ["R1","R2","both"]
            {"name": "seq_num", "type": "int", "default": 1000},  # 抽取评估序列长度
            {"name": "rfam_version", "type": "string", "default": "14.1"},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.if_gz = 0
        self.fq_type = "PE"
        self.samples = {}
        self.tools =[]
        self.gz2fastq = self.add_tool('ref_rna_v2.gzfastq2fastq')
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def ungzfastq(self):
        self.logger.info('需要解压fastq文件，正在解压')
        self.gz2fastq.set_options({'fastq_path': self.option('fastq_dir').prop["path"]})
        self.gz2fastq.on("end", self.stat_run)
        self.gz2fastq.run()
        self.logger.info('解压完成')

    def check_options(self):
        if  self.option("seq_type") == "fasta" and self.option("assess_soft") == "bowtie2":
            raise OptionError("fasta文件仅可使用blast比对")
        if not self.option("fastq_dir").is_set:
            raise OptionError("需要传入fastq文件或者文件夹")
        if self.option("fastq_dir").is_set:
            list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            if not os.path.exists(list_path):
                OptionError("缺少list文件")
            with open(list_path, 'r') as list_r:
                list_info = list_r.read()
                if u'.gz' in list_info:
                    self.if_gz = 1
                file_num = len(list_info.strip().split("\n"))
            self.samples = self.get_list()
            self.logger.info("共有{}个文件".format(str(file_num)))
            self.logger.info("共有{}个样本".format(str(len(self.samples))))
            if file_num == len(self.samples) * 2:
                self.fq_type = "PE"
            else:
                self.fq_type = "SE"
            row_num = len(open(list_path, "r").readline().split())
            self.logger.info(row_num)
            if self.fq_type == "PE" and row_num != 3:
                raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code = "25000603")
            if self.fq_type == "PE":
                for s in sorted(self.samples):
                    if self.samples[s]["l"].split(".")[-1] in ["gz"]:
                        self.samples[s]["l"] = ".".join(self.samples[s]["l"].split(".")[:-2]) + ".fastq"
                    if self.samples[s]["r"].split(".")[-1] in ["gz"]:
                        self.samples[s]["r"] = ".".join(self.samples[s]["r"].split(".")[:-2]) + ".fastq"
            if self.fq_type == "SE" and self.option("seq_type") == "fastq":
                for s in sorted(self.samples):
                    if self.samples[s]["s"].split(".")[-1] in ["gz"]:
                        self.samples[s]["s"] = ".".join(self.samples[s].split(".")[:-2]) + ".fastq"
            print(self.samples)
        return True

    def get_list(self):
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        self.logger.info(samples)
        return samples

    def run_tools(self):
        if self.option("assess_soft").lower() == "blast":
            self.rfam_run()
        if self.option("assess_soft").lower() == "bowtie2":
            self.rrna_run()

        self.logger.info('{}'.format(self.events))
        for eve in self.events.values():
            self.logger.info('{}'.format(eve.is_start))

    def run(self):
        if self.if_gz:
            self.ungzfastq()
        else:
            self.run_tools()
        super(RrnaStatWorkflow, self).run()

    def rrna_run(self):
        rfam = self.pfam_db = AnnotConfig().get_file_path(
            file="Rfam",
            db="rfam",
            version=self.option("rfam_version"))
        rfam_path = os.path.join(os.path.dirname(rfam),"Rfam_bowtie2_index/Rfam")
        n = 0
        for f in sorted(self.samples):
            options = {}
            mapping_tool = self.add_tool('tool_lab.rfam_stat.bowtie2')
            if self.fq_type == "PE":
                fq_l = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["l"])
                fq_r = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["r"])
                options = {
                    "ref_genome":rfam_path,
                    'left_reads': fq_l,
                    'right_reads': fq_r,
                    'sample': f,
                    "seq_method": self.fq_type,
                    'seq_num' : self.option('seq_num')
                }
            elif self.fq_type == "SE":
                fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["s"])
                options = {
                    "ref_genome": rfam_path,
                    'single_end_reads': fq_s,
                    'sample': f,
                    "seq_method": self.fq_type,
                    'seq_num': self.option('seq_num')
                }

            self.logger.info("tool {} prok_rna.bowtie2".format(options))
            mapping_tool.set_options(options)
            self.tools.append(mapping_tool)
            # mapping_tool.on("end", self.finish_update, "mapping_{}".format(n))
            # mapping_tool.run()
            n += 1
        self.on_rely(self.tools, self.set_output)
        for mapping_tool in self.tools:
            mapping_tool.run()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def rfam_run(self):
        for n, sample_name in enumerate(self.samples):
            k = 'l' if self.fq_type == 'PE' else 's'
            rfam = self.add_tool('tool_lab.rfam_stat.blast_to_rfam')
            if self.option("seq_type") == "fastq":
                rfam.set_options({
                    'query': os.path.join(self.option('fastq_dir').path, self.samples[sample_name][k]),
                    'sample_name': sample_name,
                    'rfam_version': self.option('rfam_version'),
                    'seq_type' : self.option('seq_type') ,
                    'seq_num' :self.option('seq_num')
                })
            elif self.option("seq_type") == "fasta":
                rfam.set_options({
                    'query': os.path.join(self.option('fastq_dir').path, self.samples[sample_name]),
                    'sample_name': sample_name,
                    'rfam_version': self.option('rfam_version'),
                    'seq_type': self.option('seq_type'),
                    'seq_num': self.option('seq_num')
                })
            # rfam.run()
            self.tools.append(rfam)
        self.on_rely(self.tools, self.set_output)
        for rfam in self.tools:
            rfam.run()

    def set_output(self):
        self.logger.info("set output")
        for f in glob.glob(r"{}/*".format(self.output_dir)):
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)
        rfam_dir = os.path.join(self.output_dir, "RfamStat")
        if self.option("assess_soft").lower() == "blast":
            os.mkdir(rfam_dir)
        rrna_dir = os.path.join(self.output_dir, "RrnaStat")
        if self.option("assess_soft").lower() == "bowtie2":
            os.mkdir(rrna_dir)
        for tool in self.tools:
            out_files = os.listdir(tool.output_dir)
            for f in out_files:
                f_path = os.path.join(tool.output_dir, f)
                if "rfam" in f:
                    target = os.path.join(rfam_dir, f)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(f_path, target)
                elif ".stat" in f:
                    target = os.path.join(rrna_dir, f)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(f_path, target)
                else:
                    pass
        rfam_stat = os.path.join(self.output_dir,"Results.xls")
        if self.option("assess_soft").lower() == "blast":
            rrna_stat = glob.glob(self.output_dir + '/RfamStat/' + "*_vs_rfam.stat.xls")
            with open(rfam_stat, "w") as w:
                w.write("Sample" + "\t" + "rRNA(%)" + "\n")
                for fs in sorted(rrna_stat):
                    f = open(fs, "r")
                    f.readline()
                    items = f.readline().strip().split("\t")
                    specimen_name = items[0]
                    r_rna_ratio = items[1]
                    w.write(specimen_name + "\t" + r_rna_ratio.split("%")[0] + "\n")
            shutil.rmtree(self.output_dir + '/RfamStat')
        elif self.option("assess_soft").lower() == "bowtie2":
            rrna_stat = glob.glob(self.output_dir + '/RrnaStat/' + "*.stat")
            with open(rfam_stat, "w") as w:
                w.write("Sample" + "\t" + "rRNA(%)" + "\n")
                for fs in sorted(rrna_stat):
                    sample_name = "".join(os.path.basename(fs).split(".")[:-1])
                    f = open(fs, "r")
                    for line in f.readlines():
                        if "overall alignment rate" in line:
                            ratio = line.strip().split()[0]
                            break
                    w.write(sample_name + "\t" + ratio.split("%")[0] + "\n")
            shutil.rmtree(self.output_dir + '/RrnaStat')
        self.logger.info("done")
        self.set_db()
        # self.end()

    def set_db(self):
        self.database = self.api.api('tool_lab.rrna_stat')
        self.database.add_rrna_stat(
            self.option("main_id"),
            os.path.join(self.output_dir, "Results.xls"),
        )
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "核糖体rna评估结果目录",0],
            ['Results.xls', '', 'Rfam核糖体RNA占比结果表', 0],
        ])
        super(RrnaStatWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.rrna_stat import RrnaStatWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            "id": "rfam_stat_" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.rrna_stat",
            "options": dict(
                fastq_dir="/mnt/ilustre/users/sanger-dev/workspace/20210203/Refrna_tsg_249609/remote_input/fastq_dir/rawdata",
                seq_type ="fastq",
                assess_soft ="blast",
                # assess_method = "",
                seq_num="1000",
            )
        }
        wsheet = Sheet(data=data)
        wf =RrnaStatWorkflow(wsheet)
        wf.sheet.id = 'pdf2image'
        wf.sheet.project_sn = 'pdf2image'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
