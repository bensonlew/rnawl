# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 20190312

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class TransgeneWorkflow(Workflow):
    """
    转基因
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TransgeneWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "ref_fa", "type": "string"},  # ref.fa
            {"name": "insertion_id", "type": "string"},  # insert.fa中的insertion_id
            {"name": "sequence", "type": "string"},  # insert.fa中的sequence
            {"name": "rawdata_path_1", "type": "string"},  # rawdata路径
            {"name": "rawdata_path_2", "type": "string"},  # rawdata路径
            {"name": "sample", "type": "string"},  # 样本名称
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.aimhii = self.add_tool("wgs_v2.aimhii")
        self.insert_result = self.add_tool("wgs_v2.insert_result")
        self.fastq_qc = self.add_tool("wgs.fastp")
        self.updown_stream_seq = self.add_tool("wgs_v2.updown_stream_seq")
        self.insert_fa = ""
        self.ref_fa = os.path.join(self.config.SOFTWARE_DIR, ("database/dna_geneome/" + self.option("ref_fa")))

    def check_options(self):
        if not self.option('ref_fa'):
            raise OptionError('必须输入ref_fa', code="14500301")
        if not self.option('insertion_id'):
            raise OptionError('必须输入insertion_id', code="14500302")
        if not self.option('sequence'):
            raise OptionError('必须输入sequence', code="14500302")
        if not self.option('rawdata_path_1'):
            raise OptionError('必须输入rawdata_path_1', code="14500302")
        if not self.option('sample'):
            raise OptionError('必须输入sample', code="14500302")
        if not self.option('main_id'):
            raise OptionError('必须输入main_id', code="14500303")
        if not self.option('task_id'):
            raise OptionError('必须输入task_id', code="14500304")
        return True

    def run_fastq_qc(self):
        """
        质控
        """
        options = {
            "fq1": self.option('rawdata_path_1'),
            "fq2": self.option('rawdata_path_2'),
            "qualified_quality_phred": "20",
            "length_required": "36",
            "cut_by_quality5": "20",
            "cut_by_quality3": "3",
            "cut_mean_quality": "20",
            "n_base_limit": "10",
            "thread": "8",
            "compression": "6",
            "sample_name": self.option('sample'),
            "adapter_sequence": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "adapter_sequence_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        }
        self.fastq_qc.set_options(options)
        self.fastq_qc.on("end", self.run_aimhii)
        self.fastq_qc.on("end", self.set_output, "fastq_qc")
        self.fastq_qc.run()

    def make_file(self):
        """
        组建insert.fa文件
        :return:
        """
        with open(os.path.join(self.work_dir, "insert.fa"), "w")as fw:
            write_lines = ">" + self.option('insertion_id') + "\n" + self.option('sequence')
            fw.write(write_lines)
        self.insert_fa = os.path.join(self.work_dir, "insert.fa")

    def run_aimhii(self):
        clean_fq_1 = os.path.join(self.fastq_qc.output_dir, (self.option("sample") + ".clean.1.fastq.gz"))
        clean_fq_2 = os.path.join(self.fastq_qc.output_dir, (self.option("sample") + ".clean.2.fastq.gz"))
        options = {
            "ref_fa": self.ref_fa,
            "insert_fa": self.insert_fa,
            "clean_fq_1": clean_fq_1,
            "clean_fq_2": clean_fq_2,
            "sample": self.option("sample")
        }
        self.aimhii.set_options(options)
        self.aimhii.on("end", self.run_insert_result)
        self.aimhii.on("end", self.set_output, "aimhii")
        self.aimhii.run()

    def run_insert_result(self):
        with open(os.path.join(self.aimhii.output_dir, (self.option('sample') + ".result.csv")), "r")as fr:
            lines = fr.readlines()
            if len(lines) < 2:
                self.set_error('转基因结果为空!')
        options = {
            "result_csv": self.aimhii.option("result_csv").prop["path"],
            "sample": self.option("sample")
        }
        self.insert_result.set_options(options)
        self.insert_result.on("end", self.run_updown_stream_seq)
        self.insert_result.on("end", self.set_output, "insert_result")
        self.insert_result.run()

    def run_updown_stream_seq(self):
        insert_xls = os.path.join(self.insert_result.output_dir, (self.option("sample") + ".insert.xls"))
        options = {
            "insert_xls": insert_xls,
            "ref_fa": self.ref_fa,
            "sample": self.option("sample")
        }
        self.updown_stream_seq.set_options(options)
        self.updown_stream_seq.on("end",  self.set_db)
        self.updown_stream_seq.on("end", self.set_output, "updown_stream_seq")
        self.updown_stream_seq.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'fastq_qc':
            self.linkdir(obj.output_dir, 'fastq_qc')
        if event['data'] == 'aimhii':
            self.linkdir(obj.output_dir, 'aimhii')
        if event['data'] == 'insert_result':
            self.linkdir(obj.output_dir, 'insert_result')
        if event['data'] == 'updown_stream_seq':
            self.linkdir(obj.output_dir, 'updown_stream_seq')

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        self.logger.info("保存结果到mongo")
        transgene_api = self.api.api('wgs_v2.transgene')
        main_id = self.option("main_id")
        file = os.path.join(self.insert_result.output_dir, (self.option("sample") + ".insert.xls"))
        ref_file = self.ref_fa
        self.logger.info("开始进行transgene导表")
        transgene_api.add_transgene(main_id, file, self.insert_fa, ref_file, os.path.join(self.output_dir, "aimhii/merged.sam"))
        self.logger.info("cnv长度统计导表成功！")
        self.end()

    def run(self):
        self.make_file()
        self.run_fastq_qc()
        super(TransgeneWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(TransgeneWorkflow, self).end()
