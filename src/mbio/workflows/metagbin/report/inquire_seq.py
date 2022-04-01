# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __modify__ = '2019.02.28'

from biocluster.workflow import Workflow
import os
import json


class InquireSeqWorkflow(Workflow):
    """
    序列查询
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(InquireSeqWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'method', 'type': 'string'},  #
            {'name': 'query', 'type': 'infile', 'format': 'metagbin.file_gz'},  # 压缩文件
            {'name': 'ref_fasta', 'type': 'string'},  # 参考fasta序列
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.new_taxon = ""
        self.inquire_seq = self.add_tool("metagbin.inquire_seq")
        self.gunzip_fasta = self.add_tool('sequence.fasta_ungz')

    def run(self):
        self.run_ungz()
        super(InquireSeqWorkflow, self).run()

    def run_ungz(self):
        fasta = self.option("query").prop['path']
        self.gunzip_fasta.set_options({
            "fasta": fasta,
            "sample_name": "all",
        })
        self.gunzip_fasta.on('end', self.run_seq)
        self.gunzip_fasta.run()

    def run_seq(self):
        self.writ_fasta()
        method = ''
        if self.option("method") in ['nul']:
            method = 'blastn'
        elif self.option("method") in ['pro']:
            method = 'blastx'
        opts =({
            "in_file": self.gunzip_fasta.option("out_fa"),
            "method": method,
            "ref": self.work_dir + "/ref.fasta",
        })
        self.inquire_seq.set_options(opts)
        self.inquire_seq.on('end',self.set_db)
        self.inquire_seq.run()

    def writ_fasta(self):
        if os.path.exists(self.work_dir + "/ref.fasta"):
            os.remove(self.work_dir + "/ref.fasta")
        file = self.work_dir + "/ref.fasta"
        with open (file,'w') as f:
            f.write(">ref\n")
            f.write(self.option('ref_fasta'))

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        self.api_path = self.api.api('metagbin.inquire_seq')
        inquire_id = self.option("main_id")
        self.api_path.add_anno_nr_detail(inquire_id,self.inquire_seq.output_dir + '/all.blast.xls')
        os.system("sed -i '1i Query-Name\tHit-Name\tIdentity\tHSP-Len\tmismatched\tgap\tQuery_start\tQuery_end\tHit_start\tHit_end\tE-value\tScore' %s" % self.inquire_seq.output_dir + '/all.blast.xls')
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.inquire_seq.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果目录"],
            ['all.blast.xls','xls','查询比对blast结果']
        ])
        super(InquireSeqWorkflow, self).end()