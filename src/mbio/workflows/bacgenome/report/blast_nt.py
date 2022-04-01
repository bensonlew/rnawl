# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __modify__ = '2020/04/02'

from biocluster.workflow import Workflow
import os
from mbio.packages.bacgenome.common import fasta_cutoff,get_scaffoldlist
import pandas as pd

class BlastNtWorkflow(Workflow):
    """
    扫描图的scaffold的注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BlastNtWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "samp", "type": "string"},
            {"name": "min_len", "type": "string", "default":"1000"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.blast = self.add_tool('align.blast')

    def run(self):
        self.run_blast()
        super(BlastNtWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api('bac_assem.blast')
        self.logger.info("开始进行导表")
        table = os.path.join(self.output_dir, self.option('samp') + ".nt.xls")
        listseqid = get_scaffoldlist(self.work_dir + "/" + "seq.fa")
        api_path.add_blast_nt_detail("blast_nt_detail", self.option("main_id"), table, listseqid)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果目录"],
        ])
        super(BlastNtWorkflow, self).end()

    def run_blast(self):
        fasta_cutoff(self.option('fa').prop['path'], self.option("min_len"), self.work_dir + "/" + "seq.fa")
        opts = {
            'query': self.work_dir + "/" + "seq.fa",
            'query_type': "nucl",
            'database': 'nt_bac',
            'outfmt': 6,
            'blast': 'blastn',
            'memory': 50,
            'num_alignment': 1
        }
        self.blast.set_options(opts)
        self.blast.on("end", self.set_output)
        self.blast.run()

    def set_output(self):
        if os.path.exists(self.output_dir + "/" + self.option('samp') + ".nt.xls"):
            os.remove(self.output_dir + "/" + self.option('samp') + ".nt.xls")
        data = pd.read_table(self.blast.output_dir + "/seq_vs_nt_bac.xls", header=0,
                             names=["Score", "Evalue", "align_length",
                                    "Identity (%)", "Similarity (%)", "Query_id", "q.length", "q.start", "q.end",
                                    "q.direction", "hit_name","s.length", "s.start", "s.end", "s.direction", "Description"])
        data["Coverage (%)"] = (data["q.end"] - data["q.start"]) / data["q.length"].astype("float") * 100
        data["Subject ID"] = data.apply(lambda x: x["Description"].split(" ")[0], axis=1)
        data["Sample Name"] = self.option("samp")
        self.blast_data = data.reindex(
            columns=["Sample Name", "Query_id", "align_length", "q.length", "q.start", "q.end",
                     "q.direction", "Subject ID", "s.length", "s.start", "s.end", "s.direction",
                     "Description", "Identity (%)", "Coverage (%)", "Evalue", "Score"])
        table_result = os.path.join(os.path.join(self.output_dir, self.option("samp") + ".nt.xls"))
        self.blast_data.to_csv(table_result, sep="\t", index=0)
        self.set_db()
