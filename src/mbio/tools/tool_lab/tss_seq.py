# coding=utf-8
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import glob
import shutil
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
__author__ = 'liubinxu'


class TssSeqAgent(Agent):
    """
    predict small_rna TF binding site
    db: vertebrates,plants,insects,fungi
    species: species for trans mir
    pro_dis: promoter distance before translation start site
    """
    def __init__(self, parent):
        super(TssSeqAgent, self).__init__(parent)
        options = [
            {'name': 'tss_up', 'type': 'int', 'default': 5000},
            {'name': 'tss_down', 'type': 'int', 'default': 1000},
            {"name": "ref", "type": "infile", "format": "sequence.fasta"},
            {"name": "gtf", "type": "infile", "format": "ref_rna_v2.gtf"},
            {"name": "bed", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "id_file", "type": "infile", "format": "ref_rna_v2.common"}
        ]
        self.add_option(options)



    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = "{}G".format('30')

    def end(self):
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(TssSeqAgent, self).end()

class TssSeqTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(TssSeqTool, self).__init__(config)
        self.python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/python'
        self.python = '/miniconda2/bin/python'
        self.bedtool_path = 'bioinfo/seq/bedtools-2.25.0/bin/bedtools'

        self.filter_list =None
    def gtf2bed(self):
        self.option("gtf").to_bed()
        bed_path = os.path.split(self.option("gtf").prop['path'])[0]
        bed = os.path.join(bed_path, os.path.split(self.option("gtf").prop['path'])[1] + ".bed")
        self.option("bed", bed)


    def get_filter_list(self, bed_in):
        if self.option("id_file").is_set:
            with open(self.option("id_file").prop['path'], 'r') as f:
                self.filter_list = [l.strip() for l in f]

        with open(bed_in, 'r') as f_in, open(bed_in + 'filter', 'w') as f_out:
            for line in f_in:
                cols = line.split("\t")
                if cols[3] in self.filter_list:
                    f_out.write(line)
        return bed_in + 'filter'



    def get_tss_bed(self):
        """
        获取启动子区域bed文件
        """

        tss_path = self.work_dir + '/tss_range.bed'
        with open(tss_path, 'w') as tss_w:
            with open(self.option("bed").prop['path'], 'r') as known_f:
                for line in known_f:
                    cols = line.strip().split("\t")
                    chro = cols[0]
                    seq_id = cols[3]
                    strand = cols[5]
                    start = cols[1]
                    end = cols[2]
                    if strand == "+":
                        tss_start = str(max(int(start) - int(self.option("tss_up")), 0))
                        tss_end = str(int(start) + int(self.option("tss_down")))

                    else:
                        tss_start = str(max(int(end) - int(self.option("tss_down")), 0))
                        tss_end = str(int(end) + int(self.option("tss_up")))
                        tss_w.write("\t".join([chro, tss_start, tss_end, seq_id, "0", strand]) + "\n")
        return tss_path



    def get_fasta_from_bed(self, bed):
        cmd = "{bedtool_path} getfasta -fi {fna} -bed {tss_bed} -s -name -fo {work_dir}/tss.fa".format(
            bedtool_path=self.bedtool_path, fna=self.option("ref").prop['path'], work_dir = self.work_dir, tss_bed=bed)
        os.system(cmd)
        bed2fa = self.add_command("bedtools", cmd).run()
        self.wait(bed2fa)
        if bed2fa.return_code == 0:
            self.logger.info("bed2fa")
        else:
            self.set_error("bed2fa 运行错误")



    def set_output(self):
        if os.path.exists(self.output_dir + '/tss.fa'):
            os.remove(self.output_dir + '/tss.fa')
        os.link(self.work_dir + '/tss.fa', self.output_dir + '/tss.fa')

    def run(self):
        super(TssSeqTool, self).run()
        if self.option('gtf').is_set:
            self.gtf2bed()
        bed = self.get_tss_bed()
        if self.option("id_file").is_set:
            bed = self.get_filter_list(bed)
        self.get_fasta_from_bed(bed)
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/data5'
        data = {
            "id": "tf_binding" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "tool_lab.tss_seq",
            "instant": False,
            "options": dict(
                ref = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/Ensemble_release_36/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa",
                known_pre = test_dir + "/" + "known_mirna_detail",
                novol_pre = test_dir + "/" + "novel_miR_mature_infor.xls",
                species = "ath",
                db = "plant"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
