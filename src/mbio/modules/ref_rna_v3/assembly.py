# -*- coding: utf-8 -*-
# __author__ = "wangzhaoyue,shicaiping,qinjincheng"

import glob
import os
import shutil
import unittest

from biocluster.module import Module


class AssemblyModule(Module):
    def __init__(self, work_id):
        super(AssemblyModule, self).__init__(work_id)
        options = [
            {"name": "sample_bam_dir", "type": "infile", "format": "align.bwa.bam_dir"},
            {"name": "ref_fa", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "fr_stranded", "type": "string", "default": "fr-unstranded"},
            {"name": "strand_direct", "type": "string", "default": "none"},
            {"name": "assemble_method", "type": "string", "default": "stringtie"},

            {"name": "F", "type": "float", "default": 0.1},
            {"name": "fpkm_cut", "type": "float", "default": 1},
            {"name": "min_isoform_fraction", "type": "float", "default": 0.1},

            {"name": "min_coverage", "type": "int", "default": 3},
            {"name": "min_read", "type": "int", "default": 5},
            {"name": "min_iso", "type": "float", "default": 0.1},
            {"name": "min_cov", "type": "int", "default": 5},
            {"name": "min_tpm", "type": "int", "default": 1},

            {"name": "new_transcripts_gtf", "type": "outfile", "format": "ref_rna_v2.gtf"},
            {"name": "ref_and_new_gtf", "type": "outfile", "format": "ref_rna_v2.gtf"},
            {"name": "trans2gene", "type": "outfile", "format": "ref_rna_v2.common"},
            {"name": "new_transcripts_fa", "type": "outfile", "format": "ref_rna_v2.fasta"},
            {"name": "all_transcripts_fa", "type": "outfile", "format": "ref_rna_v2.fasta"},
        ]
        self.add_option(options)
        self.tools = []
        self.tool_dict = dict()

    def check_options(self):
        pass

    def run(self):
        super(AssemblyModule, self).run()
        if self.option("assemble_method").lower() == "cufflinks":
            self.run_cufflinks()
        elif self.option("assemble_method").lower() == "stringtie":
            self.run_stringtie()
        else:
            self.end()

    def run_stringtie(self):
        for sample_bam in glob.glob(os.path.join(self.option("sample_bam_dir").path, "*.bam")):
            stringtie = self.add_tool("ref_rna_v2.assembly.stringtie")
            stringtie.set_options({
                "sample_bam": sample_bam,
                "ref_gtf": self.option("ref_gtf"),
                "fr_stranded": self.option("fr_stranded"),
                "strand_direct": self.option("strand_direct"),
                "min_coverage": self.option("min_coverage"),
                "min_read": self.option("min_read"),
            })
            self.tools.append(stringtie)
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
        else:
            self.on_rely(self.tools, self.run_stringtie_merge)
        for tool in self.tools:
            tool.run()

    def run_stringtie_merge(self):
        gtf_list_fp = os.path.join(self.work_dir, "gtf_list.txt")
        open(gtf_list_fp, "w").writelines("{}\n".format(tool.option("sample_gtf").path) for tool in self.tools)
        stringtie_merge = self.add_tool("ref_rna_v2.assembly.stringtie_merge")
        stringtie_merge.set_options({
            "gtf_list_fp": gtf_list_fp,
            "ref_fa": self.option("ref_fa"),
            "ref_gtf": self.option("ref_gtf"),
            "min_iso": self.option("min_iso"),
            "min_tpm": self.option("min_tpm"),
            "min_cov": self.option("min_cov"),
        })
        stringtie_merge.on("end", self.run_gffcompare)
        stringtie_merge.run()
        self.tool_dict["stringtie_merge"] = stringtie_merge

    def run_gffcompare(self):
        merged_gtf = self.tool_dict["stringtie_merge"].option("merged_gtf")
        gffcompare = self.add_tool("ref_rna_v2.assembly.gffcompare")
        gffcompare.set_options({
            "merged_gtf": merged_gtf,
            "ref_gtf": self.option("ref_gtf"),
        })
        gffcompare.on("end", self.run_new_transcripts)
        gffcompare.run()
        self.tool_dict["gffcompare"] = gffcompare

    def run_cufflinks(self):
        pass

    def run_new_transcripts(self):
        tmap = self.tool_dict["gffcompare"].option("tmap")
        merged_gtf = self.tool_dict["stringtie_merge"].option("merged_gtf")
        new_transcripts = self.add_tool("ref_rna_v2.assembly.new_transcripts")
        new_transcripts.set_options({
            "tmap": tmap,
            "merged_gtf": merged_gtf,
            "ref_gtf": self.option("ref_gtf"),
            "ref_fa": self.option("ref_fa")
        })
        new_transcripts.on("end", self.run_statistics)
        new_transcripts.run()
        self.tool_dict["new_transcripts"] = new_transcripts

    def run_statistics(self):
        statistics = self.add_tool("ref_rna_v2.assembly.statistics")
        statistics.set_options({
            "all_transcripts": self.tool_dict["new_transcripts"].option("all_transcripts"),
            "add_code_merged": self.tool_dict["new_transcripts"].option("add_code_merged"),
            "new_transcripts": self.tool_dict["new_transcripts"].option("new_transcripts"),
            "new_genes": self.tool_dict["new_transcripts"].option("new_genes"),
            "old_transcripts": self.tool_dict["new_transcripts"].option("old_transcripts"),
            "old_genes": self.tool_dict["new_transcripts"].option("old_genes"),
        })
        statistics.on("end", self.set_output)
        statistics.run()
        self.tool_dict["statistics"] = statistics

    def set_output(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        shutil.copytree(self.tool_dict["stringtie_merge"].output_dir, os.path.join(self.output_dir, 'StringtieMerge'))
        shutil.copytree(self.tool_dict["gffcompare"].output_dir, os.path.join(self.output_dir, 'Gffcompare'))
        shutil.copytree(self.tool_dict["new_transcripts"].output_dir, os.path.join(self.output_dir, 'NewTranscripts'))
        shutil.copytree(self.tool_dict["statistics"].output_dir, os.path.join(self.output_dir, 'Statistics'))
        # self.option('new_transcripts_gtf').set_path()
        # self.option('ref_and_new_gtf').set_path()
        # self.option('trans2gene').set_path()
        # self.option('new_transcripts_fa').set_path()
        # self.option('all_transcripts_fa').set_path()
        self.end()

    def end(self):
        super(AssemblyModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data5 = {
            "id": "assembly_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "module",
            "name": "ref_rna_v3.assembly",
            "instant": False,
            "options": dict(
                sample_bam_dir="/mnt/ilustre/users/sanger-dev/workspace/20200313/Refrna_tsg_37039/RnaseqMapping"
                               "/output/bam",
                ref_fa="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana"
                       "/TAIR10_Ensembl_43/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
                ref_gtf="/mnt/ilustre/users/sanger-dev/workspace/20200313/Refrna_tsg_37039/FileCheck"
                        "/Arabidopsis_thaliana.TAIR10.43.gtf",
                fr_stranded="fr-unstranded",
                strand_direct="none",
                assemble_method="stringtie",
            )
        }
        wsheet = Sheet(data=data5)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
