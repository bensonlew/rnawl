# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from __future__ import division
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import subprocess

class GenomeStructureAgent(Agent):
    """
    获得染色体上基因数量、GC含量信息
    """
    def __init__(self, parent):
        super(GenomeStructureAgent, self).__init__(parent)
        options = [
            {"name": "in_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "in_gff", "type": "infile", "format": "gene_structure.gff3"},
            {"name": "in_gtf", "type": "infile", "format":"gene_structure.gtf"},
            {"name": "ref_genome", "type": "string", "default": ""}
        ]
        self.add_option(options)
        self.step.add_steps("genome_structure")
        self.on('start', self.start_genome_structure)
        self.on("end", self.end_genome_structure)

    def start_genome_structure(self):
        self.step.genome_structure.start()
        self.step.update()

    def end_genome_structure(self):
        self.step.genome_structure.finish()
        self.step.update()

    def check_options(self):
        if not self.option("in_fasta").is_set:
            raise OptionError("参数in_fasta不能为空")

    def set_resource(self):
        self._cpu = 4
        self._memory = "6G"


class GenomeStructureTool(Tool):
    def __init__(self, config):
        super(GenomeStructureTool, self).__init__(config)
        self.seqkit_path = self.config.SOFTWARE_DIR + "/bioinfo/seq"
        self.bedtools = self.config.SOFTWARE_DIR + "/bioinfo/rna/bedtools2-master/bin/bedtools"
        self.scripts = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts"
        self.bioawk_path = self.config.SOFTWARE_DIR + "/bioinfo/seq/bioawk"

    def cmd1(self):
        cmd = "less %s | grep \"^#\" -v | awk -F '\\t' " \
              "-vOFS='\\t' '$9 ~ /ID=transcript/ {print $1,$9}' |  sed 's/\\t.*biotype=/\\t/g' | " \
              "cut -d ';' -f1 | awk -F '\\t' -vOFS='\\t' '{chr[$1];Gene[$1]++; " \
              "if($2 == \"protein_coding\"){Pod[$1]++}else if($2 ~/RNA/){Rna[$1]++}" \
              "else if($2 ~ /pseudogene/){Pse[$1]++}}END{for(i in Gene) " \
              "{print i,Gene[i],Pod[i],Rna[i],Pse[i]}}'  |sort -k1,1 | " \
              "awk -F '\\t' -vOFS='\\t' 'BEGIN{print \"Chr\",\"Gene\",\"ProteinCoding\",\"OtherRNA\"," \
              "\"Pseudogene\"}1' > %s/gene.content.tab.xls" % (self.option("in_gff").prop["path"], self.output_dir)
        os.system(cmd)

    def cmd2(self):
        cmd = "%s/seqkit fx2tab -n -i -g -l  %s" \
              " | awk  -vOFS='\\t' 'BEGIN{print \"Chr\",\"Size(Mb)\",\"GC\"}{print $1,int($2/10000+0.5)/100,$3}' " \
              "| perl %s/tabletools_add.pl -i %s/gene.content.tab.xls -t - -n 1 > %s/gene.stat.xls" % \
              (self.seqkit_path, self.option("in_fasta").prop["path"], self.scripts, self.output_dir, self.output_dir)
        self.logger.info(cmd)
        os.system(cmd)

    def cmd3(self):
        cmd = "%s/bioawk  -c fastx 'BEGIN{print \"Ref seq\tLength\tGC\"};{print $name\"\\t\"length($seq)\"\\t\"gc($seq)*100}' %s > %s/gene.stat.xls" % (self.bioawk_path, self.option("in_fasta").prop["path"], self.output_dir)
        self.logger.info(cmd)
        os.system(cmd)

    def run(self):
        super(GenomeStructureTool, self).run()
        if self.option("ref_genome") != "customer_mode":
            self.cmd1()
            self.cmd2()
        else:
            self.cmd3()
        self.end()
