# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/1/9'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import check_command
from mbio.packages.align.pocp import Pocp
import pandas as pd
from mbio.packages.metagbin.common_function import get_ani_species,get_pocp_genus

class TaxStatAgent(Agent):
    """
    version 1.0
    """

    def __init__(self, parent):
        super(TaxStatAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "ref", "type": "infile", "format": "sequence.fasta,sequence.fasta_dir"},
            {"name": "query_blast", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "ref_blast", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "blasr_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "tax_file", "type": "string"},
            {"name": "ani_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "sample", "type": "string"},
            {"name": "analysis", "type": "string"},
            {"name": "task_id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("analysis"):
            raise OptionError("需定义调用的分析内容")
        if self.option("analysis") not in ['blasr', 'ani', 'pocp']:
            raise OptionError("分析不在范围内")
        if self.option("analysis") == "blasr":
            if not self.option("blasr_table").is_set:
                raise OptionError("blasr输入未定义")
            if not self.option("tax_file"):
                raise OptionError("对应物种未定义")
        elif self.option("analysis") == "ani" and not self.option("ani_table").is_set:
            raise OptionError("ani输入未定义")
        elif self.option("analysis") == "pocp":
            if not self.option("query").is_set:
                raise OptionError("pocp分析未定义query参数")
            if not self.option("ref").is_set:
                raise OptionError("pocp分析未定义ref参数")
            if not self.option("query_blast").is_set:
                raise OptionError("pocp分析未定义query_blast参数")
            if not self.option("ref_blast").is_set:
                raise OptionError("pocp分析未定义ref_blast参数")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"


class TaxStatTool(Tool):
    def __init__(self, config):
        super(TaxStatTool, self).__init__(config)

    def run(self):
        super(TaxStatTool, self).run()
        if self.option("analysis") == "pocp":
            self.run_pocp_stat()
        elif self.option("analysis") == "ani":
            self.run_ani_stat()
        elif self.option("analysis") == "blasr":
            self.run_blasr_stat()
        self.set_output()
        self.end()

    def run_pocp_stat(self):
        result = os.path.join(self.output_dir, "pocp_result.xls")
        pocp_obj = Pocp()
        pocp_obj.parse_fa(self.option("query").path)
        pocp_obj.parse_fa(self.option("ref").path)
        pocp_obj.parse_table(self.option("query_blast").path, "m5")
        pocp_obj.parse_table(self.option("ref_blast").path, "m5")
        with open(result, "w") as file:
            file.write("Genome\tReference\tReference_Taxonomy\tPOCP(%)\n")
            genome = os.path.basename(self.option("query").path)
            genome = genome.split("_CDS")[0]
            reference = os.path.basename(self.option("ref").path)
            reference = os.path.splitext(reference)[0]
            reference = reference.split("_CDS")[0]
            genus = get_pocp_genus(self.option("task_id"), reference)
            pocp_value = pocp_obj.get_pocp_value()
            file.write(genome + "\t" + reference + "\t" + genus + "\t%.3f" % pocp_value + "\n")

    def run_blasr_stat(self):
        result = os.path.join(self.output_dir, "blasr_result.xls")
        db_path = os.path.join(self.config.SOFTWARE_DIR, "database/taxon_db")
        if self.option("tax_file") == "silva128":
            tax_path = os.path.join(db_path, self.option("tax_file"), "silva.16s.tax")
        elif self.option("tax_file") == "silva132":
            tax_path = os.path.join(db_path, self.option("tax_file"), "silva.16s.tax")
        elif self.option("tax_file") == "greengene":
            tax_path = os.path.join(db_path, "Greengenes135/greengenes.16s.tax")
        elif self.option("tax_file") == "nt":
            tax_path = os.path.join(db_path, self.option("tax_file"), "nt.tax")
        tax_info = pd.read_table(tax_path, header=None, index_col=0)
        tax_info.index = tax_info.index.astype("str")
        with open(self.option("blasr_table").path, "r") as file, open(result, "w") as file2:
            lines = file.readlines()
            if len(lines) == 0:
                return True
            file2.write("Genome\tTaxonomy\tIdentity(%)\n")
            for line in lines:
                line = line.strip().split()
                self.logger.info(line)
                try:
                    tax = tax_info.loc[line[1]].tolist()[0]
                except:
                    tax = "unknown"
                file2.write("%s\t%s\t%s\n" % (self.option("sample"), tax, line[3]))

    def run_ani_stat(self):
        result = os.path.join(self.output_dir, "ani_result.xls")
        with open(self.option("ani_table").path, "r") as file1, open(result, "w") as file2:
            lines = file1.readlines()
            if len(lines) == 0:
                return True
            file2.write("Genome\tReference\tReference_Taxonomy\tANI(%)\n")
            for line in lines:
                line = line.strip().split()
                query = os.path.basename(line[0])
                reference = os.path.basename(line[1])
                # species = get_ani_species(reference)
                species = get_pocp_genus(self.option("task_id"), reference)
                reference = os.path.splitext(reference)[0]
                reference = reference.split("_CDS")[0]
                ani_value = line[2]
                file2.write(query + "\t" + reference + "\t" + species + "\t%.3f\n" % float(ani_value))

    def set_output(self):
        pass