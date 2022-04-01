# -*- coding: utf-8 -*-
# __author__ = 'zhengyuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import os

class CgcTaxonAgent(Agent):

    def __init__(self, parent):
        super(CgcTaxonAgent, self).__init__(parent)
        options = [
            {"name": "genetable", "type": "infile", "format": "sequence.profile_table"},    # 输入文件，gene列表
            {"name": "blastoutput", "type": "infile", "format": "align.blast.blast_table"},   # 输入文件，blast的输出结果表
            {"name": "taxonfile", "type": "outfile", "format": "annotation.nr.nr_taxon"}   # 输出文件，带上注释信息的gene表
        ]
        self.add_option(options)
        self.step.add_steps("cgc_taxon")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def check_options(self):
        if self.option("genetable").is_set or self.option("blastoutput").is_set:
            pass
        else:
            raise OptionError("必须输入gene表或blast结果表")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def stepstart(self):
        self.step.cgc_taxon.start()
        self.step.update()

    def stepfinish(self):
        self.step.cgc_taxon.finish()
        self.step.update()

    def end(self):
        '''
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['annotation_result.xls', 'xls', 'gene对应注释文件']
        ])
        '''
        super(CgcTaxonAgent, self).end()

class CgcTaxonTool(Tool):

    def __init__(self, config):
        super(CgcTaxonTool, self).__init__(config)
        self.ref_db = os.path.join(self.config.SOFTWARE_DIR, "database/Cgc_taxon/reference_cgc.xls")

    def get_gene(self):
        genefile = self.option("genetable").prop['path']
        # gene_lst = []
        gene_lst = {}
        with open(genefile, "rb") as queryin:
            for gene in queryin:
                gene = gene.strip().split("\t")
                if gene[0] not in gene_lst:
                    # gene_lst.append(gene)
                    gene_lst[gene[0]] = gene[1]
        return gene_lst

    def blast_gene(self):
        genefile = self.option("blastoutput").prop['path']
        blast_line = defaultdict(list)
        with open(genefile, "rb") as blastout:
            for line in blastout:
                line = line.strip()
                tmp = line.split("\t")
                gene = tmp[10]
                line = [tmp[5], tmp[10]]
                line.extend(tmp[0:5])
                line.extend(tmp[6:10])
                line.extend(tmp[11:14])
                blast_line[gene].append("\t".join(line))
        return blast_line

    def gene_annotation(self):
        outfile = os.path.join(self.output_dir, "annotation_result.xls")
        gene_lst = self.get_gene()
        with open(self.ref_db, "rb") as ref:
            head = ref.next()
            head = "\t".join(head.split("\t")[1:])
            head = "Gene_id\tAbund\t" + head
            with open(outfile, "wb") as output:
                output.write(head)
                for line in ref:
                    line = line.strip()
                    attrs = line.split("\t")
                    gene = attrs[0]
                    content = "\t".join(attrs[1:])
                    # if gene in gene_lst:
                    #     output.write(gene + "\t" + content + "\n")
                    if gene_lst.has_key(gene):
                        output.write(gene + "\t" + gene_lst[gene] + "\t" + content + "\n")

    def blast_annotation(self):
        outfile = os.path.join(self.output_dir, "annotation_result.xls")
        gene_set = self.blast_gene()
        gene_lst = gene_set.keys()
        with open(self.ref_db, "rb") as ref:
            head = ref.next()
            head = "\t".join(head.split("\t")[1:])
            head = gene_set['Hit-Name'][0] + "\t" + head
            with open("annotation_tmp", "wb") as output:
                output.write(head)
                for line in ref:
                    line = line.strip()
                    attrs = line.split("\t")
                    gene = attrs[0]
                    content = "\t".join(attrs[1:])
                    if gene in gene_lst:
                        for bline in gene_set[gene]:
                            output.write(bline + "\t" + content + "\n")
        self.sort_query("annotation_tmp", outfile)
        os.remove("annotation_tmp")

    def sort_query(self, infile, outfile):
        store = defaultdict(list)
        with open(infile, "r") as inf:
            head = inf.next()
            for line in inf:
                line = line.strip()
                tmp = line.split("\t")
                query = tmp[0]
                store[query].append(line)
        with open(outfile, "w") as outf:
            outf.write(head)
            for q in sorted(store.keys()):
                for l in store[q]:
                    outf.write(l + "\n")

    def run(self):
        super(CgcTaxonTool, self).run()
        if self.option("genetable").is_set:
            self.gene_annotation()
        else:
            self.blast_annotation()
        self.end()




