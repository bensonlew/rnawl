# -*- coding: utf-8 -*-
# __author__ = "qinjincheng"

from biocluster.workflow import Workflow
import os
import time
import pandas as pd
import io

class GenesetSelfWorkflow(Workflow):
    """
    Read gene_list file, creat document in sg_geneset and sg_geneset_detail
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetSelfWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name="genes", type="infile", format="lnc_rna.common"),
            dict(name="file_path", type="infile", format="lnc_rna.common"),
            dict(name="name", type="string", default=None),
            dict(name="gene_type", type="string"),
            dict(name="task_id", type="string"),
            dict(name="geneset_id", type="string"),
            dict(name="update_info", type="string"),
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.start_listener()
        if self.option("genes").is_set:
            self.match_db(self.option('file_path').prop['path'], self.option('genes').prop['path'], self.option("gene_type"))
        self.set_db()

    def match_db(self, file, genes, gene_type):
        """
        Get geneset intersection between uploaded file and the specific sg_annotation_query related set.
        :param file: file contains geneset list uploaded by client (header: None)
        :param genes: file contains geneset list built by to_file (header: mirna/gene)
        :param gene_type: M or G
        :return: intersection geneset (list)
        """
        query = pd.read_table(genes, header=0,sep="\t")

        query_list = list(query.iloc[:,0])
        while "" in query_list:
             query_list.remove("")

        with open(os.path.join(self.work_dir, "geneset_self.txt"), "w") as f, io.open(file, "r",encoding='UTF-8-sig') as genef:
            f.write("gene_list\n")
            gene_list=list()
            for gene in genef.readlines():
                gene = gene.strip()
                if gene in query_list:
                    if gene not in gene_list:
                        f.write("{}\n".format(gene))
                        gene_list.append(gene)
        self.set_output()
        return gene_list

    def set_output(self):
        """
        Link geneset_self.txt to output dir.
        :return: None
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        geneset_file = os.path.join(self.work_dir, "geneset_self.txt")
        os.link(geneset_file, os.path.join(self.output_dir, "geneset_self.txt"))
        self.logger.info("Finishing to link geneset_self.txt to output dir")

    def set_db(self):
        """
        Export geneset geneset_self.txt to mongo.
        """
        geneset_self = self.api.api("lnc_rna.geneset_self")
        if self.option("genes").is_set:
            geneset_table = os.path.join(self.output_dir, "geneset_self.txt")
        else:
            geneset_table = self.option("file_path").prop['path']
        geneset_self.add_geneset(
            geneset_table=geneset_table,
            task_id=self.option("task_id"),
            name=self.option("name"),
            gene_type=self.option("gene_type"),
            main_id=self.option("geneset_id")
        )
        self.end()

    def end(self):
        super(GenesetSelfWorkflow, self).end()
