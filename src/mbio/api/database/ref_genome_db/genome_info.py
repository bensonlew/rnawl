# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
# __date__ = '20200320

from biocluster.api.database.base import Base, report_check
import datetime
import pandas as pd
import pymongo
from mbio.api.database.ref_rna_v2.api_base import ApiBase


class GenomeInfo(ApiBase):
    def __init__(self, bind_object=None):
        super(GenomeInfo, self).__init__(bind_object)
        self._project_type = 'ref_genome_db'

    @report_check
    def add_genome_info(self, client, genome_stat, trinity_stat, level_file, g2t2p, biomart, biomart_type, species_name, ref_anno_version, hyperlink, release_date, remark, genome_id):
        g2t2p_pd = pd.read_table(g2t2p, header=None, usecols=[0,1])
        gene_num = len(set(g2t2p_pd[0]))
        transcript_num = len(set(g2t2p_pd[1]))
        chromosome_num = 0
        contig_num = 0
        scaffold_num = 0
        protein_coding_num = 0
        other_rna_num = 0
        pseudogene_num = 0
        n50 = 0
        gc = 0
        size = 0
        chr = 0
        gene_name_num = 0
        gene_description_num = 0

        biomart_pd = pd.read_table(biomart, header=None)
        if biomart_type == "type1":
            NONE_VIN = (biomart_pd[2].isnull()) | (biomart_pd[2].apply(lambda x: str(x).isspace() or str(x) == "-"))
            gene_name_num = len(set(biomart_pd[~NONE_VIN][0]))
            NONE_VIN = (biomart_pd[7].isnull()) | (biomart_pd[7].apply(lambda x: str(x).isspace() or str(x) == "-"))
            gene_description_num = len(set(biomart_pd[~NONE_VIN][0]))
        elif biomart_type == "type2":
            NONE_VIN = (biomart_pd[2].isnull()) | (biomart_pd[2].apply(lambda x: str(x).isspace() or str(x) == "-"))
            gene_name_num = len(set(biomart_pd[~NONE_VIN][0]))
            NONE_VIN = (biomart_pd[5].isnull()) | (biomart_pd[5].apply(lambda x: str(x).isspace() or str(x) == "-"))
            gene_description_num = len(set(biomart_pd[~NONE_VIN][0]))
        elif biomart_type == "type3":
            NONE_VIN = (biomart_pd[3].isnull()) | (biomart_pd[3].apply(lambda x: str(x).isspace() or str(x) == "-"))
            gene_description_num = len(set(biomart_pd[~NONE_VIN][0]))

        with open(trinity_stat, "r") as f:
            for line in f:
                if line.startswith("Contig N50"):
                    n50 = line.strip().split(":")[1]
                elif line.startswith("Percent GC"):
                    gc = line.strip().split(":")[1]
                elif line.startswith("Total assembled bases"):
                    size = line.strip().split(":")[1]
                elif line.startswith("Total trinity genes"):
                    chr = line.strip().split(":")[1]

        level_pd = pd.read_table(level_file, header=None)
        level_pd_lev = list(level_pd[1])
        for level_l in [string.lower() for string in level_pd_lev]:
            if level_l == "chromosome":
                chromosome_num += 1
            elif level_l == "scaffold":
                scaffold_num += 1
            elif level_l == "contig":
                contig_num += 1
            else:
                pass

        with open(genome_stat, "r") as fr:
            fr.readline()
            for line in fr:
                tmp = line.strip().split("\t")
                if tmp[4] != "_":
                    protein_coding_num += int(tmp[4])
                if tmp[5] != "_":
                    other_rna_num += int(tmp[5])
                if tmp[6] != "_":
                    pseudogene_num += int(tmp[6])

        main_insert_data = {
            'project_sn': self.bind_object.sheet.project_sn,
            'task_id': self.bind_object.sheet.id,
            'common_name': self.bind_object.sheet.common_name if self.bind_object.sheet.common_name else "",
            'species_name': species_name if species_name else "",
            'ref_anno_version': ref_anno_version if ref_anno_version else "",
            "hyperlink": hyperlink if hyperlink else "",
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'assembly_level': 'Chromosome: ' + str(chromosome_num) + ', Scaffold: ' + str(scaffold_num) + ', Contig: ' + str(contig_num),
            'genome_size': str(round(float(float(size)/1024/1024), 4)) + " Mb",
            'n50': int(n50),
            'protein': protein_coding_num,
            'noncoding': other_rna_num,
            'pseudogene': pseudogene_num,
            'gene': gene_num,
            'transcript': transcript_num,
            'name': gene_name_num,
            'description': gene_description_num,
            'release_date': release_date if release_date else "",
            'chr_num': int(chr),
            'gc': float(gc),
            'version': 'v2',
            'status': 'end',
            'remark': remark,
            'genome_id': genome_id
        }

        if client == "client03":
            self._db = self.db
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            self._db = tsg_client["sanger_ref_genome_db"]
        results = self._db['sg_genome_info'].find({"task_id": self.bind_object.sheet.id})
        for result in results:
            self._db['sg_genome_info'].delete_one(result)

        try:
            main_collection = self._db['sg_genome_info']
            main_collection.insert_one(main_insert_data)
        except Exception, e:
            self.bind_object.set_error("导入参考基因组信息sg_genome_info失败")
        else:
            self.bind_object.logger.info("导入参考基因组信息sg_genome_info成功")