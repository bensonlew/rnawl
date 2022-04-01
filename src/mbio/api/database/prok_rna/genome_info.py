# -*- coding: utf-8 -*-
# __author__ = 'shijin'
# last modified by shicaiping at 20180510

from biocluster.api.database.base import Base, report_check
import os
import datetime
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class GenomeInfo(ApiBase):
    def __init__(self, bind_object=None):
        super(GenomeInfo, self).__init__(bind_object)
        self._project_type = 'prok_rna'

    @report_check
    def add_genome_info(self, file_path, species_name, accession, genome_size, gc_content, major = True):
        main_insert_data = {
            'project_sn': self.bind_object.sheet.project_sn,
            'task_id': self.bind_object.sheet.id,
            'species_name': species_name if species_name else "",
            'accession': accession if accession else "",
            'genome_size': genome_size,
            "gc_content": gc_content,
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        main_collection = self.db['sg_species_information']
        main_id = main_collection.insert_one(main_insert_data).inserted_id
        self.update_db_record('sg_species_information', main_id, status="end", main_id=main_id)
        if major:
            self.add_genome_info_detail(file_path=file_path, inserted_id=main_id)


    @report_check
    def add_genome_info_detail(self, file_path, inserted_id):
        insert_data_list = []
        with open(file_path, "r") as fr:
            fr.readline()
            for line in fr:
                tmp = line.strip().split("\t")
                chr = tmp[0]
                size = tmp[1]
                gc = tmp[2]
                gene = tmp[3]
                protein_coding = tmp[4]
                other_rna = tmp[5]
                pseudo = tmp[6]
                insert_data = {
                    "species_id": inserted_id,
                    "ref_seq":chr,
                    "length":size,
                    "qc_percent": gc,
                    "gene": gene,
                    "proteincoding": protein_coding,
                    "other_rna": other_rna,
                    "pseudogene": pseudo
                }
                insert_data_list.append(insert_data)
        collections = self.db['sg_species_information_detail']
        collections.insert_many(insert_data_list)
        self.bind_object.logger.info("导入sg_species_information_detail细节表成功")
