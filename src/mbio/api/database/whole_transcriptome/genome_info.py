# -*- coding: utf-8 -*-
# __author__ = 'shijin'
# last modified by shicaiping at 20180510

from biocluster.api.database.base import Base, report_check
import os
import datetime
import unittest
from biocluster.config import Config
from mbio.api.database.whole_transcriptome.api_base import ApiBase

class GenomeInfo(ApiBase):
    def __init__(self, bind_object=None):
        super(GenomeInfo, self).__init__(bind_object)
        self._project_type = 'whole_transcriptome'

    @report_check
    def add_genome_info(self, file_path, species_name, ref_anno_version, hyperlink, major = True):
        main_insert_data = {
            'project_sn': self.bind_object.sheet.project_sn,
            'task_id': self.bind_object.sheet.id,
            'species_name': species_name if species_name else "",
            'ref_anno_version': ref_anno_version if ref_anno_version else "",
            "hyperlink": hyperlink if hyperlink else "",
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'version': 'v1'
        }
        main_collection = self.db['species_information']
        main_id = main_collection.insert_one(main_insert_data).inserted_id
        if major:
            self.add_genome_info_detail(file_path=file_path, inserted_id=main_id)
            self.update_db_record('species_information', main_id, status="end", main_id=main_id)

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
                    "chromosome":chr,
                    "length":size,
                    "qc_percent": gc,
                    "gene": gene,
                    "proteincoding": protein_coding,
                    "other_rna": other_rna,
                    "pseudogene": pseudo
                }
                insert_data_list.append(insert_data)
        collections = self.db['species_information_detail']
        collections.insert_many(insert_data_list)
        self.bind_object.logger.info("导入species_information_detail细节表成功")

class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
      from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
      from biocluster.wsheet import Sheet
      import random

      data = {
            "id": "whole_transcriptome",
            "project_sn": "whole_transcriptome",
            "type": "workflow",
            "name": "whole_transcriptome_test_api",
            "options": {
            },
        }
      wsheet = Sheet(data=data)
      wf = RefrnaTestApiWorkflow(wsheet)
      wf.IMPORT_REPORT_DATA = True
      wf.IMPORT_REPORT_AFTER_END = False
      wf.test_api = wf.api.api("whole_transcriptome.genome_info")
      genome_info = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf.genome_stat.xls"
      wf.test_api.add_genome_info(file_path=genome_info, species_name="Homo_sapiens", ref_anno_version="GRCh38", hyperlink="http://asia.ensembl.org/Homo_sapiens/Info/Index", major = True)

if __name__ == '__main__':
    unittest.main()
