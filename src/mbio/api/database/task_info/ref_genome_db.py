# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import datetime
from biocluster.api.database.base import Base, report_check
from bson import SON
import pymongo


class RefGenomeDb(Base):
    def __init__(self, bind_object):
        super(RefGenomeDb, self).__init__(bind_object)
        self._db = None
        self._project_type = 'ref_genome_db'

    #@report_check
    def add_task_info(self, client, cluster, genome_id=None, result_dir=None, raw_input=None, db_name=None, annot_version_dict=None, medical=None):
        if client == "client03":
            if db_name:
                self._db = db_name
            else:
                self._db = self.db
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            self._db = tsg_client["sanger_ref_genome_db"]
        if client == "client03":
            status = 'tsg_finish'
        elif cluster == "sanger":
            status = 'sanger_finish'
        else:
            status = 'isanger_finish'

        database_list = [
            {
                "software_database" : "NR",
                "source" : "ftp://ftp.ncbi.nlm.nih.gov/blast/db/",
                "version" : "Version {}".format(annot_version_dict.get("nr", "2019.6.26"))
            },
            {
                "software_database" : "GO",
                "source" : "http://www.geneontology.org/",
                "version" : "Version {}".format(annot_version_dict.get("go", "2019.7.1"))
            },
            {
                "software_database" : "KEGG",
                "source" : "http://www.genome.jp/kegg/",
                "version" : "Version {}".format(annot_version_dict.get("kegg", "2020.03"))
            },
            {
                "software_database" : "Rfam",
                "source" : "http://rfam.janelia.org/",
                "version" : "Version Rfam v{}".format(annot_version_dict.get("rfam", "14.1"))
            },
            {
                "software_database" : "Swiss-Prot",
                "source" : "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
                "version" : "Version {}".format(annot_version_dict.get("swissprot", "2019.7.1"))
            },
            {
                "software_database" : "eggNOG",
                "source" : "http://eggnogdb.embl.de/#/app/home",
                "version" : "Version {}".format(annot_version_dict.get("eggnog", "5.0"))
            }
        ]

        if medical:
            database_list += [
                {
                    "software_database" : "Uniprot",
                    "source" : "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete",
                    "version" : "Version {}".format(annot_version_dict.get("uniprot", "2020.9.1"))
                },
                {
                    "software_database" : "Reactome",
                    "source" : "https://reactome.org/",
                    "version" : "Version {}".format(annot_version_dict.get("reactome", "72"))
                },
                {
                    "software_database" : "DO",
                    "source" : "https://disease-ontology.org/",
                    "version" : "Version {}".format(annot_version_dict.get("do", "2020.9.3"))
                },
                {
                    "software_database" : "Disgenet",
                    "source" : "https://www.disgenet.org/",
                    "version" : "Version {}".format(annot_version_dict.get("disgenet", ""))
                }
            ]
        json_data = [
            ('task_id', self.bind_object.sheet.id),
            ('member_id', self.bind_object.sheet.member_id),
            ('project_sn', self.bind_object.sheet.project_sn),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('is_demo', 0),
            ('demo_id', self.bind_object.sheet.id),
            ('member_type', int(self.bind_object.sheet.member_type) if self.bind_object.sheet.member_type else self.bind_object.sheet.member_type),
            ('cmd_id', int(self.bind_object.sheet.cmd_id) if self.bind_object.sheet.cmd_id else self.bind_object.sheet.cmd_id),
            ('raw_input', raw_input),
            ('genome_id', genome_id),
            ('result_dir', result_dir),
            ('version', 'v2'),
            ('database_list', database_list),
            ('task_status', status)  # 用于记录任务的状态，分别为tsg_finish, sanger_finish, isanger_finish
        ]

        results = self._db['sg_task'].find({"task_id": self.bind_object.sheet.id})
        for result in results:
            self._db['sg_task'].delete_one(result)
        self._db['sg_task'].insert_one(SON(json_data))
        self.bind_object.logger.info('任务信息导入sg_task成功。')

    def find_task_info(self):
        client = self.bind_object.sheet.client
        if client == "client03":
            self._db = self.db
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            self._db = tsg_client["sanger_ref_genome_db"]
        results = self._db['sg_task'].find({"task_id": self.bind_object.sheet.id})
        genome_ids = list()
        if results.count() != 0:
            for result in results:
                if 'genome_id' in result and result['genome_id'] not in genome_ids:
                    genome_ids.append(result['genome_id'])
        return genome_ids