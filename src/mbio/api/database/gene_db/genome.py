# !/usr/bin/python
# -*- coding: utf-8 -*-
import types
import re
import os
import json
import math
from collections import OrderedDict
import unittest
import datetime
import glob
from bson.objectid import ObjectId
from api_base import ApiBase
import sys
from mbio.packages.gene_db.accession2taxon import taxon


class Genome(ApiBase):
    def __init__(self, bind_object):
        self.taxon = taxon()
        super(Genome, self).__init__(bind_object)

    def add_genome_info(self, main_info={}, genome_report=None):
        """
        add transcript and gene mapping file to database
        :param t2g_file: mapping file with two columns: transcript gene. No header line.
        :param project_sn: project number
        :param task_id: task id
        :return: main table id.
        """

        with open(genome_report, 'r') as f:
            for line in f:
                if line.strip() == "#":
                    break
                else:
                    cols = line.strip().split(":")
                    main_info.update({
                        re.sub(r'[^\da-zA-Z_]', "_", cols[0].lstrip("# ")).lower() : cols[1].strip()
                    })
        '''
        main_info = {
            'assembly_level': 'Chromosome',
            'assembly_name': 'GRCh38.p13',
            'assembly_type': 'haploid-with-alt-loci',
            'bioproject': 'PRJNA31257',
            'date': '2019-02-28',
            'description': 'Genome Reference Consortium Human Build 38 patch release 13 (GRCh38.p13)',
            'genbank_assembly_accession': 'GCA_000001405.28',
            'genome_representation': 'full',
            'organism_name': 'Homo sapiens (human)',
            'refseq_category': 'Reference Genome',
            'release_type': 'patch',
            'submitter': 'Genome Reference Consortium',
            'taxid': '9606'
        }
        '''

        time_now = datetime.datetime.now()
        main_info.update({
            "created_ts" : time_now.strftime('%Y-%m-%d %H:%M:%S'),
            "status" : 'start'
        })

        main_id = self.create_db_table('sgdb_genome', [main_info])
        self.add_genome_chr(main_info.get("genbank_assembly_accession", "unknown") , genome_report)
        self.update_db_record('sgdb_genome', main_id, status="end", main_id=main_id)
        return main_id

    def get_genome_info(self, species):
        genome_dict = self.db['sgdb_genome'].find_one({"species": species})
        return genome_dict

    def get_genome_chr(self, **kwargs):
        """
        :acc_id acession_id
        """
        results = self.get_tables_by_main_record("sgdb_genome_chr", **kwargs)
        return results


    def add_genome_chr(self, acc_id, genome_report):
        """
        :version:version info
        """

        headers = []
        chrs = []
        with open(genome_report, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    headers.append(line.strip())
                else:
                    chrs.append(line.strip())

        name_old = headers[-1].lstrip("# ").split("\t")
        name_std = [re.sub(r'[^\da-zA-Z_]', "_", name).lower() for name in name_old]

        data_list = []
        for chrline in chrs:
            chr_list = chrline.strip().split("\t")
            chr_dict = dict(zip(name_std, chr_list))
            chr_dict.update({
                "acc_id": acc_id
            })
            data_list.append(chr_dict)

        self.create_db_table('sgdb_genome_chr', data_list)
        return True

    def update_taxon(self):
        col = self.db['sgdb_genome']
        genome_records = col.find()
        for genome_record in genome_records:
            if 'taxid' in genome_record:
                linkage = self.taxon.get_linkage_name(genome_record['taxid'])
                self.update_db_record('sgdb_genome', genome_record["_id"], linkage=linkage)



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def __init__(self, method_name, report_file):
        self.report_file = report_file
        super(TestFunction, self).__init__(methodName=method_name)
        
    toolbox = Genome(None)

    def add_genome_info(self):

        main_info = dict(
        )
        # test_dir = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/gene_db"
        self.toolbox.add_genome_info(main_info, self.report_file)

    def get_genome_info(self):
        records = self.toolbox.get_genome_chr(acc_id=ObjectId("5e65e83c17b2bf70c6338aff"), assembly_unit={"$in": ["Primary Assembly", "non-nuclear"]})
        for record in records:
            print record

    def update_taxon(self):
        self.toolbox.update_taxon()



if __name__ == '__main__':
    suite = unittest.TestSuite()
    '''
    if sys.argv[1] in ["-h", "-help", "--h", "--help"]:
        print "\n".join(["add_genome_info", "get_genome_info"])
        if len(sys.argv) == 3:
    '''
    run_fun = sys.argv[1]
    
    # report_file = sys.argv[1]
    suite.addTest(TestFunction(run_fun, None))
    unittest.TextTestRunner(verbosity=2).run(suite)


