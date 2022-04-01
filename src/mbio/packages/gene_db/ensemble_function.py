# -*- coding: utf-8 -*-
import os
import math
import subprocess
import numpy as np
import pandas as pd
import argparse
from multiprocessing import Pool
import glob
import matplotlib
import gzip
import json
import copy
import pandas as pd
from collections import defaultdict
from collections import OrderedDict
import xml.etree.ElementTree as ET
from mbio.api.database.gene_db import genome
import csv


class EnsembleFunction(object):
    def __init__(self, genome_acc):
        self.table = OrderedDict()
        self.table_show = OrderedDict()
        # self.genome_api = genome.Genome(None)
        # self.chr_dict = self.get_genome_acc(genome_acc)

    def merge(self, mart_sql, conf, data_dir, species_abr, main_table):
        self.get_table_schema(mart_sql, species_abr, data_dir)
        conf_tab = self.get_schema_info(conf, species_abr)
        merged_go_file = main_table + ".go.tsv"
        merged_interpro_file = main_table + ".interpro.tsv"
        merged_domain_file = main_table + ".domain.tsv"
        # self.merge_go_result(conf_tab, data_dir,  main_table, merged_go_file)
        # self.merge_interpro_result(conf_tab, data_dir,  main_table, merged_interpro_file)
        self.merge_domain_result(conf_tab, data_dir,  main_table, merged_domain_file)

    def get_genome_acc(self, genome_acc):
        # print "genome acc is {}".format(genome_acc)
        chr_records = self.genome_api.get_genome_chr(acc_id=genome_acc,
                                                     assembly_unit={"$in": ["Primary Assembly",  "non-nuclear", 'C57BL/6J']})
        '''
        chr_records = self.genome_api.get_genome_chr(acc_id=genome_acc,
                                                     assembly_unit={"$in": ["Primary Assembly", "non-nuclear"]})
        '''
        chr_list = list()
        for record in chr_records:
            # print record
            if record['sequence_role'] == 'assembled-molecule':
                # print record
                if record['sequence_name'].startswith("chr"):
                    record['sequence_name'] = record['sequence_name'].split("chr")[1]
                chr_list.append((record['refseq_accn'], record['sequence_name']))
            elif record['sequence_role'] in ['unlocalized-scaffold', 'unplaced-scaffold']:
                chr_list.append((record['refseq_accn'], record['genbank_accn']))
        return dict(chr_list)


    def get_schema_info(self, conf, species_abr):
        xml_str = None
        with gzip.open(conf, 'r') as conf_f:
            for line in conf_f:
                cols = line.strip().split("\t")
                if cols[1].startswith("<DatasetConfig dataset=\"" + species_abr):
                    xml_str = cols[1].replace("\\n", "\n")
                    break
        if xml_str:
            tree = ET.fromstring(xml_str)
        else:
            raise "conf 文件没有该物种的信息"

        feature_page =  tree.find("AttributePage")
        group_attrs = ['internalName', 'displayName', 'description']
        collection_attrs = ['internalName', 'displayName', 'description']
        desc_attrs = ['internalName', 'displayName', 'description', 'key', 'field', 'linkoutURL', 'tableConstraint', 'useDefault']

        with open(species_abr + ".tab", 'w') as feature_f:
            headers = [a + "_g" for a in group_attrs] + \
                      [a + "_c" for a in collection_attrs] + desc_attrs
            feature_f.write("\t".join(headers) + "\n")

            for attr_group in feature_page.findall("AttributeGroup"):
                group_dict = dict(attr_group.items())
                for attr_collection in attr_group.findall("AttributeCollection"):
                    collection_dict = dict(attr_collection.items())
                    '''
                    if collection_dict.get("internalName", "") in ["ensembl_attributes", "xrefs"]:
                        pass
                    else:
                        continue
                    '''

                    for attr_desc in attr_collection.findall("AttributeDescription"):
                        desc_dict = dict(attr_desc.items())

                        attrs = [group_dict.get(a, "") for a in group_attrs]
                        attrs += [collection_dict.get(a, "") for a in collection_attrs]
                        attrs += [desc_dict.get(a, "") for a in desc_attrs]
                        feature_f.write("\t".join(attrs) + "\n")

        return species_abr + ".tab"


    def get_id_mapping_table(self, conf_tab, intername=None):
        db_dict = dict()
        table = pd.read_table(conf_tab, header=0)

        if type(intername) == str:
            table_xref = table[table['internalName_c'] == intername]
        elif type(intername) == list:
            table_xref = table[table['internalName_c'].isin(intername)]

        for row in table_xref.iterrows():
            row_dict = dict(row[1])
            if row_dict.get("tableConstraint") in db_dict:
                db_dict[row_dict.get("tableConstraint")]["fields"].append(row_dict.get("field"))
                db_dict[row_dict.get("tableConstraint")]["inter_names"].append(row_dict.get("internalName"))
            else:
                db_dict[row_dict.get("tableConstraint")] = {
                    "key": row_dict.get("key"),
                    "inter_names": [row_dict.get("internalName")],
                    "fields": [row_dict.get("field")]
                }
        return db_dict

    def read_db_to_dict(self, db_path, db_schema, get_dict):
        # 提取每一张表的内容
        field_list = get_dict["fields"]
        key = get_dict["key"]

        records_dict = dict()
        # print db_path, db_schema
        with gzip.open(db_path, 'rb') as tf:
            for line in tf:
                table_cols = line.strip("\n").split("\t")
                record = dict(zip(db_schema, table_cols))
                # print record
                if record[key] in records_dict:
                    records_dict[record[key]].append([record[field] for field in field_list])
                    # 如果一个基因对应其它数据库多个属性只取第一个
                    '''
                    add = [record[field] for field in field_list]
                    records_dict[record[key]] = [des + "|" + add[n] for n,des in enumerate(records_dict[record[key]])]
                    '''
                else:
                    records_dict[record[key]] = [[record[field] for field in field_list]]
        return records_dict

    def read_main_table(self, main_table):
        main_id_dict = OrderedDict()
        with open(main_table, 'r') as f:
            for dic in csv.DictReader(f, delimiter='\t'):
                main_id_dict[dic['sg_tran_id']] = {
                    'sg_gene_id': dic['sg_gene_id'],
                    'sg_tran_id': dic['sg_tran_id'],
                    'transcript_id_1064_key': dic['transcript_id_1064_key'],
                    'translation_id_1068_key': dic['translation_id_1068_key'],
                    'ensembl_peptide_id': dic['ensembl_peptide_id']
                }
        return main_id_dict

    def merge_go_result(self, conf_tab, data_dir, main_table, merged_go_file):
        '''
        合并基因id表格
        '''
        go_db_dict = self.get_id_mapping_table(conf_tab, intername="go_closure")
        go_db_list = go_db_dict.keys()

        xref_db_dict = dict()

        main_id_dict = self.read_main_table(main_table)
        for db_name in go_db_list:
            db_schema = self.table[db_name]
            db_path = os.path.join(data_dir, db_name + ".txt.gz")
            get_dict = go_db_dict[db_name]
            db_dict = self.read_db_to_dict(db_path, db_schema, get_dict)

            with open(merged_go_file, 'w') as fo:
                main_list = ['sg_gene_id','sg_tran_id','transcript_id_1064_key','translation_id_1068_key','ensembl_peptide_id']
                fields_list = get_dict['inter_names']
                fo.write("\t".join(main_list + fields_list) + "\n")
                for sg_t, t_dict in main_id_dict.items():
                    if t_dict['transcript_id_1064_key'] in db_dict:
                        main_value_list = [t_dict.get(m, "") for m in main_list]
                        for gos in db_dict[t_dict['transcript_id_1064_key']]:
                            fo.write("\t".join(main_value_list + gos) + "\n")

    def merge_interpro_result(self, conf_tab, data_dir, main_table, merged_interpro_file):
        '''
        合并基因id表格
        '''
        interpro_db_dict = self.get_id_mapping_table(conf_tab, intername="prot_interpro")
        interpro_db_list = interpro_db_dict.keys()

        xref_db_dict = dict()

        main_id_dict = self.read_main_table(main_table)
        for db_name in interpro_db_list:
            db_schema = self.table[db_name]
            db_path = os.path.join(data_dir, db_name + ".txt.gz")
            get_dict = interpro_db_dict[db_name]
            db_dict = self.read_db_to_dict(db_path, db_schema, get_dict)

            with open(merged_interpro_file, 'w') as fo:
                main_list = ['sg_gene_id','sg_tran_id','transcript_id_1064_key','translation_id_1068_key','ensembl_peptide_id']
                fields_list = get_dict['inter_names']
                fo.write("\t".join(main_list + fields_list) + "\n")
                for sg_t, t_dict in main_id_dict.items():
                    if t_dict['transcript_id_1064_key'] in db_dict:
                        main_value_list = [t_dict.get(m, "") for m in main_list]
                        for interpros in db_dict[t_dict['transcript_id_1064_key']]:
                            fo.write("\t".join(main_value_list + interpros) + "\n")


    def merge_domain_result(self, conf_tab, data_dir, main_table, merged_domain_file):
        '''
        合并基因id表格
        '''
        domain_db_dict = self.get_id_mapping_table(conf_tab, intername=['protein_domain', 'protein_feature'])
        domain_db_list = domain_db_dict.keys()

        xref_db_dict = dict()

        main_id_dict = self.read_main_table(main_table)
        main_list = ['sg_gene_id','sg_tran_id','transcript_id_1064_key','translation_id_1068_key','ensembl_peptide_id']
        fields_list = ['domain_id', 'domain_start', 'domain_end', 'domain_db']
        # fields_list = get_dict['fields']
        
        with open(merged_domain_file, 'w') as fo:
            fo.write("\t".join(main_list + fields_list) + "\n")
        for db_name in domain_db_list:
            print db_name
            db_schema = self.table[db_name]
            db_path = os.path.join(data_dir, db_name + ".txt.gz")
            get_dict = domain_db_dict[db_name]
            db_dict = self.read_db_to_dict(db_path, db_schema, get_dict)

            with open(merged_domain_file, 'aw') as fo:
                # main_list = ['sg_gene_id','sg_tran_id','transcript_id_1064_key','translation_id_1068_key','ensembl_peptide_id']
                fields_list = get_dict['inter_names']
                print fields_list
                # fo.write("\t".join(main_list + fields_list) + "\n")
                for sg_t, t_dict in main_id_dict.items():
                    main_value_list = [t_dict.get(m, "") for m in main_list]
                    for domains in db_dict[t_dict['translation_id_1068_key']]:
                        fo.write("\t".join(main_value_list + domains + [fields_list[0]]) + "\n")


    def merge_external_tran_id(self, conf_tab):
        pass


    def get_table_schema(self, mart_sql, species_abr, data_dir):
        with gzip.open(mart_sql, 'rb') as mart_sql_f:
            start = False
            table_name = None
            table_cols = None
            table_cols2 = None

            for line in mart_sql_f:
                # print "one line"

                # print line
                if line.startswith(")"):
                    start = False

                if start == True:
                    field = line.strip().split()[0]
                    # print field
                    if field.startswith("`"):
                        if table_name in self.table:
                            self.table[table_name].append(field[1:-1])

                        else:
                            self.table[table_name] = [field[1:-1]]

                    else:
                        pass

                if line.startswith("CREATE TABLE `" + species_abr):
                    if "_homolog_" in line:
                        pass
                    else:
                        start = True
                        table_name = line.split()[2][1:-1]
                        '''
                        with gzip.open(os.path.join(data_dir, table_name + ".txt.gz"), 'rb') as tf:
                            line = tf.readline()
                            table_cols = line.strip().split("\t")
                            table_cols.reverse()
                        '''
        with open(species_abr + 'table.json', 'w') as f:
            f.write(json.dumps(self.table, indent=4))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-acc', type=str, required=True,
                        help='genome acc')
    parser.add_argument('-mart_sql', type=str, required=True,
                        help='sql database file create for biomart')
    parser.add_argument('-conf', type=str, required=True,
                        help='conf for biomart')
    parser.add_argument('-data_dir', type=str, default=True,
                        help="database directory")
    parser.add_argument('-species_abr', type=str, default=True,
                        help="species abr string")
    parser.add_argument('-main_table', type=str, default=True,
                        help="merged id mapping file")

    # ----------------------------------------------------------------------------------------------
    args = parser.parse_args()
    ef = EnsembleFunction(args.acc)
    ef.merge(mart_sql=args.mart_sql, conf=args.conf, data_dir=args.data_dir, species_abr=args.species_abr, main_table=args.main_table)
    # ef_j = json(ef.table)
    # print json.dumps(ef.table, indent=4)


