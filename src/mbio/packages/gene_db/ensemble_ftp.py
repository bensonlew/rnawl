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


class EnsembleFtp(object):
    def __init__(self, genome_acc):
        self.table = OrderedDict()
        self.table_show = OrderedDict()
        self.genome_api = genome.Genome(None)
        self.chr_dict = self.get_genome_acc(genome_acc)

    def merge(self, mart_sql, conf, data_dir, species_abr):
        self.get_table_schema(mart_sql, species_abr, data_dir)
        conf_tab = self.get_schema_info(conf, species_abr)
        # self.chr_list = self.get_genome_acc(self)
        merged_gene_file = species_abr + "merged_gene_table.tsv"
        main_table = species_abr + "__gene__main"
        main_table_dict = OrderedDict((
            ('gene_id_1020_key', 'gene_id_1020_key'),
            ('biotype_1020', 'gene_biotype'),
            ('stable_id_1023', 'ensembl_gene_id'),
            ('display_label_1074',  'external_gene_name'),
            ('name_1059', 'chromosome_name'),
            ('seq_region_start_1020', 'start_position'),
            ('seq_region_end_1020', 'end_position'),
            ('seq_region_strand_1020', 'strand'),
            ('gene__main_stable_id_version', 'ensembl_gene_id_version')
        ))

        
        main_field_retain = [k for k in main_table_dict if k in self.table[main_table]]
        main_field_internames = [v for k,v in main_table_dict.items() if k in self.table[main_table]]
        
        self.merge_external_gene_id(conf_tab, data_dir, merged_gene_file, main_table, main_field_retain, main_field_internames, level="G")


        main_table_dict.update(OrderedDict((
            ('transcript_id_1064_key', 'transcript_id_1064_key'),
            ('biotype_1064', 'transcript_biotype'),
            ('stable_id_1066', 'ensembl_transcript_id'),
            ('display_label_1074_r1', 'external_transcript_name'),
            ('transcript__main_stable_id_version', 'ensembl_transcript_id_version'),
            ('translation_id_1068_key', 'translation_id_1068_key'),
            ('stable_id_1070', 'ensembl_peptide_id'),
            ('translation__main_stable_id_version', 'ensembl_peptide_id_version')
        )))

        main_table = species_abr + "__translation__main"
        main_field_retain = [k for k in main_table_dict if k in self.table[main_table]]
        main_field_internames = [v for k,v in main_table_dict.items() if k in self.table[main_table]]
        merged_gene_file = species_abr + "merged_tran_table.tsv"
        exon_table = species_abr + "__exon_transcript__dm"
        self.merge_external_gene_id(conf_tab, data_dir, merged_gene_file, main_table, main_field_retain, main_field_internames, exon_table=exon_table)

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


    def get_id_mapping_table(self, conf_tab):
        db_dict = dict()
        table = pd.read_table(conf_tab, header=0)

        table_xref = table[table['internalName_c'] == "xrefs"]

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
                    pass
                    # 如果一个基因对应其它数据库多个属性只取第一个
                    '''
                    add = [record[field] for field in field_list]
                    records_dict[record[key]] = [des + "|" + add[n] for n,des in enumerate(records_dict[record[key]])]
                    '''
                else:
                    records_dict[record[key]] = [record[field] for field in field_list]
        return records_dict

    def get_trans_exon(self, data_dir, exon_table):
        tran_exon =dict()
        exon_file = os.path.join(data_dir, exon_table + ".txt.gz")
        exon_schema = self.table[exon_table]
        with gzip.open(exon_file, 'rb') as fin:
            for line in fin:
                cols = line.strip().split("\t")
                record = dict(zip(exon_schema, cols))
                region = record['seq_region_start_1015'] + "." + record['seq_region_end_1015']
                if record["transcript_id_1064_key"] in tran_exon:
                    tran_exon[record["transcript_id_1064_key"]].append(region)
                else:
                    tran_exon[record["transcript_id_1064_key"]] = [region]
        return tran_exon
    def merge_external_gene_id(self, conf_tab, data_dir, merged_gene_file, main_table, main_field_retain, main_field_internames, level="T", exon_table=None):
        '''
        合并基因id表格
        '''
        db_dict = self.get_id_mapping_table(conf_tab)
        if level == "G":
            gene_db_list = [db_name for db_name in db_dict.keys() if db_dict[db_name]['key'].startswith("gene_id_")]
        else:
            gene_db_list = db_dict.keys()

        xref_db_dict = dict()
        for db_name in gene_db_list:
            db_schema = self.table[db_name]
            db_path = os.path.join(data_dir, db_name + ".txt.gz")
            get_dict = db_dict[db_name]

            xref_db_dict[db_name] = self.read_db_to_dict(db_path, db_schema, get_dict)

        main_file = os.path.join(data_dir, main_table + ".txt.gz")
        main_schema = self.table[main_table]
        inter2field = dict(zip(main_field_internames, main_field_retain))
        with gzip.open(main_file, 'rb') as fin, open(merged_gene_file, 'w') as fo:
            if level == "G":
                fo.write("sg_gene_id" + "\t" +  "\t".join(main_field_internames))
            else:
                tran_exon = self.get_trans_exon(data_dir, exon_table)
                fo.write("sg_gene_id\tsg_tran_id" + "\t" +  "\t".join(main_field_internames))

            for db_name in gene_db_list:
                fo.write("\t" + "\t".join(db_dict[db_name]["inter_names"]))
            fo.write("\n")
            for line in fin:
                cols = line.strip().split("\t")
                record = dict(zip(main_schema, cols))
                if record[inter2field['strand']] == "1":
                    record[inter2field['strand']] = "+"
                elif record[inter2field['strand']] == "-1":
                    record[inter2field['strand']] = "-"

                # print "chr_dict is {}".format(self.chr_dict)
                if record[inter2field['chromosome_name']] not in self.chr_dict.values():
                    continue
                strand_id = {"+": "1", "-": "0"}
                std_id = record[inter2field['chromosome_name']] + '_' + record[inter2field['start_position']] + '_' + record[inter2field['end_position']] + '_' + strand_id[record[inter2field['strand']]]
                if level == "G":
                    fo.write(std_id + "\t" + "\t".join([record[field] for field in main_field_retain]))
                else:
                    std_tran_id = "_".join([record[inter2field['chromosome_name']]]
                                           + tran_exon[record['transcript_id_1064_key']]
                                           + [strand_id[record[inter2field['strand']]]])
                    fo.write(std_id + "\t" + std_tran_id + "\t" + "\t".join([record[field] for field in main_field_retain]))

                for db_name in gene_db_list:
                    get_dict = db_dict[db_name]
                    key_field = record[get_dict["key"]]
                    # print xref_db_dict[db_name]
                    if key_field in xref_db_dict[db_name]:
                        fo.write("\t" + "\t".join(xref_db_dict[db_name][key_field]))
                    else:
                        fo.write("\t" + "\t".join(len(db_dict[db_name]["inter_names"]) * ['\N']))
                fo.write("\n")


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
    # ----------------------------------------------------------------------------------------------
    args = parser.parse_args()
    ef = EnsembleFtp(args.acc)
    ef.merge(mart_sql=args.mart_sql, conf=args.conf, data_dir=args.data_dir, species_abr=args.species_abr)
    # ef_j = json(ef.table)
    # print json.dumps(ef.table, indent=4)


