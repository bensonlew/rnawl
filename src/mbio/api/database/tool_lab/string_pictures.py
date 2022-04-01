# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20210202
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
import glob
import re
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class StringPictures(ApiBase):
    def __init__(self, bind_object):
        super(StringPictures, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_string(self, main_id, string_dir, s3_path=None, params=None, project_sn='tool_lab', task_id='tool_lab',):
        # add main table info
        if main_id is None:
            name = "String_Picture" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='String picture',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_string_pictures', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        annotation_file_path = glob.glob(os.path.join(string_dir, '*annotation.xls'))[0]
        annotation_file = pd.read_table(annotation_file_path, header=0, index_col=None, sep='\t')
        annotation_file.rename(columns={'#node': 'node'}, inplace=True)
        annotation_columns_list = list()
        annotation_columns_list.append({'field': 'node', 'filter': False, 'sort': False, 'title': 'Node', 'type': 'string'})
        annotation_columns_list.append({'field': 'accession', 'filter': False, 'sort': False, 'title': 'Node Accession ID', 'type': 'string'})
        annotation_columns_list.append({'field': 'string_id', 'filter': False, 'sort': False, 'title': 'String ID ', 'type': 'string'})
        annotation_columns_list.append({'field': 'annotation', 'filter': False, 'sort': False, 'title': 'Annotation ', 'type': 'string'})
        annotation_data_columns = {'column': annotation_columns_list, 'condition': {}}
        annotation_columns_data = json.dumps(annotation_data_columns)
        annotation_file['string_picture_id'] = main_id
        annotation_detail = annotation_file.to_dict('r')
        self.create_db_table('sg_string_picture_annotation_detail', annotation_detail)
        self.update_db_record('sg_string_picture', main_id, main_id=main_id, annotation_column_data_detail=annotation_columns_data)

        bitscore_file_path = glob.glob(os.path.join(string_dir, '*bitscore.xls'))[0]
        bitscore_file = pd.read_table(bitscore_file_path, header=0, index_col=None, sep='\t')
        bitscore_columns_list = list()
        bitscore_columns_list.append({'field': 'ncbiTaxonId_A', 'filter': False, 'sort': False, 'title': 'NCBI Taxon ID (A)', 'type': 'string'})
        bitscore_columns_list.append({'field': 'stringId_A', 'filter': False, 'sort': False, 'title': 'String ID (A)', 'type': 'string'})
        bitscore_columns_list.append({'field': 'accession_A', 'filter': False, 'sort': False, 'title': 'Accession (A)', 'type': 'string'})
        bitscore_columns_list.append({'field': 'ncbiTaxonId_B', 'filter': False, 'sort': False, 'title': 'NCBI Taxon ID (B)', 'type': 'string'})
        bitscore_columns_list.append({'field': 'stringId_B', 'filter': False, 'sort': False, 'title': 'String ID (B)', 'type': 'string'})
        bitscore_columns_list.append({'field': 'accession_B', 'filter': False, 'sort': False, 'title': 'Accession (B)', 'type': 'string'})
        bitscore_columns_list.append({'field': 'bitscore', 'filter': False, 'sort': False, 'title': 'Bitscore', 'type': 'float'})
        bitscore_data_columns = {'column': bitscore_columns_list, 'condition': {}}
        bitscore_columns_data = json.dumps(bitscore_data_columns)
        bitscore_file['string_picture_id'] = main_id
        bitscore_detail = bitscore_file.to_dict('r')
        self.create_db_table('sg_string_picture_bitscore_detail', bitscore_detail)
        self.update_db_record('sg_string_picture', main_id, main_id=main_id, bitscore_column_data_detail=bitscore_columns_data)

        interaction_file_path = glob.glob(os.path.join(string_dir, '*interaction.xls'))[0]
        interaction_file = pd.read_table(interaction_file_path, header=0, index_col=None, sep='\t')

        columns_dict = {'Accession (protein A)': 'accession(protein_A)', 'Accession (protein B)': 'accession(protein_B)',
                        'STRING identifier (protein A)': 'string_identifier(protein_A)',
                        'STRING identifier (protein B)': 'string_identifier(protein_B)',
                        'common protein name (protein A)': 'common_protein_name(protein_A)',
                        'common protein name (protein B)': 'common_protein_name(protein_B)',
                        'NCBI taxon identifier': 'ncbi_taxon_identifier',
                        'combined score': 'combined_score',
                        'gene neighborhood score': 'gene_neighborhood_score',
                        'gene fusion score': 'gene_fusion_score',
                        'phylogenetic profile score': 'phylogenetic_profile_score',
                        'coexpression score': 'coexpression_score',
                        'experimental score': 'experimental_score',
                        'database score': 'database_score',
                        'textmining score': 'textmining_score'}
        interaction_file.rename(columns=columns_dict, inplace=True)
        interaction_file.fillna('', inplace=True)
        interaction_columns_list = list()
        interaction_columns_list.append({'field': 'accession(protein_A)', 'filter': False, 'sort': False, 'title': 'Protein ID(A)', 'type': 'string'})
        interaction_columns_list.append({'field': 'accession(protein_B)', 'filter': False, 'sort': False, 'title': 'Protein ID (B)', 'type': 'string'})
        interaction_columns_list.append({'field': 'string_identifier(protein_A)', 'filter': False, 'sort': False, 'title': 'String ID (A)', 'type': 'string'})
        interaction_columns_list.append({'field': 'string_identifier(protein_B)', 'filter': False, 'sort': False, 'title': 'String ID (B)', 'type': 'string'})
        interaction_columns_list.append({'field': 'common_protein_name(protein_A)', 'filter': False, 'sort': False, 'title': 'Common protein name (A)', 'type': 'string'})
        interaction_columns_list.append({'field': 'common_protein_name(protein_B)', 'filter': False, 'sort': False, 'title': 'Common protein name (B)', 'type': 'string'})
        interaction_columns_list.append({'field': 'ncbi_taxon_identifier', 'filter': False, 'sort': False, 'title': 'Ncbi Taxon Identifier', 'type': 'string'})
        interaction_columns_list.append({'field': 'combined_score', 'filter': False, 'sort': False, 'title': 'Combined Score', 'type': 'float'})
        interaction_columns_list.append({'field': 'gene_neighborhood_score', 'filter': False, 'sort': False, 'title': 'Gene Neighborhood Score', 'type': 'float'})
        interaction_columns_list.append({'field': 'gene_fusion_score', 'filter': False, 'sort': False, 'title': 'Gene Fusion Score', 'type': 'float'})
        interaction_columns_list.append({'field': 'phylogenetic_profile_score', 'filter': False, 'sort': False, 'title': 'Phylogenetic Profile Score', 'type': 'float'})
        interaction_columns_list.append({'field': 'coexpression_score', 'filter': False, 'sort': False, 'title': 'Coexpression Score', 'type': 'float'})
        interaction_columns_list.append({'field': 'experimental_score', 'filter': False, 'sort': False, 'title': 'Experimental Score', 'type': 'float'})
        interaction_columns_list.append({'field': 'database_score', 'filter': False, 'sort': False, 'title': 'Database Score', 'type': 'float'})
        interaction_columns_list.append({'field': 'textmining_score', 'filter': False, 'sort': False, 'title': 'Textmining Score', 'type': 'float'})

        interaction_data_columns = {'column': interaction_columns_list, 'condition': {}}
        interaction_columns_data = json.dumps(interaction_data_columns)
        interaction_file['string_picture_id'] = main_id
        interaction_detail = interaction_file.to_dict('r')
        self.create_db_table('sg_string_picture_interaction_detail', interaction_detail)
        self.update_db_record('sg_string_picture', main_id, main_id=main_id, interaction_column_data_detail=interaction_columns_data)
        self.update_db_record('sg_string_picture', main_id, status='end', main_id=main_id, pictures_file=s3_path)



