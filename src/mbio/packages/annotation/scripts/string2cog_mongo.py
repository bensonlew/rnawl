# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 2016.11.18

import sys
from pymongo import MongoClient
from Bio.Blast import NCBIXML
import xml.etree.ElementTree as ET


class string2cog(object):
    def __init__(self):
        self.client = MongoClient('mongodb://10.100.200.129:27017')
        self.cog_string = self.client.sanger_biodb.COG_String
        self.cog = self.client.sanger_biodb.COG
        self.func_type = {
            'INFORMATION STORAGE AND PROCESSING': sorted(['J', 'A', 'K', 'L', 'B']),
            'CELLULAR PROCESSES AND SIGNALING': sorted(['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O']),
            'METABOLISM': sorted(['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q']),
            'POORLY CHARACTERIZED': sorted(['R', 'S']),
        }
        self.func_decs = {
            'J': 'Translation, ribosomal structure and biogenesis',
            'A': 'RNA processing and modification', 'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'B': 'Chromatin structure and dynamics',
            'D': 'Cell cycle control, cell division, chromosome partitioning',
            'Y': 'Nuclear structure', 'V': 'Defense mechanisms', 'T': 'Signal transduction mechanisms',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility', 'Z': 'Cytoskeleton', 'W': 'Extracellular structures',
            'U': 'Intracellular trafficking, secretion, and vesicular transport',
            'O': 'Posttranslational modification, protein turnover, chaperones',
            'C': 'Energy production and conversion', 'G': 'Carbohydrate transport and metabolism',
            'E': 'Amino acid transport and metabolism', 'F': 'Nucleotide transport and metabolism',
            'H': 'Coenzyme transport and metabolism', 'I': 'Lipid transport and metabolism',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
            'R': 'General function prediction only', 'S': 'Function unknown'
        }

    def string2cog_by_mongo(self, string_xml, out_dir):
        cog_list_path = out_dir + '/cog_list.xls'
        cog_table_path = out_dir + '/cog_table.xls'
        cog_summary_path = out_dir + '/cog_summary.xls'
        with open(cog_list_path, 'wb') as listfile, open(cog_table_path, 'wb') as tablefile, open(cog_summary_path, 'wb') as summaryfile:
            listfile.write('Query_name\tCOG\tNOG\n')
            tablefile.write('#Query_name\tQuery_length\tHsp_start_of_query\tHsp_end_of_query\tHsp_strand_of_query\tHit_name\tHit_description\tHit_length\tHsp_start_of_hit\tHsp_end_of_hit\tCOG/NOG_group\tCOG/NOG_group_description\tCOG/NOG_group_categories\tCOG/NOG_region_start\tCOG/NOG_region_end\tCoverage_of_COG/NOG_region\tIdentities_of_COG/NOG_region\tPositives_Identities_of_COG/NOG_region\n')
            c = self.cog_string
            totalseq = 0  # ע��COG/NOG����������
            funcount = {'COG': {}, 'NOG': {}}
            document = ET.parse(string_xml)
            root = document.getroot()
            identations = root.find('BlastOutput_iterations')
            for identation in identations.findall('Iteration'):
                querydef = identation.find('Iteration_query-def')
                queryname = querydef.text.split(' ')[0]  # the query name
                query_length = identation.find('Iteration_query-len').text  # the query length
                iter_hits = identation.find('Iteration_hits')
                hits = iter_hits.findall('Hit')
                if len(hits) > 0:  # check if such name has hits
                    cog_list = []
                    duplicate = False
                    for hit in hits:
                        hit_id = hit.find('Hit_id').text
                        hit_length = hit.find('Hit_len').text
                        OG = c.find({'string_id': hit_id})
                        if OG.count() >= 1:
                            if duplicate is False:
                                totalseq += 1
                                duplicate = True
                            for og in OG:
                                ortho_group = og['orthologous_group']  # COG/NOG group
                                if len(og['function_categories']) > 1:
                                    func_cat = og['function_categories']  # modify zhouxuan 20170606
                                    # func_cat = ','.split(og['function_categories'])
                                else:
                                    func_cat = [og['function_categories']]
                                try:
                                    hit_desc = og['protein_annotation']
                                except KeyError:
                                    hit_desc = ''
                                region_start = og['start_position']
                                region_end = og['start_position']
                                if ortho_group not in cog_list:
                                    if ortho_group.startswith('COG'):
                                        # deal with cog_summary.xls###
                                        for item in func_cat:
                                            if item in funcount['COG']:
                                                funcount['COG'][item] += 1
                                            else:
                                                funcount['COG'][item] = 1
                                        listfile.write('{}\t{}\t\n'.format(queryname, ortho_group))
                                        cog_list.append(ortho_group)
                                    else:
                                        for item in func_cat:
                                            if item in funcount['NOG']:
                                                funcount['NOG'][item] += 1
                                            else:
                                                funcount['NOG'][item] = 1
                                        listfile.write('{}\t\t{}\n'.format(queryname, ortho_group))
                                        cog_list.append(ortho_group)
                            hsp = hit.find('Hit_hsps').find('Hsp')
                            Hsp_start_of_query = hsp.find('Hsp_query-from').text
                            Hsp_end_of_query = hsp.find('Hsp_query-to').text
                            calcstrand = int(Hsp_end_of_query) - int(Hsp_start_of_query)
                            if calcstrand >= 0:
                                Hsp_strand_of_query = '+'
                            else:
                                Hsp_strand_of_query = '-'
                            Hsp_start_of_hit = hsp.find('Hsp_hit-from').text
                            Hsp_end_of_hit = hsp.find('Hsp_hit-to').text
                            c2 = self.cog
                            DES = c2.find({'cog_id': ortho_group})
                            group_description = ''
                            for desresult in DES:
                                group_description = desresult['cog_description']
                            coverage = (abs(float(Hsp_end_of_hit) - float(Hsp_start_of_hit)) + 1) / (abs(float(region_end) - float(region_start)) + 1)
                            coverage = str("%.2f" % (coverage * 100)) + '%'
                            hsp_identity = hsp.find('Hsp_identity').text
                            calculen = abs(float(Hsp_end_of_hit) -
                                           float(Hsp_start_of_hit)) + 1
                            identity = str("%.2f" %
                                           ((float(hsp_identity) / calculen) * 100)) + '%'
                            hsp_positives = hsp.find('Hsp_positive').text
                            positives = str("%.2f" %
                                            ((float(hsp_positives) / calculen) * 100)) + '%'
                            tablefile.write(queryname + '\t' + query_length + '\t' + Hsp_start_of_query + '\t' + Hsp_end_of_query + '\t' + Hsp_strand_of_query + '\t' + hit_id + '\t' + hit_desc + '\t' + hit_length + '\t' + Hsp_start_of_hit + '\t' + Hsp_end_of_hit + '\t' + ortho_group + '\t' + group_description + '\t' + ','.join(func_cat) + '\t' + region_start + '\t' + region_end + '\t' + coverage + '\t' + identity + '\t' + positives + '\n')
            summaryfile.write('#Total seqs with COG/NOG:' + str(totalseq) + '\n')
            summaryfile.write('#Type\tfunctional_categories\tCOG\tNOG\n')
            for thekey in ['INFORMATION STORAGE AND PROCESSING', 'CELLULAR PROCESSES AND SIGNALING', 'METABOLISM', 'POORLY CHARACTERIZED']:
                for g in self.func_type[thekey]:
                    detail = self.func_decs[g]
                    category = '[' + g + ']' + ' ' + detail
                    try:
                        cogcount = funcount['COG'][g]
                    except KeyError:
                        cogcount = 0
                    try:
                        nogcount = funcount['NOG'][g]
                    except KeyError:
                        nogcount = 0
                    summaryfile.write(thekey + '\t' + category + '\t' +
                                      str(cogcount) + '\t' + str(nogcount) + '\n')

if __name__ == '__main__':
    string2cog = string2cog()
    string2cog.string2cog_by_mongo(string_xml=sys.argv[1], out_dir=sys.argv[2])


#test = string2cog()
#test.string2cog_by_mongo(string_xml='/mnt/ilustre/users/sanger-dev/sg-users/zengjing/denovo_rna/anno_rewrite/test.xml', out_dir='/mnt/ilustre/users/sanger-dev/sg-users/zengjing/denovo_rna/anno_rewrite/out_file/cog')
