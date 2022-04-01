# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 2017.05.03

import sys
from pymongo import MongoClient


class CogAnnot(object):
    def __init__(self):
        self.client = MongoClient('mongodb://10.100.200.129:27017')
        self.cog_string = self.client.sanger_biodb.COG_String_V9
        self.cog = self.client.sanger_biodb.COG_V9
        self.func_type = {
            'INFORMATION STORAGE AND PROCESSING': ['A', 'B', 'J', 'K', 'L'],
            'CELLULAR PROCESSES AND SIGNALING': ['D', 'M', 'N', 'O', 'T', 'U', 'V', 'W', 'Y', 'Z'],
            'METABOLISM': ['C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q'],
            'POORLY CHARACTERIZED': ['R', 'S'],
        }
        self.func_decs = {
            'A': 'RNA processing and modification',
            'B': 'Chromatin structure and dynamics',
            'C': 'Energy production and conversion',
            'D': 'Cell cycle control, cell division, chromosome partitioning',
            'E': 'Amino acid transport and metabolism',
            'F': 'Nucleotide transport and metabolism',
            'G': 'Carbohydrate transport and metabolism',
            'H': 'Coenzyme transport and metabolism',
            'I': 'Lipid transport and metabolism',
            'J': 'Translation, ribosomal structure and biogenesis',
            'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility',
            'O': 'Posttranslational modification, protein turnover, chaperones',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
            'R': 'General function prediction only',
            'S': 'Function unknown',
            'T': 'Signal transduction mechanisms',
            'U': 'Intracellular trafficking, secretion, and vesicular transport',
            'V': 'Defense mechanisms',
            'W': 'Extracellular structures',
            'Y': 'Nuclear structure',
            'Z': 'Cytoskeleton'
        }

    def cog_annotation(self, string_table, out_dir):
        cog_list_path = out_dir + '/cog_list.xls'
        cog_table_path = out_dir + '/cog_table.xls'
        cog_summary_path = out_dir + '/cog_summary.xls'
        with open(string_table, "rb") as f, open(cog_list_path, 'wb') as listfile, open(cog_table_path, 'wb') as tablefile, open(cog_summary_path, 'wb') as summaryfile:
            listfile.write('Query_name\tCOG\tNOG\tKOG\n')
            tablefile.write('#Query_name\tQuery_length\tHsp_start_of_query\tHsp_end_of_query\tHsp_strand_of_query\tHit_name\tHit_description\tHit_length\tHsp_start_of_hit\tHsp_end_of_hit\tCOG/NOG/KOG_group\tCOG/NOG/KOG_group_description\tCOG/NOG/KOG_group_categories\tCOG/NOG/KOG_region_start\tCOG/NOG/KOG_region_end\tCoverage_of_COG/NOG/KOG_region\tIdentities_of_COG/NOG/KOG_region\tPositives_Identities_of_COG/NOG/KOG_region\n')
            c = self.cog_string
            totalseq = 0  # 注释COG/NOG/KOG的序列总数
            funcount = {'COG': {}, 'NOG': {}, 'KOG': {}}
            funlist = {'COG': {}, 'NOG': {}, 'KOG': {}}
            lines = f.readlines()
            query_list = []
            for line in lines[1:]:
                line = line.strip().split("\t")
                queryname = line[5]
                if queryname not in query_list:
                    query_list.append(queryname)
            for i in range(len(query_list)):
                queryname = query_list[i]
                cog, nog, kog, cog_list = [], [], [], []
                for line in lines[1:]:
                    line = line.strip().split("\t")
                    if queryname == line[5]:
                        hit_id = line[10]
                        query_length = line[6]
                        hit_length = line[11]
                        OG = c.find({'string_id': hit_id})
                        duplicate = False
                        if OG.count() >= 1:
                            if duplicate is False:
                                totalseq += 1
                                duplicate = True
                            for og in OG:
                                ortho_group = og['orthologous_group']  # COG/NOG/KOG group
                                if len(og['function_categories']) > 1:
                                    func_cat = ','.split(og['function_categories'])
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
                                        cog.append(ortho_group)
                                        for item in func_cat:
                                            if item in funcount['COG']:
                                                funcount['COG'][item] += 1
                                                funlist['COG'][item].append(queryname)
                                            else:
                                                funcount['COG'][item] = 1
                                                funlist['COG'][item] = []
                                                funlist['COG'][item].append(queryname)
                                        cog_list.append(ortho_group)
                                    elif ortho_group.startswith('NOG'):
                                        nog.append(ortho_group)
                                        for item in func_cat:
                                            if item in funcount['NOG']:
                                                funcount['NOG'][item] += 1
                                                funlist['NOG'][item].append(queryname)
                                            else:
                                                funcount['NOG'][item] = 1
                                                funlist['NOG'][item] = []
                                                funlist['NOG'][item].append(queryname)
                                        cog_list.append(ortho_group)
                                    elif ortho_group.startswith('KOG'):
                                        kog.append(ortho_group)
                                        for item in func_cat:
                                            if item in funcount['KOG']:
                                                funcount['KOG'][item] += 1
                                                funlist['KOG'][item].append(queryname)
                                            else:
                                                funcount['KOG'][item] = 1
                                                funlist['KOG'][item] = []
                                                funlist['KOG'][item].append(queryname)
                                        cog_list.append(ortho_group)
                                Hsp_start_of_query = line[7]
                                Hsp_end_of_query = line[8]
                                calcstrand = int(Hsp_end_of_query) - int(Hsp_start_of_query)
                                if calcstrand >= 0:
                                    Hsp_strand_of_query = '+'
                                else:
                                    Hsp_strand_of_query = '-'
                                Hsp_start_of_hit = line[12]
                                Hsp_end_of_hit = line[13]
                                c2 = self.cog
                                DES = c2.find({'cog_id': ortho_group})
                                group_description = ''
                                for desresult in DES:
                                    group_description = desresult['cog_description']
                                coverage = (abs(float(Hsp_end_of_hit) - float(Hsp_start_of_hit)) + 1) / (abs(float(region_end) - float(region_start)) + 1)
                                coverage = str("%.2f" % (coverage * 100)) + '%'
                                identity = line[3]
                                align_length = line[2]
                                similarity = line[4]
                                hsp_positives = float(similarity) / 100 * float(align_length)
                                hsp_identity = float(identity) / 100 * float(align_length)
                                calculen = abs(float(Hsp_end_of_hit) - float(Hsp_start_of_hit)) + 1
                                identity = str("%.2f" % ((float(hsp_identity) / calculen) * 100)) + '%'
                                positives = str("%.2f" % ((float(hsp_positives) / calculen) * 100)) + '%'
                                tablefile.write(queryname + '\t' + query_length + '\t' + Hsp_start_of_query + '\t' + Hsp_end_of_query + '\t' + Hsp_strand_of_query + '\t' + hit_id + '\t' + hit_desc + '\t' + hit_length + '\t' + Hsp_start_of_hit + '\t' + Hsp_end_of_hit + '\t' + ortho_group + '\t' + group_description + '\t' + ','.join(func_cat) + '\t' + region_start + '\t' + region_end + '\t' + coverage + '\t' + identity + '\t' + positives + '\n')
                if cog or nog or kog:
                    listfile.write('{}\t{}\t{}\t{}\n'.format(queryname, ';'.join(cog), ';'.join(nog), ';'.join(kog)))
            summaryfile.write('#Total seqs with COG/NOG/KOG:' + str(totalseq) + '\n')
            summaryfile.write('#Type\tfunctional_categories\tCOG\tNOG\tKOG\tcog_list\tnog_list\tkog_list\n')
            for thekey in ['INFORMATION STORAGE AND PROCESSING', 'CELLULAR PROCESSES AND SIGNALING', 'METABOLISM', 'POORLY CHARACTERIZED']:
                for g in self.func_type[thekey]:
                    detail = self.func_decs[g]
                    category = '[' + g + ']' + ' ' + detail  # A ""
                    try:
                        # cogcount = funcount['COG'][g]
                        coglist = list(set(funlist['COG'][g]))
                        cogcount = len(coglist)
                    except KeyError:
                        cogcount = 0
                        coglist = []
                    try:
                        # nogcount = funcount['NOG'][g]
                        noglist = list(set(funlist['NOG'][g]))
                        nogcount = len(noglist)
                    except KeyError:
                        nogcount = 0
                        noglist = []
                    try:
                        # kogcount = funcount['KOG'][g]
                        koglist = list(set(funlist['KOG'][g]))
                        kogcount = len(koglist)
                    except KeyError:
                        kogcount = 0
                        koglist = []
                    summaryfile.write(thekey + '\t' + category + '\t' +
                                      str(cogcount) + '\t' + str(nogcount) + '\t' + str(kogcount) + '\t' + ';'.join(coglist) + '\t' + ';'.join(noglist) + '\t' + ';'.join(koglist) + '\n')


if __name__ == '__main__':
    a = CogAnnot()
    a.cog_annotation(string_table=sys.argv[1], out_dir=sys.argv[2])
