# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'

import sys
from pymongo import MongoClient
from Bio.Blast import NCBIXML
import xml.etree.ElementTree as ET
import sqlite3

outpath=sys.argv[2]

dbfile = '/mnt/ilustre/users/sanger-dev/app/database/COG/cog_temporary_new_index_id.db'
#dbfile = '../cog_temporary_new_index_id.db'
conn = sqlite3.connect(dbfile)
c = conn.cursor()

# get cog_list.xls
'''
listfile = open('/cog_list.xls', 'w')
listfile.write('Query_name' + '\t' + 'COG' + '\t' + 'NOG' + '\n')

summaryfile = open('/cog_summary.xls', 'w')
tablefile = open('/cog_table.xls', 'w')
'''
listfile = open(outpath+'/cog_list.xls', 'w')
listfile.write('Query_name' + '\t' + 'COG' + '\t' + 'NOG' + '\n')

summaryfile = open(outpath+'/cog_summary.xls', 'w')
tablefile = open(outpath+'/cog_table.xls', 'w')

tablefile.write('#Query_name' + '\t' + 'Query_length' + '\t' + 'Hsp_start_of_query' + '\t' + 'Hsp_end_of_query' + '\t' + 'Hsp_strand_of_query' + '\t'
                'Hit_name' + '\t' + 'Hit_description' + '\t' + 'Hit_length' + '\t' + 'Hsp_start_of_hit' + '\t' +
                'Hsp_end_of_hit' + '\t' + 'COG/NOG_group' + '\t' + 'COG/NOG_group_description' + '\t' +
                'COG/NOG_group_categories' + '\t' + 'COG/NOG_region_start' + '\t' + 'COG/NOG_region_end' +
                '\t' + 'Coverage_of_COG/NOG_region' + '\t' + 'Identities_of_COG/NOG_region' + '\t' +
                'Positives_Identities_of_COG/NOG_region' + '\n')
'''
client = MongoClient('mongodb://192.168.10.189:27017')
db = client.sanger_biodb
coll = db.COG_mapping
coll2 = db.COG_description
'''

funccats = {'J': 'Translation, ribosomal structure and biogenesis',
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
            'R': 'General function prediction only', 'S': 'Function unknown'}

functype = {}
functype['INFORMATION STORAGE AND PROCESSING'] = sorted(
    ['J', 'A', 'K', 'L', 'B'])
functype['CELLULAR PROCESSES AND SIGNALING'] = sorted(
    ['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O'])
functype['METABOLISM'] = sorted(['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'])
functype['POORLY CHARACTERIZED'] = sorted(['R', 'S'])
# print functype
totalseq = 0
funcount = {}
#document = ET.parse('../transcript.fa_vs_string.xml')
document=ET.parse(sys.argv[1])
root = document.getroot()
identations = root.find('BlastOutput_iterations')
j = 0
for identation in identations.findall('Iteration'):
    j += 1
    querydef = identation.find('Iteration_query-def')
    queryname = querydef.text  # the query name
    queryname=queryname.split(' ')[0]
    querylen = identation.find('Iteration_query-len')
    query_length = querylen.text  # the query length
    iter_hits = identation.find('Iteration_hits')
    hits = iter_hits.findall('Hit')
    if len(hits) > 0:  # check if such name has hits
        l = []
        duplicate=False
        for hit in hits:
            hitdef = hit.find('Hit_id').text
            qid = hitdef.split(' ')[0]  # also the hit name
            hit_length = hit.find('Hit_len').text  # the hit length
            # find related COG
            #OG = coll.find({'_id': qid})
            OG=c.execute("SELECT * FROM COG_mapping WHERE _id = ?",(qid,)).fetchall()
            #if OG.count() == 1:
            if len(OG)>=1:
                if duplicate==False:
                    totalseq+=1
                    duplicate=True
                for og in OG:
                    ortho_group = og[3]  # COG/NOG group
                    if len(og[-1]) > 1:
                        func_cat = ','.split(og[-1])
                    else:
                        func_cat = [og[-1]]
                    try:
                        hit_desc = og[-2]
                    except KeyError:
                        hit_desc = ''
                    region_start = og[1]
                    region_end = og[2]
                    if ortho_group not in l:
                        if ortho_group.startswith('COG'):
                            ###deal with cog_summary.xls###
                            if funcount.has_key('COG'):
                                for item in func_cat:
                                    if funcount['COG'].has_key(item):
                                        funcount['COG'][item] += 1
                                    else:
                                        funcount['COG'][item] = 1
                            else:
                                newdic = {}
                                for item1 in func_cat:
                                    newdic[item1] = 1
                                funcount['COG'] = newdic
                            ###summary ends###
                            ###deal with cog_list.xls###
                            listfile.write(queryname + '\t' +
                                           ortho_group + '\t\t' + '\n')
                            l.append(ortho_group)
                        ###list ends###
                        else:
                            # if ortho_group not in l:
                            ###deal with cog_summary.xls###
                            if funcount.has_key('NOG'):
                                for oitem in func_cat:
                                    if funcount['NOG'].has_key(oitem):
                                        funcount['NOG'][oitem] += 1
                                    else:
                                        funcount['NOG'][oitem] = 1
                            else:
                                ndic = {}
                                for oitem1 in func_cat:
                                    ndic[oitem1] = 1
                                funcount['NOG'] = ndic
                            ###summary ends###
                            ###deal with cog_list.xls###
                            listfile.write(queryname + '\t\t' +
                                           ortho_group + '\n')
                            l.append(ortho_group)
                        ###list end###
                hsps = hit.find('Hit_hsps')
                hsp = hsps.find('Hsp')
                Hsp_start_of_query = hsp.find('Hsp_query-from').text
                Hsp_end_of_query = hsp.find('Hsp_query-to').text
                calcstrand = int(Hsp_end_of_query) - int(Hsp_start_of_query)
                if calcstrand >= 0:
                    Hsp_strand_of_query = '+'
                else:
                    Hsp_strand_of_query = '-'
                Hsp_start_of_hit = hsp.find('Hsp_hit-from').text
                Hsp_end_of_hit = hsp.find('Hsp_hit-to').text
                #DES = coll2.find({'_id': ortho_group})
                DES=c.execute("SELECT * FROM COG_description WHERE cog_id = ?",(ortho_group,)).fetchall()
                for desresult in DES:
                    group_description = desresult[1]
                # print Hsp_end_of_hit,Hsp_start_of_hit,region_end,region_start
                coverage = (abs(float(Hsp_end_of_hit) - float(Hsp_start_of_hit)) +
                            1) / (abs(float(region_end) - float(region_start)) + 1)
                #print coverage
                coverage = str("%.2f" % (coverage * 100)) + '%'
                #print coverage
                hsp_identity = hsp.find('Hsp_identity').text
                calculen = abs(float(Hsp_end_of_hit) -
                               float(Hsp_start_of_hit)) + 1
                identity = str("%.2f" %
                               ((float(hsp_identity) / calculen) * 100))+ '%'
                #print identity
                hsp_positives = hsp.find('Hsp_positive').text
                positives = str("%.2f" %
                                ((float(hsp_positives) / calculen) * 100)) + '%'
                #print positives
                tablefile.write(queryname + '\t' + query_length + '\t' + Hsp_start_of_query + '\t' +
                                Hsp_end_of_query + '\t' + Hsp_strand_of_query + '\t' + qid + '\t' +
                                hit_desc + '\t' + hit_length + '\t' + Hsp_start_of_hit + '\t' +
                                Hsp_end_of_hit + '\t' + ortho_group + '\t' + group_description + '\t' +
                                ','.join(func_cat) + '\t' + region_start + '\t' + region_end + '\t' +
                                coverage + '\t' + identity + '\t' + positives + '\n')


print funcount
###write whole cog_summary.xls###
summaryfile.write('#Total seqs with COG/NOG:' + str(totalseq) + '\n')
summaryfile.write('#Type' + '\t' + 'functional_categories' +
                  '\t' + 'COG' + '\t' + 'NOG' + '\n')
for thekey in ['INFORMATION STORAGE AND PROCESSING', 'CELLULAR PROCESSES AND SIGNALING', 'METABOLISM', 'POORLY CHARACTERIZED']:
    for g in functype[thekey]:
        detail = funccats[g]
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
###cog_summary.xls finished###
