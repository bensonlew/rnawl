# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

import xml.etree.ElementTree as ET
import sqlite3
from biocluster.config import Config
import re
import sys
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas


class KeggAnnotation(object):
    def __init__(self):

        self.mongodb = Config().get_mongo_client(mtype="ref_rna", ref=True)[Config().get_mongo_dbname("ref_rna", ref=True)]
        # self.mongodb = Config().mongo_client.sanger_biodb
        self.sqlitedb_path = Config().SOFTWARE_DIR + '/database/KEGG/ko/ko.db'
        self.sqlitedb = sqlite3.connect(self.sqlitedb_path)
        self.stat_info = dict()  # {'path_id': {'seqs': set([gene1,gene2,...]), 'koids': set(['K04345', 'K04432',...]), 'seqlist': set(['gene1(K02183)', 'gene2(K04345)',...])},...}
        self.cursor = self.sqlitedb.cursor()
        self.pathway_pic = Config().SOFTWARE_DIR + '/database/KEGG/pathway_map2/'

    def transversion(self, givenlist):
        return map(lambda x: x[0], givenlist)

    def kegg_by_sqlite(self, kegg_xml, kegg_table):
        """
        传入blast比对kegg的xml文件，获取kegg注释信息：kegg_table.xls，并获取统计pathway信息的self.stat_info,用于生成pathway_table.xls,pid.txt
        kegg_xml:blast比对kegg的xml文件
        kegg_table:kegg_table.xls结果文件
        """
        with open(kegg_table, 'wb') as tablefile:
            tablefile.write('#Query\tKO_ID (Gene id)\tKO_name (Gene name)\tHyperlink\tPaths\n')
            document = ET.parse(kegg_xml)
            root = document.getroot()
            identations = root.find('BlastOutput_iterations')
            for identation in identations.findall('Iteration'):
                query = identation.find('Iteration_query-def').text.split()[0]
                iter_hits = identation.find('Iteration_hits')
                hits = iter_hits.findall('Hit')
                duplicated = set()
                if len(hits) > 0:
                    for hit in hits:
                        gid = hit.find('Hit_id').text
                        koids = self.cursor.execute("SELECT koid FROM KoGenes WHERE gid = '%s'" % gid).fetchall()
                        koids = self.transversion(koids)
                        if len(koids) >= 1:
                            for theid in koids:
                                koresults = self.cursor.execute("SELECT * FROM Kos WHERE koid = '%s'" % theid).fetchall()
                                if len(koresults) >= 1:
                                    for ko in koresults:
                                        koid = ko[0]
                                        koname = ko[1]
                                        kohlink = 'http://www.genome.jp/dbget-bin/www_bget?ko:{}'.format(koid)
                                        pids = self.cursor.execute(
                                            "SELECT pid FROM KoPathways WHERE koid = '%s'" % koid)
                                        pids = self.transversion(pids)
                                        if len(pids) >= 1:
                                            newp = []
                                            for p in pids:
                                                newp.append('path:' + p)
                                                if p in self.stat_info:
                                                    self.stat_info[p]['seqs'].add(query)
                                                    self.stat_info[p]['seqlist'].add('{}({})'.format(query, koid))
                                                    self.stat_info[p]['koids'].add(koid)
                                                else:
                                                    s1 = set()
                                                    s1.add(query)
                                                    l1 = set()
                                                    l1.add('{}({})'.format(query, koid))
                                                    k1 = set()
                                                    k1.add(koid)
                                                    d1 = {}
                                                    d1['seqs'] = s1
                                                    d1['seqlist'] = l1
                                                    d1['koids'] = k1
                                                    self.stat_info[p] = d1
                                            pathinfo = ';'.join(newp)
                                        else:
                                            pathinfo = ''
                                        if koid not in duplicated:
                                            tablefile.write(
                                                query + '\t' + koid + '\t' + koname + '\t' + kohlink + '\t' + pathinfo + '\n')
                                            duplicated.add(koid)

    def get_pathway_result(self, stat_info, pathway, pidpath, db='sqlite'):
        """
        传入KEGG_pathway的相关统计信息的字典（即self.stat_info），生成pathway_table.xls,pid.txt
        stat_info：字典，存放着pathway的相关统计信息，见self.stat_info
        pathway:pathway_table.xls输出文件结果路径
        pidpath：pid.txt输出文件结果路径
        db:查询的数据库，可选sqlite或mongodb
        """
        # pathway = out_dir + '/pathway_table.xls'
        # pidpath = out_dir + '/pid.txt'
        with open(pathway, 'wb') as pathfile, open(pidpath, 'wb') as pidfile:
            pathfile.write('Pathway\tPathway_definition\tnumber_of_seqs\tseqs_kos/gene_list\tpathway_imagename\n')
            sdic = {}
            for thekey in stat_info:
                sdic[thekey] = len(stat_info[thekey]['seqs'])
            sortkey = list(reversed(sorted(sdic, key=lambda x: sdic[x])))
            for dickey in sortkey:
                pway = 'path:' + dickey
                if db == 'sqlite':
                    defseqrch = self.cursor.execute("SELECT name FROM Pathways WHERE pid = '%s'" % dickey).fetchone()
                if defseqrch:
                    pway_def = defseqrch[0]
                num_seq = len(stat_info[dickey]['seqs'])
                seq_list = ';'.join(list(stat_info[dickey]['seqlist']))
                pathfile.write(pway + '\t' + pway_def + '\t' + str(num_seq) + '\t' + seq_list + '\t' + dickey + '.png' + '\n')
                pidfile.write(dickey + '\t' + ';'.join(list(stat_info[dickey]['koids'])) + '\n')

    def get_kegg_layer(self, pathwayfile, layerfile, taxonomyfile, db='sqlite'):
        """
        传入pathway_table.xls统计信息，生成kegg_layer.xls,kegg_taxonomy.txt
        pathwayfile：pathway_table.xls统计信息
        layerfile:kegg_layer.xls结果文件
        taxonomyfile：kegg_taxonomy.txt结果文件
        db:查询的数据库，可选sqlite或mongodb
        """
        d = {}
        with open(pathwayfile, 'rb') as f:
            f.readline()
            for line in f:
                line = line.strip('\n').split('\t')
                seqnum = int(line[2])
                tlayer = line[1]
                if db == 'sqlite':
                    result = self.cursor.execute("SELECT * FROM kegg_layers WHERE third_layer = ?", (tlayer,))
                result = result.fetchall()
                if len(result) >= 1:
                    for document in result:
                        first_layer = document[0]
                        second_layer = document[1]
                        if first_layer in d:
                            if second_layer in d[first_layer]:
                                d[first_layer][second_layer] += seqnum
                            else:
                                d[first_layer][second_layer] = seqnum
                        else:
                            newdic = {}
                            newdic[second_layer] = seqnum
                            d[first_layer] = newdic
        with open(layerfile, 'wb') as l, open(taxonomyfile, 'wb') as t:
            for dickey in d:
                count = 0
                second_layer = []
                for skey in d[dickey]:
                    l.write('{}\t{}\t{}\n'.format(dickey, skey, d[dickey][skey]))
                    count += d[dickey][skey]
                    second_layer.append('--{}\t{}\n'.format(skey, d[dickey][skey]))
                t.write('{}\t{}\n'.format(dickey, count))
                for i in second_layer:
                    t.write(i)

    def get_pictrue(self, pidfile, out_dir):
        """
        传入pid.txt统计信息，生成pathway绘图文件夹
        pidfile：pid.txt统计信息
        out_dir:输出结果的目录
        """
        with open(pidfile, 'rb') as f:
            for line in f:
                line = line.strip('\n')
                if line:
                    line = line.split('\t')
                    pid = line[0]
                    koid = line[1].split(';')
                    l = []
                    if pid != 'ko00312' and pid != 'ko00351':
                        pathway = KGML_parser.read(open(self.pathway_pic + pid + '.kgml'))
                        pathway.image = self.pathway_pic + pid + '.png'
                        for ko in koid:
                            for degree in pathway.entries.values():
                                if re.search(ko, degree.name):
                                    l.append(degree.id)
                        for theid in l:
                            for graphic in pathway.entries[theid].graphics:
                                graphic.fgcolor = '#CC0000'
                        canvas = KGMLCanvas(pathway, import_imagemap=True)
                        canvas.draw(out_dir + '/' + pid + '.pdf')

if __name__ == '__main__':
    # python KeggAnnotation.py Trinity_vs_kegg.xml  kegg_table.xls pathway_table.xls pid.txt kegg_layer.xls kegg_taxonomy.xls ./kegg/pathway
    kegg = KeggAnnotation()
    kegg.kegg_by_sqlite(kegg_xml=sys.argv[1], kegg_table=sys.argv[2])
    kegg.get_pathway_result(stat_info=self.stat_info, pathway=sys.argv[3], pidpath=sys.argv[4], db='sqlite')
    kegg.get_kegg_layer(pathwayfile=sys.argv[3], layerfile=sys.argv[5], taxonomyfile=sys.argv[6])
    kegg.get_pictrue(pidfile=sys.argv[4], out_dir=sys.argv[7])
