# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'
# last_modified: 20160720

import sys
import MySQLdb
import MySQLdb.cursors

db = MySQLdb.connect(host=sys.argv[2], user=sys.argv[
                     3], passwd=sys.argv[4], db="b2gdb")
cur = db.cursor(MySQLdb.cursors.DictCursor)

f = open(sys.argv[1]).read().split('\n')
d = {}
allgo = set()
allgene = []
for linerecord in f:
    if len(linerecord) != 0:
        info = linerecord.split('\t')
        dickey = info[0]
        dicvalue = info[1].split(';')
        d[dickey] = dicvalue
        for item in dicvalue:
            allgo.add(item)
        if dickey not in allgene:
            allgene.append(dickey)

go_seq = {}
go_des = {}
seqlist = {}
for idkey in d:
    molecular_function = ""
    biological_process = ""
    cellular_component = ""
    for goid in d[idkey]:
        cur.execute("SELECT DISTINCT ancestor.*, graph_path.term1_id AS ancestor_id FROM term INNER JOIN graph_path ON (term.id=graph_path.term2_id) INNER JOIN term AS ancestor ON (ancestor.id=graph_path.term1_id)  WHERE term.acc='%s'" % goid)
        result = cur.fetchall()
        idlist = []
        deslist = []
        m = 0  # molecular_function
        b = 0  # biological_process
        c = 0  # cellular_component
        for row in result:
            reid = row['id']
            name = row['name']
            GO = row['acc']
            go_seq[(GO, idkey)] = 1
            go_des[GO] = name
            if GO == goid:
                des = name
            if GO == "GO:0003674":
                m = 1
            if GO == "GO:0008150":
                b = 1
            if GO == "GO:0005575":
                c = 1
            if seqlist.has_key(name):
                if seqlist[name].has_key(idkey):
                    if goid not in seqlist[name][idkey]:
                        seqlist[name][idkey].append(goid)
                else:
                    nl = [goid]
                    seqlist[name][idkey] = nl
            else:
                nd = {}
                nd[idkey] = [goid]
                seqlist[name] = nd
        if m == 1:
            molecular_function = goid + '|' + des + ';'
        if b == 1:
            biological_process = goid + '|' + des + ';'
        if c == 1:
            cellular_component = goid + '|' + des + ';'


def getDescendants(factor1):
    goterm = {}
    goname = ''
    cur.execute("SELECT DISTINCT descendant.acc, descendant.name, descendant.term_type FROM  term  INNER JOIN graph_path ON (term.id=graph_path.term1_id)  INNER JOIN term AS descendant ON (descendant.id=graph_path.term2_id) WHERE distance = 1 and term.name='%s'" % factor1)
    results = cur.fetchall()
    for row in results:
        goname = row['name']
        GOnum = row['acc']
        goterm[goname] = GOnum
    return goterm

go1234file = open('go1234level_statistics.xls', 'w')
go123file = open('go123level_statistics.xls', 'w')
go12file = open('go12level_statistics.xls', 'w')
go1234file.write(
    'GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\tGO Term (Lev4)\tGO ID (Lev4)\tSeq Number\tPercent\tSeq List\n')
go123file.write('GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\tSeq Number\tPercent\tSeq List\n')
go12file.write('GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tSeq Number\tPercent\tSeq List\n')
# x=4
godetail = open('go_detail.xls', 'w')
godetail.write(
    'GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\tGO Term (Lev4)\tGO ID (Lev4)\tSeq Number\tSeq List\n')
# x=4
tree = ''
majorterm = {'molecular_function': {1: {'molecular_function': "GO:0003674"}}, 'biological_process': {
    1: {'biological_process': "GO:0008150"}}, 'cellular_component': {1: {'cellular_component': "GO:0005575"}}}
for majorkey in majorterm:
    for i in range(1, 4):
        for tname in majorterm[majorkey][i]:
            newterm = getDescendants(tname)
            for m in newterm:
                if majorterm[majorkey].has_key(i + 1):
                    majorterm[majorkey][i + 1][m] = newterm[m]
                else:
                    newdic = {m: newterm[m]}
                    majorterm[majorkey][i + 1] = newdic
    for j in range(2, 5):
        zterm = {}
        for lterm in majorterm[majorkey][j]:
            if seqlist.has_key(lterm):
                zterm[lterm] = len(seqlist[lterm])
        genes = {}
        geneswithdetail = {}
        for nterm in zterm:
            seqs = []
            seqsdetail = []
            for x in seqlist[nterm]:
                seqs.append(x)
                newinfo = x + '(' + ','.join(seqlist[nterm][x]) + ')'
                if newinfo not in seqsdetail:
                    seqsdetail.append(newinfo)
            genes[nterm] = ';'.join(seqs)
            geneswithdetail[nterm] = ';'.join(seqsdetail)
        for ttname in majorterm[majorkey][1]:
            newterm2 = getDescendants(ttname)
            for m2 in newterm2:
                newterm3 = getDescendants(m2)
                if j == 2:
                    if zterm.has_key(m2):
                        if zterm[m2] != '':
                            pcent = float(zterm[m2]) / float(len(allgene))
                            go12file.write(ttname + '\t' + m2 + '\t' + majorterm[majorkey][2][m2] + '\t' + str(zterm[m2]) + '\t' + str("%.8f" % pcent) + '\t' + genes[m2] + '\n')
                for m3 in newterm3:
                    newterm4 = getDescendants(m3)
                    if j == 3:
                        if zterm.has_key(m3):
                            if zterm[m3] != '':
                                pcent = float(zterm[m3]) / float(len(allgene))
                                go123file.write(ttname + '\t' + m2 + '\t' + majorterm[majorkey][2][m2] + '\t' + m3 + '\t' + majorterm[majorkey][
                                                 3][m3] + '\t' + str(zterm[m3]) + '\t' + str("%.8f" % pcent) + '\t' + genes[m3] + '\n')
                    for m4 in newterm4:
                        if j == 4:
                            if zterm.has_key(m4):
                                if zterm[m4] != '':
                                    pcent = float(zterm[m4]) / float(len(allgene))
                                    go1234file.write(ttname + '\t' + m2 + '\t' + majorterm[majorkey][2][m2] + '\t' + m3 + '\t' + majorterm[majorkey][
                                                     3][m3] + '\t' + m4 + '\t' + majorterm[majorkey][4][m4] + '\t' + str(zterm[m4]) + '\t' + str("%.8f" % pcent) + '\t' +genes[m4] + '\n')
                                    godetail.write(ttname + '\t' + m2 + '\t' + majorterm[majorkey][2][m2] + '\t' + m3 + '\t' + majorterm[majorkey][3][
                                                   m3] + '\t' + m4 + '\t' + majorterm[majorkey][4][m4] + '\t' + str(zterm[m4]) + '\t' + geneswithdetail[m4] + '\n')
