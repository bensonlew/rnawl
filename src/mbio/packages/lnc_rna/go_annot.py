# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan, qinjincheng'

import sys
import MySQLdb
import MySQLdb.cursors
import os

def main(query, host, user, passwd, db, dir_out):

    print 'INFO: start reading {}'.format(query)
    lines = open(query).readlines()
    d = dict()
    all_go = set()
    all_gene = list()
    for line in lines:
        if len(line) != 0:
            info = line.strip().split('\t')
            dickey = info[0]
            dicvalue = info[1].split(';')
            d[dickey] = dicvalue
            for item in dicvalue:
                all_go.add(item)
            if dickey not in all_gene:
                all_gene.append(dickey)

    print 'INFO: start connecting MySQLdb by class {}'.format(MySQLdb)
    sql = MySQLdb.connect(host=host, user=user, passwd=passwd, db=db)
    global cursor
    cursor = sql.cursor(MySQLdb.cursors.DictCursor)

    print 'INFO: start processing MySQLdb by operating {}'.format(cursor)
    go_seq = dict()
    go_des = dict()
    des_dict = dict()
    for idkey in d:
        molecular_function = str()
        biological_process = str()
        cellular_component = str()
        for goid in d[idkey]:
            cursor.execute("SELECT DISTINCT ancestor.*, graph_path.term1_id AS ancestor_id FROM term INNER JOIN graph_path ON (term.id=graph_path.term2_id) INNER JOIN term AS ancestor ON (ancestor.id=graph_path.term1_id)  WHERE term.acc='%s'" % goid)
            result = cursor.fetchall()
            id_list = list()
            des_list = list()
            m = 0  # molecular_function
            b = 0  # biological_process
            c = 0  # cellular_component
            for row in result:
                name = row['name']
                GO = row['acc']
                go_seq[(GO, idkey)] = 1
                go_des[GO] = name
                if GO == goid:
                    des = name
                if GO == 'GO:0003674':
                    m = 1
                if GO == 'GO:0008150':
                    b = 1
                if GO == 'GO:0005575':
                    c = 1
                if des_dict.has_key(name):
                    if des_dict[name].has_key(idkey):
                        if goid not in des_dict[name][idkey]:
                            des_dict[name][idkey].append(goid)
                    else:
                        nl = [goid]
                        des_dict[name][idkey] = nl
                else:
                    nd = {}
                    nd[idkey] = [goid]
                    des_dict[name] = nd
            if m == 1:
                molecular_function = goid + '|' + des + ';'
            if b == 1:
                biological_process = goid + '|' + des + ';'
            if c == 1:
                cellular_component = goid + '|' + des + ';'

    print 'INFO: start iterating and exporting files to {}'.format(dir_out)
    go1234file = open(os.path.join(dir_out, 'go_level1234_statistics.xls'), 'w')
    go123file = open(os.path.join(dir_out, 'go_level123_statistics.xls'), 'w')
    go12file = open(os.path.join(dir_out, 'go_level12_statistics.xls'), 'w')
    godetailfile = open(os.path.join(dir_out, 'go_detail.xls'), 'w')
    go1234file.write('GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\tGO Term (Lev4)\tGO ID (Lev4)\tSeq Number\tPercent\tSeq List\n')
    go123file.write('GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\tSeq Number\tPercent\tSeq List\n')
    go12file.write('GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tSeq Number\tPercent\tSeq List\n')
    godetailfile.write('GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\tGO Term (Lev4)\tGO ID (Lev4)\tSeq Number\tSeq List\n')
    majorterm = {
        'molecular_function': {1: {'molecular_function': 'GO:0003674'}},
        'biological_process': {1: {'biological_process': 'GO:0008150'}},
        'cellular_component': {1: {'cellular_component': 'GO:0005575'}},
    }
    for majorkey in majorterm:
        for i in range(1, 4):
            for term_name in majorterm[majorkey][i]:
                new_term = get_descendants(term_name)
                for n in new_term:
                    if majorterm[majorkey].has_key(i + 1):
                        majorterm[majorkey][i + 1][n] = new_term[n]
                    else:
                        majorterm[majorkey][i + 1] = {n: new_term[n]}
        for j in range(2, 5):
            zterm = dict()
            for lterm in majorterm[majorkey][j]:
                if des_dict.has_key(lterm):
                    zterm[lterm] = len(des_dict[lterm])
            genes = dict()
            genes_detail = dict()
            for nterm in zterm:
                seqs = list()
                seqs_detail = list()
                for l in des_dict[nterm]:
                    seqs.append(l)
                    new_info = l + '(' + ','.join(des_dict[nterm][l]) + ')'
                    if new_info not in seqs_detail:
                        seqs_detail.append(new_info)
                genes[nterm] = ';'.join(seqs)
                genes_detail[nterm] = ';'.join(seqs_detail)
            for tname in majorterm[majorkey][1]:
                yterm = get_descendants(tname)
                for y in yterm:
                    xterm = get_descendants(y)
                    if j == 2:
                        if zterm.has_key(y):
                            if zterm[y] != '':
                                pct = float(zterm[y]) / float(len(all_gene))
                                go12file.write(tname + '\t' + y + '\t' + majorterm[majorkey][2][y] + '\t' + str(zterm[y]) + '\t' + str('%.8f' % pct) + '\t' + genes[y] + '\n')
                    for x in xterm:
                        wterm = get_descendants(x)
                        if j == 3:
                            if zterm.has_key(x):
                                if zterm[x] != '':
                                    pct = float(zterm[x]) / float(len(all_gene))
                                    go123file.write(tname + '\t' + y + '\t' + majorterm[majorkey][2][y] + '\t' + x + '\t' + majorterm[majorkey][3][x] + '\t' + str(zterm[x]) + '\t' + str('%.8f' % pct) + '\t' + genes[x] + '\n')
                        for w in wterm:
                            if j == 4:
                                if zterm.has_key(w):
                                    if zterm[w] != '':
                                        pct = float(zterm[w]) / float(len(all_gene))
                                        go1234file.write(tname + '\t' + y + '\t' + majorterm[majorkey][2][y] + '\t' + x + '\t' + majorterm[majorkey][3][x] + '\t' + w + '\t' + majorterm[majorkey][4][w] + '\t' + str(zterm[w]) + '\t' + str('%.8f' % pct) + '\t' + genes[w] + '\n')
                                        godetailfile.write(tname + '\t' + y + '\t' + majorterm[majorkey][2][y] + '\t' + x + '\t' + majorterm[majorkey][3][x] + '\t' + w + '\t' + majorterm[majorkey][4][w] + '\t' + str(zterm[w]) + '\t' + genes_detail[w] + '\n')
    else:
        go1234file.close()
        go123file.close()
        go12file.close()
        godetailfile.close()
        print 'INFO: succeed in exporting 4 files to {}'.format(dir_out)

def get_descendants(term_name):
    go_term = dict()
    go_name = str()
    cursor.execute("SELECT DISTINCT descendant.acc, descendant.name, descendant.term_type FROM term INNER JOIN graph_path ON (term.id=graph_path.term1_id) INNER JOIN term AS descendant ON (descendant.id=graph_path.term2_id) WHERE distance = 1 and term.name='%s'" % term_name)
    results = cursor.fetchall()
    for row in results:
        go_name = row['name']
        go_acc = row['acc']
        go_term[go_name] = go_acc
    return go_term

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])