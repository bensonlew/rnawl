import sys


def goSplit(infile):
    f = open(infile).read().split('\n')
    body = f[1:]
    lev2 = open('go2level.xls', 'w')
    lev3 = open('go3level.xls', 'w')
    lev4 = open('go4level.xls', 'w')
    allnum = set()
    for i in range(1, 4):
        d = {}
        for record in body:
            item = record.split('\t')
            if item != ['']:
                Type = item[0]
                GOid = item[i * 2 - 1]
                GOterm = item[i * 2]
                seqnum = int(item[7])
                seqlist = item[8].split(';')
                #seqlist += '(' + GOid + ')'
                # if i == 1:
                #allnum += seqnum
                if d.has_key(GOid):
                    #d[GOid]['number'] += seqnum
                    for seq in seqlist:
                        posi1=seq.find('(')
                        query=seq[:posi1]
                        allnum.add(query)
                        gos=seq[posi1+1:-1].split(',')
                        if d[GOid]['SeqList'].has_key(query):
                            for goitem in gos:
                                d[GOid]['SeqList'][query].add(goitem)
                        else:
                            goset=set()
                            for gitem in gos:
                                goset.add(gitem)
                            d[GOid]['SeqList'][query]=goset
                        #d[GOid]['SeqList'].add(seq)
                else:
                    newquerydic={}
                    for qseq in seqlist:
                        posi2=qseq.find('(')
                        query2=qseq[:posi2]
                        qgos=qseq[posi2+1:-1].split(',')
                        newgoset=set()
                        allnum.add(query2)
                        for gid in qgos:
                            newgoset.add(gid)
                        newquerydic[query2]=newgoset
                    d[GOid] = {'term_type': Type,
                               'term': GOterm, 'SeqList': newquerydic}
        if i == 1:
            lev2.write('term_type' + '\t' + 'term' + '\t' + 'GO' + '\t' +
                       'number' + '\t' + 'percent' + '\t' + 'sequence' + '\n')
        if i == 2:
            lev3.write('term_type' + '\t' + 'term' + '\t' + 'GO' + '\t' +
                       'number' + '\t' + 'percent' + '\t' + 'sequence' + '\n')
        if i == 3:
            lev4.write('term_type' + '\t' + 'term' + '\t' + 'GO' + '\t' +
                       'number' + '\t' + 'percent' + '\t' + 'sequence' + '\n')
        for theid in d:
            s = d[theid]['SeqList']
            '''
            newdic = {}
            for t in s:
                posi1 = t.find('(')
                sid = t[:posi1]
                gofs = t[posi1 + 1:-1].split(',')
                if newdic.has_key(sid):
                    for g in gofs:
                        if g not in newdic[sid]:
                            newdic[sid].append(g)
                else:
                    newlist = []
                    for k in gofs:
                        if k not in newlist:
                            newlist.append(k)
                    newdic[sid] = newlist
            seqinfo = []
            for transkey in newdic:
                if transkey != '':
                    seqinfo.append(
                        transkey + '(' + ','.join(newdic[transkey]) + ')')
            '''
            pcent = float(len(s)) / len(allnum)
            qlist=[]
            for t in s:
                allgo=','.join(list(s[t]))
                qlist.append(t+'('+allgo+')')
            if i == 1:
                lev2.write(d[theid]['term_type'] + '\t' + theid + '\t' + d[theid]['term'] + '\t' + str(
                    len(s)) + '\t' + str("%.8f" % pcent) + '\t' + ';'.join(qlist) + '\n')
            if i == 2:
                lev3.write(d[theid]['term_type'] + '\t' + theid + '\t' + d[theid]['term'] + '\t' + str(
                    len(s)) + '\t' + str("%.8f" % pcent) + '\t' + ';'.join(qlist) + '\n')
            if i == 3:
                lev4.write(d[theid]['term_type'] + '\t' + theid + '\t' + d[theid]['term'] + '\t' + str(
                    len(s)) + '\t' + str("%.8f" % pcent) + '\t' + ';'.join(qlist) + '\n')
            '''
            if i == 1:
                lev2.write(d[theid]['term_type'] + '\t' + d[theid]['term'] + '\t' + ';'.join(seqinfo) +
                           '\t' + theid + '\t' + str(d[theid]['number']) + '\t' + str("%.3f" % pcent) + '\n')
            if i == 2:
                lev3.write(d[theid]['term_type'] + '\t' + d[theid]['term'] + '\t' + ';'.join(seqinfo) +
                           '\t' + theid + '\t' + str(d[theid]['number']) + '\t' + str("%.3f" % pcent) + '\n')
            if i == 3:
                lev4.write(d[theid]['term_type'] + '\t' + d[theid]['term'] + '\t' + ';'.join(seqinfo) +
                           '\t' + theid + '\t' + str(d[theid]['number']) + '\t' + str("%.3f" % pcent) + '\n')
            '''

goSplit(sys.argv[1])
