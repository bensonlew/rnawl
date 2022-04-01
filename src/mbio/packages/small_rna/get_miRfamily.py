from collections import defaultdict

m2p = defaultdict(set)
with open('mirna_2_prefam.txt', 'r') as m2p_r:
    for line in m2p_r:
        line = line.strip().split('\t')
        m2p[line[1]].add(line[0])

id2mir = dict()
with open('mirna.txt', 'r') as mir:
    for line in mir:
        line = line.strip().split('\t')
        id2mir[line[0]] = line[2]

with open('mirna_prefam.txt', 'r') as mp_r, open('miRfamily.dat', 'w') as dat_w:
    for line in mp_r:
        line = line.strip().split('\t')
        if line[0] in m2p:
            for i in sorted(list(m2p[line[0]])):
                try:
                    dat_w.write(id2mir[i] + '\t' + '\t'.join(line[1:]) + '\n')
                except:
                    pass
