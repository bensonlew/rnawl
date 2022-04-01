import sys

if len(sys.argv) != 3:
    exit('USAGE: %s infile outfile' %sys.argv[0])

protein_list = list()
species_dict = dict()

with open(sys.argv[1], 'r') as in_r:
    in_r.readline()
    for line in in_r:
        line = line.strip('\n').split('\t')
        if line[5] in protein_list:
            continue
        else:
            protein_list.append(line[5])
            species = '___'
            if u'[' in line[-1]:
                species = line[-1].split('[')[1].strip(']')
            if u'var.' in species:
                species = species.split('var.')[0].strip()
            if not species in species_dict:
                species_dict[species] = 0
            species_dict[species] += 1

with open(sys.argv[2], 'w') as o_w:
    o_w.write('species\thit_num\n')
    for spe in species_dict:
        o_w.write(spe + '\t' + str(species_dict[spe]) + '\n')