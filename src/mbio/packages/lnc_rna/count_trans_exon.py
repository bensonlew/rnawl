# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, qinjincheng'

from optparse import OptionParser
import os
import fileinput

parser = OptionParser(description='Export files containing number of downstream products')
parser.add_option('-i', '--dir_in', dest='dir_in', help='Input directory containing files about counts relationship')
parser.add_option('-s', '--steps', dest='steps', help='Step size of the counts')
parser.add_option('-o', '--dir_out', dest='dir_out', help='Output directory for results')
(opts, args) = parser.parse_args()

def main(dir_in, steps, dir_out):
    for file_in, step, file_tmp, file_out in transfer(dir_in, steps, dir_out):
        print 'INFO: start computing distribution of transcripts and exons number by {} step'.format(step)
        producer(file_in, step, file_tmp, file_out)
        if os.path.isfile(file_out):
            print 'INFO: succeed in exporting number of downstream products in {}'.format(file_out)
        else:
            raise Exception('ERROR: fail to count export number of downstream products in {}'.format(file_out))
    else:
        print 'INFO: succeed in exporting files containing number of downstream products'

def transfer(dir_in, steps, dir_out):
    files = ['old_genes.gtf.trans', 'old_transcripts.gtf.exon', 'new_genes.gtf.trans', 'new_transcripts.gtf.exon']
    for file_in in os.listdir(dir_in):
        if file_in in files:
            file_path = os.path.join(dir_out, file_in)
            for step in steps.split(','):
                middle_txt = os.path.join(dir_out, '{}_{}.test.txt'.format(file_in, step))
                final_txt = os.path.join(dir_out, '{}_{}.txt'.format(file_in, step))
                yield file_path, int(step), middle_txt, final_txt

def producer(file_in, step, file_tmp, file_out):
    if os.path.getsize(file_in) == 0:
        open(file_out, 'w').close()
    with open(file_in) as r:
        lines = r.readlines()
        if len(lines) >= 3:
            f1 = open(file_tmp, 'w')
            f2 = open(file_out, 'w')
            dic = dict()
            group = 0
            list_set = set()
            for line in fileinput.input(file_in):
                lines = line.strip().split('\t')
                dic[lines[0]] = lines[1]
                list_set.add(int(lines[1]))
            sort_set = list(list_set)
            for i in sorted(sort_set):
                value = 0
                ids_list = list()
                for key in dic.keys():
                    if int(dic[key]) == i:
                        value += 1
                        ids_list.append(key)
                if int(step) == 1:
                    f1.write('{}\t{}\t{}\n'.format(i, value, ids_list))
                    f2.write('{}\t{}\t{}\n'.format(i, value, ids_list))
                else:
                    f1.write('{}\t{}\t{}\n'.format(i, value, ids_list))
            f1.close()
            if int(step) != 1:
                max_num = sorted(sort_set)[-1]
                group = max_num / step
                mod = max_num % step
                if mod != 0:
                    group += 1
            for n in range(group):
                f3 = open(file_tmp)
                final_value = 0
                final_ids = list()
                for line in f3:
                    num = line.strip().split('\t')[0]
                    value = line.strip().split('\t')[1]
                    id_s = line.strip().split('\t')[2]
                    ids_list = id_s.strip('[').strip(']').strip("'").split(',')
                    if (int(num) > (n * step)) and (int(num) <= ((n + 1) * step)):
                        final_value += int(value)
                        for id_num in ids_list:
                            final_ids.append(id_num)
                if n == 0:
                    area_line = '{}~{}\t{}\t{}\n'.format(n * step, (n + 1) * step, final_value, final_ids)
                else:
                    area_line = '{}~{}\t{}\t{}\n'.format(n * step + 1, (n + 1) * step, final_value, final_ids)
                f2.write(area_line)
                f3.close()
            f2.close()

if __name__ == '__main__':
    if opts.dir_in and opts.steps and opts.dir_out:
        main(opts.dir_in, opts.steps, opts.dir_out)
    else:
        parser.print_help()