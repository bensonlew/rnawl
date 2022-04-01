# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# __version__ = 'v1.0'
# __last_modified__ = '20160112'
from __future__ import division
import argparse


parser = argparse.ArgumentParser(description='calculate q20 and q30 of a fastq')
parser.add_argument('-p', '--phred', help='phred of the input fastq, default:33')
parser.add_argument('-i', '--input', help='input fastq file', required=True)
parser.add_argument('-o', '--output', help='output stat file', required=True)
args = vars(parser.parse_args())

if args['phred']:
    phred = int(args['phred'])
else:
    phred = 33
infile = args['input']
outfile = args['output']

try:
    with open(infile, 'rb') as r:
        pass
except IOError:
    raise IOError('无法打开输入的fastq文件，请检查文件路径是否正确')

try:
    with open(outfile, 'wb') as w:
        pass
except IOError:
    raise IOError('无法生成输出文件，请检查是否有输出路径的写入权限')

with open(infile, 'rb') as r:
    line = r.next()
    if line[0] != "@":
        raise ValueError("fastq 文件格式不正确")
    line = r.next()
    line = r.next()
    if line[0] != '+':
        raise ValueError("fastq 文件格式不正确")

total_base = 0
q20_base = 0
q30_base = 0
count = 0
with open(infile, 'rb') as r:
    for line in r:
        line = r.next()
        line = r.next()
        line = r.next().rstrip('\n')
        count += 1
        if count % 10000 == 0:
            print "processing seq " + str(count)
        length = len(line)
        for i in xrange(length):
            value = ord(line[i])
            total_base += 1
            if value - phred >= 20:
                q20_base += 1
            if value - phred >= 30:
                q30_base += 1
q20_rate = (q20_base / total_base) * 100
q30_rate = (q30_base / total_base) * 100
q20_rate = "{:.2f}".format(q20_rate)
q30_rate = "{:.2f}".format(q30_rate)
with open(outfile, 'wb') as w:
    w.write(infile + "\t" + str(total_base) + '\t' + str(q20_base) + "\t" + str(q20_rate) + "\n")
    w.write(infile + "\t" + str(total_base) + '\t' + str(q30_base) + "\t" + str(q30_rate) + "\n")
