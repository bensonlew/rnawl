# -*- coding: utf-8 -*-
# __author__ = 'gaohao'


import pandas as pd
from optparse import OptionParser

def file_config(input,res,output):
    file = open(output,'w')
    f= pd.read_csv(input, sep='\t', header=None,lineterminator='\n')
    num = f.iloc[:, 5].max()
    with open(input,'r') as f:
         lines = f.readlines()
         for i in range(0, len(lines)):
             if i ==0:
                 line = lines[i].strip().split('\t')
                 file.write('max_rd_len=%s \n'%num)
                 file.write('[LIB]\n')
                 file.write('avg_ins=%s \n' % line[4])
                 file.write('reverse_seq=0\n')
                 file.write('asm_flags=3\n')
                 file.write('rank=1\n')
                 file.write('q1=%s\n' %line[1])
                 file.write('q2=%s\n' %line[2])
                 #file.write('q=%s\n' % line[3])
             else:
                 file.write('[LIB]\n')
                 file.write('avg_ins=%s \n' % line[4])
                 file.write('reverse_seq=0\n')
                 file.write('asm_flags=3\n')
                 file.write('rank=1\n')
                 file.write('q1=%s\n' % line[1])
                 file.write('q2=%s\n' % line[2])
                 #file.write('q=%s\n' % line[3])

    with open(res,'r') as g:
         lines = g.readlines()
         for i in range(0, len(lines)):
             line = lines[i].strip().split('\t')
             file.write('[LIB]\n')
             file.write('avg_ins=%s \n' % line[3])
             file.write('reverse_seq=1\n')
             file.write('asm_flags=2\n')
             file.write('rank=2\n')
             file.write('q1=%s\n' %line[1])
             file.write('q2=%s\n' %line[2])
    file.close()

def file_config2(input,output):
    file = open(output,'w')
    f = pd.read_csv(input, sep='\t',header=None,lineterminator='\n')
    num =f.iloc[:,5].max()
    with open(input,'r') as f:
         lines = f.readlines()
         for i in range(0, len(lines)):
             if i ==0:
                 line = lines[i].strip().split('\t')
                 file.write('max_rd_len=%s \n'%num)
                 file.write('[LIB]\n')
                 file.write('avg_ins=%s \n' % line[4])
                 file.write('reverse_seq=0\n')
                 file.write('asm_flags=3\n')
                 file.write('rank=1\n')
                 file.write('q1=%s\n' %line[1])
                 file.write('q2=%s\n' %line[2])
                 #file.write('q=%s\n' % line[3])
             else:
                 file.write('[LIB]\n')
                 file.write('avg_ins=%s \n' % line[4])
                 file.write('reverse_seq=0\n')
                 file.write('asm_flags=3\n')
                 file.write('rank=1\n')
                 file.write('q1=%s\n' % line[1])
                 file.write('q2=%s\n' % line[2])
                 #file.write('q=%s\n' % line[3])
    file.close()

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser()
    parser.add_option('-p', dest ='pe',metavar='[PE_list.xls]')
    parser.add_option('-m', dest ='mp',metavar='[MP_list.xls]')
    parser.add_option('-o', dest ='output',metavar='[OUT FILE]')
    (options,args) = parser.parse_args()
    if options.mp != None:
        file_config(options.pe, options.mp, options.output)
    elif options.mp == None:
        file_config2(options.pe, options.output)
    else:
        print "python config_file.py -p file -m file -o outfile \npython config_file.py -p file -o outfile \n"


if __name__=='__main__':
    main()
