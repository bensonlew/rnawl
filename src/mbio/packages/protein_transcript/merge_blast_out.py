# -*- coding: utf-8 -*-
# author fengyitong 2019-01


import sys
import re


def merge_diomand_out(file_list, type_='xml', out=None):
    if type_ == 'xml':
        w = open(out, 'wb')
        head = False
        for f in file_list:
            with open(f, 'rb') as r:
                flag = False
                for line in r:
                    if not flag:
                        if not head:
                            w.write(line)
                        if re.match(r'^<Iteration>', line):
                            if head:
                                w.write(line)
                            else:
                                head = True
                            flag = True
                    elif not re.match(r'</BlastOutput', line) and line != '\n':
                        w.write(line)
                    else:
                        pass
        w.write(line_end)
        w.close()
    else:
        w = open(out, 'wb')
        head = True
        for f in file_list:
            with open(f, 'rb') as r:
                for line in r:
                    if not head:
                        w.write(line)
                        head = True
                    else:
                        w.write(line)
        w.close()

file_list = sys.argv[1]
type_ = sys.argv[2]
out = sys.argv[3]

line_end = '</BlastOutput_iterations>\n</BlastOutput>'

merge_diomand_out(file_list.split(';'), type_=type_, out=out)