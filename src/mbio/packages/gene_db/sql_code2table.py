# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import xml.etree.ElementTree as ET
import gzip


sql_file = sys.argv[1]
with open(sql_file, 'r') as mart_sql_f:
    schema = False
    table_name=None
    header = list()
    for line in mart_sql_f:
        # print "one line"
        # print line
        if line.startswith("CREATE TABLE"):
            schema = True
            header = list()
            table_name = line.split('`')[1]

        if line.startswith(")"):
            schema = False
            with open(table_name + ".tsv", 'w') as f:
                f.write("\t".join(header) + "\n")

        if schema:
            field = line.split('`')
            # print field
            if len(field) > 2 and field[0] == '  ':
                header.append(field[1])

        if line.startswith("INSERT INTO"):
            cols = line.split(' ')
            table_name = cols[2].strip("`")
            with open(table_name + ".tsv", 'aw') as f:
                data_list = cols[4][1:].split("),(")
                for data in data_list:
                    # data_set = eval(data)
                    # data_str = map(str, data_set)
                    f.write(data.replace(",", "\t") + "\n")
