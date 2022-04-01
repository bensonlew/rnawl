# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modified: liubinxu 20180901

import sys
import os
import re
from biocluster.config import Config


if __name__ == '__main__':
    f = open(sys.argv[1]).read().split('\n')
    out_file = sys.argv[2]

    reactome_path_file = open(out_file, 'w')

    for linerecord in f:
        if len(linerecord) != 0:
            info = linerecord.split('\t')
            gene_id = info[0]
            path_link = info[1]
            paths = [l.split("&")[0] for l in path_link.split(";")]
            reactome_path_file.write("{}\t{}\n".format(gene_id, ";".join(list(set(paths)))))
    reactome_path_file.close
