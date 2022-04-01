import json
import logging
import re
import sys

import regex
from biocluster.config import Config


if __name__ == '__main__':
    hmm_file = sys.argv[1]
    with open(hmm_file, 'r') as f, open(hmm_file + "name_acc_des.txt", 'w') as fo:
        for line in f:
            if line.startswith("HMMER3/f"):
                name = ""
                acc = ""
                desc = ""
            if line.startswith("NAME"):
                name = line.rstrip("\n")[6:]

            if line.startswith("ACC"):
                acc = line.rstrip("\n")[6:]
            if line.startswith("DESC"):
                desc = line.rstrip("\n")[6:]

                fo.write("\t".join([name, acc, desc]) + "\n")
