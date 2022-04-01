# # -*- coding: utf-8 -*-
# # __author__ = 'qindanhua'

import subprocess
import os
import argparse

def primer(p, i, o):
    cmd = '%sprimer3_core < %s > %s' % (p, i, o)
    print cmd
    # os.system(cmd)
    try:
        subprocess.check_output(cmd, shell=True)
        return True
    except subprocess.CalledProcessError:
        print("运行primer3出错")
        return False

def main():
    _i = ''
    _o = ''
    parse = argparse.ArgumentParser()
    parse.add_argument('-p','--path', type=str, help='primer3_core path')
    parse.add_argument('-i','--infile',help='input file')
    parse.add_argument('-o','--outfile',help='outfile')
    args = parse.parse_args()
    _p = args.path
    _i = args.infile
    _o = args.outfile
    primer(_p, _i, _o) 

if __name__ == '__main__':
    main()