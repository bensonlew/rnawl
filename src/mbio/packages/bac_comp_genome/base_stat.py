# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
#20190902

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from Bio import SeqIO

def break_seq(infile, outfile):
    """
    根据输入序列，将所有的N打断
    :param infile:
    :param outfile:
    :return:
    """





if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', metavar='[cluster_file]',required=True, help='input cluster file')
    parse.add_argument('-o', metavar='[output_file]', required=True, help='input output_file')
    args = parse.parse_args()
    infile = args.i
    outfile = args.o
    break_seq(infile,outfile)