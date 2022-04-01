# -*- coding: utf-8 -*-
from collections import OrderedDict
from biocluster.iofile import File
from biocluster.config import Config
import os
from biocluster.core.exceptions import FileError
__author__ = "gdq"


class FastqListFile(File):
    """
    check fastq list file for quant toolbox
    """
    def __init__(self):
        super(FastqListFile, self).__init__()

    def to_dict(self):
        fastq = self.prop['path']
        fastq_info = OrderedDict()
        with open(fastq) as f:
            for line in f:
                if line.startswith('#') or (not line.strip()):
                    pass
                tmp_list = line.strip().split('\t')
                sample, fqs = tmp_list[0], tmp_list[1:]
                fastq_info.setdefault(sample, list())
                read1_list = [x.strip() for x in fqs[0].split(';')]
                fastq_info[sample].append(read1_list)
                if len(fqs) >= 2:
                    read2_list = [x.strip() for x in fqs[1].split(';')]
                    fastq_info[sample].append(read2_list)
        return fastq_info
