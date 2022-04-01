# -*- coding: utf-8 -*-
# __author__ = 'gudeqing'

from biocluster.iofile import File
from collections import OrderedDict
from biocluster.core.exceptions import FileError

class FastqListFile(File):
    '''
    check fastq list file for quant toolbox
    '''
    def __init__(self):
        super(FastqListFile, self).__init__()

    def check(self):
        super(FastqListFile, self).check()

    def to_dict(self):
        fastq_info = OrderedDict()
        with open(self.prop['path']) as f:
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

    def prepare(self):
        self.pe1_fastq = dict()
        self.pe2_fastq = dict()
        self.se_fastq = dict()
        self.samples = set()
        for n, line in enumerate(open(self.path)):
            items = line.strip().split('\t')
            if len(items) == 3:
                self.pe1_fastq[items[0]] = items[1]
                self.pe2_fastq[items[0]] = items[2]
            elif len(items) == 2:
                self.se_fastq[items[0]] = items[1]
            else:
                continue
            self.samples.add(items[0])
        if len(self.pe1_fastq) == n + 1 and len(self.pe2_fastq) == n + 1 and not self.se_fastq:
            return 'PE'
        elif not self.pe1_fastq and not self.pe2_fastq and len(self.se_fastq) == n + 1:
            return 'SE'
        else:
            raise FileError('can not determine sequencing type, check {}'.format(self.path))
