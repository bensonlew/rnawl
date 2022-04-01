# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os
import re

parser = OptionParser(description='Process fasta file of referential miRNA')
parser.add_option('-i', '--input', dest='input' ,help='Input fasta file')
parser.add_option('-o', '--output', dest='output', help='Output fasta file')
parser.add_option('-s', '--species', dest='species', help='Abbreviation of target species')
parser.add_option('-t', '--transform', action='store_true', default=False, help='Make base substitution from U to T')
(opts, args) = parser.parse_args()

class MirRefProcess(object):
    '''
    Process fasta file of referential miRNA
    '''
    def __init__(self, fasta, target='', transform=False, drop_info=True, lowercase=True):
        self.check(fasta=fasta)
        self.fasta = fasta
        self.target = target
        self.transform = transform
        self.drop_info = drop_info
        self.lowercase = lowercase

    def check(self, fasta):
        if not os.path.isfile(fasta):
            raise Exception('Error: {} is not a file'.format(fasta))
        try:
            first_line = open(fasta).readline()
            if first_line[0] != '>':
                print 'Warning: find format error in {}'.format(fasta)
        except:
            raise Exception('Error: {} is not readable'.format(fasta))

    def fa_to_dict(self, fa):
        genome = dict()
        seq_name_tmp = str()
        with open(fa) as f:
            for line in f:
                if line[0] == '>':
                    seq_name_now = line.strip()
                    if self.drop_info:
                        m = re.match(r'(>\S+)', seq_name_now)
                        if m:
                            seq_name_now = m.group(1)
                            if self.lowercase:
                                seq_name_now = seq_name_now.lower()
                            genome[seq_name_now] = str()
                        else:
                            print 'Warning: find illegal seq name in {}'.format(line)
                            genome[seq_name_now] = str()
                    else:
                        genome[seq_name_now] = str()
                else:
                    read = line.strip()
                    if seq_name_now:
                        genome[seq_name_now] += read
                    else:
                        genome[seq_name_tmp] += read
        if genome.has_key(''):
            print 'Warning: find not included reads with {} length'.format(genome[''])
        return genome

    def pick_target(self, genome, abbr):
        ret = dict()
        for k, v in genome.iteritems():
            m = re.match(r'(>\S+)', k)
            if m and (abbr in m.group(1)):
                ret[k] = v
        return ret

    def base_substitution(self, genome, old, new):
        ret = dict()
        for k, v in genome.iteritems():
            ret[k] = v.replace(old, new)
        return ret

    def dict_to_fa(self, genome, fa):
        with open(fa, 'w') as f:
            for k, v in genome.iteritems():
                f.write('{}\n{}\n'.format(k, v))

    def main(self, output):
        self.genome = self.fa_to_dict(fa=self.fasta)
        if self.target:
            self.genome = self.pick_target(genome=self.genome, abbr=self.target)
        else:
            print 'Warning: select all seq while no target'
        if self.transform:
            self.genome = self.base_substitution(genome=self.genome, old='U', new='T')
        self.dict_to_fa(genome=self.genome, fa=output)

if __name__ == '__main__':
    if opts.input and opts.output:
        if opts.transform:
            inst = MirRefProcess(fasta=opts.input, target=opts.species, transform=opts.transform)
            inst.main(output=opts.output)
        else:
            inst = MirRefProcess(fasta=opts.input, target=opts.species)
            inst.main(output=opts.output)
    else:
        parser.print_help()