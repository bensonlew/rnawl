#!/usr/bin/env python

# Copyright 2013-2014 Mitchell Stanton-Cook Licensed under the
# Educational Community License, Version 2.0 (the "License"); you may
# not use this file except in compliance with the License. You may
# obtain a copy of the License at
#
#     http://www.osedu.org/licenses/ECL-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an "AS IS"
# BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
# or implied. See the License for the specific language governing
# permissions and limitations under the License.

"""
Generate gff file from EMBL/Genbank for QUAST
"""

__title__        = 'to_gff'
__version__      = '0.1.1'
__description__  = "Generate gff file from EMBL/Genbank for QUAST"
__author__       = 'Mitchell Stanton-Cook'
__license__      = 'ECL 2.0'
__author_email__ = "m.stantoncook@gmail.com"
__url__         = 'http://github.com/mscook/to_gff'


import sys, os, traceback, argparse
import time

from BCBio import GFF
from Bio import SeqIO

epi = "Licence: %s by %s <%s>" % (__license__,
                                  __author__,
                                  __author_email__)
__doc__ = " %s v%s - %s (%s)" % (__title__,
                                 __version__,
                                 __description__,
                                 __url__)

def rev(seq):
    base_trans = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N', 'a':'t', 'c':'g', 't':'a', 'g':'c', 'n':'n'}
    rev_seq = list(reversed(seq))
    rev_seq_list = [base_trans[k] for k in rev_seq]
    rev_seq = ''.join(rev_seq_list)
    return(rev_seq)

def to_GFF(args):
    """
    Convert a GenBank or EMBL file to GFF

    Mainly useful for QUAST (Quality Assessment Tool for Genome Assemblies)

    :param args: an argparse args list
    """
    args.in_file  = os.path.expanduser(args.in_file)
    args.out_file = args.in_file

    in_type = "genbank"
    base =  './'
    gff = os.path.splitext(os.path.basename(args.out_file))[0]+'.gff'
    gff_out =  os.path.join(base, gff)
    Pfasta = os.path.splitext(os.path.basename(args.out_file))[0]+'.faa'
    Pfasta_out = os.path.join(base, Pfasta)
    Nfasta = os.path.splitext(os.path.basename(args.out_file))[0]+'.cds'
    Nfasta_out = os.path.join(base, Nfasta)
    fasta = os.path.splitext(os.path.basename(args.out_file))[0]+'.fasta'
    fasta_out =  os.path.join(base, fasta)

    if args.gff:
        with open(args.in_file) as fin:
            with open(gff_out, 'w') as gff_out:
                GFF.write(SeqIO.parse(fin, in_type), gff_out)

    if args.genome:
        with open(args.in_file) as fin, open(fasta_out, 'w') as opt_out:
            SeqIO.write(SeqIO.parse(fin, in_type), opt_out, "fasta")

    if args.prot:
        # get gene protein sequences
        all_entries = []
        with open(args.in_file) as fin, open(Pfasta_out, 'w') as prot_out:
            GBcds = SeqIO.InsdcIO.GenBankCdsFeatureIterator(fin)
            for cds in GBcds:
                if cds.seq is not None:
                    cds.id = cds.name
                    cds.description = ''
                    all_entries.append(cds)
            SeqIO.write(all_entries, prot_out, 'fasta')
        # get gene nucleotide sequences
    if args.cds:
        with open(args.in_file) as fin, open(Nfasta_out, 'w') as cds_out:
            for seq in SeqIO.parse(fin,in_type):
                for features in seq.features:
                    if features.type == 'CDS':
                        gene_info = ''
                        locus_tag = features.qualifiers['locus_tag'][0]
                        if 'protein_id' in features.qualifiers:
                            prot_id = 'protein_id=' + features.qualifiers['protein_id'][0]
                            product = features.qualifiers['product'][0]
                            gene_info = prot_id + ' ' + product
                        elif 'note' in features.qualifiers:
                            gene_info = features.qualifiers['note'][0]
                        cds = features.extract(seq.seq)
                        cds_out.write('>{} {}\n{}\n'.format(locus_tag, gene_info, cds))


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = argparse.ArgumentParser(description=__doc__, epilog=epi)
        parser.add_argument('-v', '--verbose', action='store_true',
                                default=False, help='verbose output')
        parser.add_argument('-genome',action='store_true',
                                default=False, help=('Get a FASTA file '
                                    '(default = no)'))
        parser.add_argument('-cds',action='store_true',
                                default=False, help=('Get a prot fasta file '
                                    '(default = no)'))
        parser.add_argument('-prot',action='store_true',
                                default=False, help=('Get a CDS fasta file '
                                    '(default = no)'))
        parser.add_argument('-gff',action='store_true',
                                default=False, help=('Get a gff file '
                                    '(default = no)'))
        parser.add_argument('in_file', action='store',
                                help=('Full path to the input .embl/.gbk'))
        parser.set_defaults(func=to_GFF)
        args = parser.parse_args()
        if args.verbose:
            print "Executing @ " + time.asctime()
        args.func(args)
        if args.verbose:
            print "Ended @ " + time.asctime()
            print 'Exec time minutes %f:' % ((time.time() - start_time) / 60.0)
        sys.exit(0)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)
