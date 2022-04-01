#!/usr/bin/python
# -*- coding: utf-8 -*-
# __author__: konghualei 20170418

from mako.template import Template
from biocluster.config import Config
import os

def ExpressVenn(rfile, fpkm, outfile_path):
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    print this_file_dir+"/genesetVenn.r"
    f = Template(filename = this_file_dir+"/genesetVenn.r")
    print f
    venn = f.render(fpkm = fpkm, outfile_path = outfile_path)
    with open(rfile, 'w') as rf:
        rf.write("%s"%(venn))
    
if __name__ == "__main__":
    #rfile = "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tofiles/newvenn.r"
    rfile = 'venn_new_test.r'
    # fpkm = "/mnt/ilustre/users/sanger-dev/workspace/20170602/GenesetVenn_tsg_2000_9965_3754/GenesetVenn/geneset_file_geneset_venn"
    #fpkm = "/mnt/ilustre/users/sanger-dev/workspace/20170602/GenesetVenn_tsg_2000_1077_77/geneset_file_geneset_venn"
    fpkm = "/mnt/ilustre/users/sanger-dev/workspace/20170602/GenesetVenn_tsg_2000_2210_9365/geneset_file_geneset_venn"
    #fpkm = "/mnt/ilustre/users/sanger-dev/workspace/20170413/Single_rsem_stringtie_zebra_9/Express/output/oldrsem/genes.TMM.fpkm.matrix"
    #outfile_path = "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tofiles"
    outfile_path = os.getcwd()
    ExpressVenn(rfile,fpkm,outfile_path)
