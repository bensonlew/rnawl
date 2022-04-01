#!/usr/bin/python
#-*- coding: utf-8 -*-

import re
import os
import subprocess

def gtf(infile, outfile):
    import re,subprocess,os
    awk = """awk -F "\\t" '{print $9}' %s""" %(infile)
    info = subprocess.check_output(awk, shell=True)
    data={}
    with open(outfile,'w+') as f1:
        i=0
        for line in info.strip().split("\n"):
            i+=1
            if i>=2:
                break
            else:
                try:
                    gene_id = re.search(r'gene_id\s*"(\S+)";', line).group(1)
                    transcript_id = re.search(r'transcript_id\s*"(\S+)";', line).group(1)
                except Exception:
                    print line
                    print 'haha'
                    break
                if transcript_id:
                    try:
                        if transcript_id not in data.keys():
                            if gene_id:
                                data[transcript_id] = gene_id
                                f1.write(transcript_id+"\t"+gene_id+"\n")
                            else:
                                print '{}没有找到对应的gene_id'.format(transcript_id)
                                pass
                        else:
                            pass
                    except Exception:
                        print line
                        print 'heihei'
                        break
                
    print 'end'
    
if __name__ == "__main__":
    infile = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/Vertebrates/Fish/Danio_rerio/Danio_rerio.GRCz10.87.gff3.gtf"
    outfile = "/mnt/ilustre/users/sanger-dev/workspace/20170329/Single_rsem_zebra_tools_1/Rsem/gene2transcript"
    gtf(infile, outfile)