import re
import argparse

parser = argparse.ArgumentParser(description="abstract the longest seq")
parser.add_argument("-i", "--input_fasta", help="input fa", required = True)
args = vars(parser.parse_args())

input = args["input_fasta"]
print input
i = -1
Eid = []
gene = []
seq = []
with open(input,"r") as f:
    for line in f:
        line = line.rstrip() 
        m = re.match(r">(.*) (gene=.*)$",line) 
        if m:
            Eid.append(m.group(1))
            gene.append(m.group(2))
            i += 1
            seq.append("")
        else:
            if line.find(">") != -1:
                n = re.match(r">(.*)",line)
                if n:
                    Eid.append(n.group(1))
                    gene.append(n.group(1))
                    i += 1
                    seq.append("")
            else:
                seq[i] += line 
        
if len(gene) == len(seq):
    print "yes,all fasta are prepared to go to the next step."
else:
    print "there are some errors in the program."

# print len(gene)
# print "the length of the gene list is {}".format(len(Eid))
Eid_unrepeat = []
gene_unrepeat = []
seq_unrepeat = []
cds_unrepeat = []
for i in range(len(gene)): 
    n = re.match(r"(gene=.+) (CDS=.+)",gene[i])
    if n:
        # print n.group(1)
        if n.group(1) not in gene_unrepeat:
            Eid_unrepeat.append(Eid[i])
            gene_unrepeat.append(n.group(1))
            seq_unrepeat.append(seq[i]) 
            cds_unrepeat.append(n.group(2))
        else:
            idx = gene_unrepeat.index(n.group(1))
            if len(seq[i]) > len(seq_unrepeat[idx]):
                seq_unrepeat[idx] = seq[i]
                Eid_unrepeat[idx] = Eid[i]
                cds_unrepeat[idx] = n.group(2)               
    else:
        # print gene[i]
        if gene[i] not in gene_unrepeat:
            Eid_unrepeat.append(Eid[i])
            gene_unrepeat.append(gene[i])
            seq_unrepeat.append(seq[i])
            cds_unrepeat.append("")
        else:
            idx = gene_unrepeat.index(gene[i])
            if len(seq[i]) > len(seq_unrepeat[idx]):
                seq_unrepeat[idx] = seq[i]
                Eid_unrepeat[idx] = Eid[i]
                cds_unrepeat[idx] = ""

                    
with open("the_longest_exons.fa","w") as o:
    for i in range(len(seq_unrepeat)): 
        o.write(">" + Eid_unrepeat[i] + " " + gene_unrepeat[i] + " " + cds_unrepeat[i] + "\n")
        o.write(seq_unrepeat[i] + "\n")
