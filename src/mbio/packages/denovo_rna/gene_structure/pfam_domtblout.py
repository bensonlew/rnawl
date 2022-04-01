# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import re


def pfam_out(domtblout,score=False):
    with open(domtblout, "r") as f, open("pfam_domain", "w") as p:
        if score:
            p.write("Seq_id\tProtein_id\tPfam_id\tDomain\tDomainDescription\tProteinStart\tProtein\tPfamStart\tPfamEnd\tDomainE-Value\tscore\n")
            for line in f:
                if re.match(r"#", line):
                    continue
                else:
                    line = line.strip("\n").split()
                    seq_id = line[3].split("|")[0]
                    if re.match(r".*::.*::.*", line[3]):
                        seq_id = line[3].split("::")[1]
                    description = line[22:]
                    description = " ".join(description)
                    write_line = seq_id + "\t" + line[3] + "\t" + line[1] + "\t" + line[0] + "\t" + description + "\t" + line[17]
                    write_line = write_line + "\t" + line[18] + "\t" + line[15] + "\t" + line[16] + "\t" + line[11] + "\t"+ line[12] + "\n"  #增加 score
                    p.write(write_line)
        else:
            p.write("Seq_id\tProtein_id\tPfam_id\tDomain\tDomainDescription\tProteinStart\tProtein\tPfamStart\tPfamEnd\t"
                "DomainE-Value\n")
            for line in f:
                if re.match(r"#", line):
                    continue
                else:
                    line = line.strip("\n").split()
                    seq_id = line[3].split("|")[0]
                    if re.match(r".*::.*::.*", line[3]):
                        seq_id = line[3].split("::")[1]
                    description = line[22:]
                    description = " ".join(description)
                    write_line = seq_id + "\t" + line[3] + "\t" + line[1] + "\t" + line[0] + "\t" + description + "\t" + line[17]
                    write_line = write_line + "\t" + line[18] + "\t" + line[15] + "\t" + line[16] + "\t" + line[11] + "\n"
                    p.write(write_line)
