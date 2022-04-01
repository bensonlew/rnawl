import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-before",type=str,required=True,help="please input correct before stat info")
parser.add_argument("-after",type=str,required=True,help="please input correct after stat info")
parser.add_argument("-out",type=str,required=True,help="please input correct out name")
args=parser.parse_args()

before_dir=args.before
after_dir =args.after
out_dir=args.out
before_totalreads={}
before_totalBases={}
novel_info_dict={}

with open(before_dir,"r")as before_infos:
    frist_info_bebore=before_infos.readline()
    for before_info in before_infos.readlines():
        before_info=before_info.split()
        before_totalreads[before_info[0]]=before_info[1]
        before_totalBases[before_info[0]]=before_info[2]
with open(after_dir, "r")as after_infos:
    with open(out_dir,"w") as out_infos:
        frist_info_bebore = after_infos.readline()
        out_infos.write("Sample"+"\t"+"Raw_reads"+"\t"+"raw_bases"+"\t"+"Clean_reads"+"\t"+"Clean_bases"+"\t"+"Clean_reads"+"\t"+"Q20(%)"+"\t"+"Q30(%)"+"GC_content(%)"+"rRNA_Ratio(%)"+"\n")
        for after_info in after_infos.readlines():
            after_info = after_info.split()
            out_infos.write(after_info[0]+"\t"+before_totalreads[after_info[0]]+"\t"+before_totalBases[after_info[0]]+"\t"+after_info[1]+"\t"+after_info[2]+"\t"+after_info[3]+"\t"+after_info[4]+"\t"+after_info[5]+"\t"+after_info[6]+"\t"+after_info[7]+"\n")
