#!/usr/bin/env python
#coding=utf-8
import argparse
import subprocess

class ChangeSign(object):
    def __init__(self,insig_list,outsig):
        self.insig_list = insig_list
        self.outsig = outsig

    def check_outsig_outfile(self,infile,outfile):
        with open(infile,"r")as inf,open(outfile+".log","w")as otlog:
            check = "ok"
            for line in inf:
                line = line.strip("\n")
                if self.outsig in line:
                    check="no"
            if check=="no":
                otlog.write("The output delimiter you selected exists in the text and may change the original field of the file, please be aware!\n")
                otlog.write("您选择的输出分隔符在文本中存在，可能会改变文件原有字段，请知悉！\n")

    def change_main(self,infile,outfile):
        job_center = ""
        for each_insig in self.insig_list:
            job_heart = ";s/%s/%s/"%(each_insig,outsig)
            job_center += job_heart
        job_main = "'"+job_center[1:]+"g'"
        jobs = "less {infile} | sed "+job_main+" > {outfile}"
        jobs = jobs.format(infile=infile,outfile=outfile)
        subprocess.call(jobs,shell=True)
        self.check_outsig_outfile(infile,outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="qu can")
    parser.add_argument('-i','--infile',required=True, help='the input file')
    parser.add_argument('-o','--outfile',required=True, help='the output file')
    parser.add_argument('--insig',required=False,default=1, help='the signs(mybe its list combine by ";+\t") which you want to change ')
    parser.add_argument('--outsig',required=False, help='the sign(just one) which you need ,please give like "\\t" ')
    args = parser.parse_args()
    insig_list = args.insig.split("+")
    outsig = args.outsig
    mychanger = ChangeSign(insig_list,outsig)
    mychanger.change_main(args.infile,args.outfile)
