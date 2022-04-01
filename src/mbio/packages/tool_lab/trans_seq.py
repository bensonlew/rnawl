#!/usr/bin/env python
import argparse
import subprocess

class TranSeq(object):
    def __init__(self,infa,outfa):
        self.infa = infa
        self.outfa = outfa

    def frame_trans(self,framenum,tablenum):
        jobs = 'transeq -sequence {infa} -outseq {outfa} -frame {framenum} -table {tablenum} '.format(infa=self.infa,outfa=self.outfa,framenum=framenum,tablenum=tablenum)
        subprocess.call(jobs,shell=True)
    
    def regions_trans(self,regionsnum_num,tablenum):
        jobs = 'transeq -sequence {infa} -outseq {outfa} -regions {regionsnum_num} -table {tablenum}'.format(infa=self.infa,outfa=self.outfa,regionsnum_num=regionsnum_num,tablenum=tablenum)
        subprocess.call(jobs,shell=True)

    def frame_trim_trans(self,framenum,tablenum):
        jobs = 'transeq -sequence {infa} -outseq {outfa} -frame {framenum} -table {tablenum} -trim '.format(infa=self.infa,outfa=self.outfa,framenum=framenum,tablenum=tablenum)
        subprocess.call(jobs,shell=True)

    def frame_clean_trans(self,framenum,tablenum):
        jobs = 'transeq -sequence {infa} -outseq {outfa} -frame {framenum} -table {tablenum} -clean '.format(infa=self.infa,outfa=self.outfa,framenum=framenum,tablenum=tablenum)
        subprocess.call(jobs,shell=True)

    def frame_trim_clean_trans(self,framenum,tablenum):
        jobs = 'transeq -sequence {infa} -outseq {outfa} -frame {framenum} -table {tablenum} -trim -clean '.format(infa=self.infa,outfa=self.outfa,framenum=framenum,tablenum=tablenum)
        subprocess.call(jobs,shell=True)

    def regions_trim_trans(self,regionsnum_num,tablenum):
        jobs = 'transeq -sequence {infa} -outseq {outfa} -regions {regionsnum_num} -table {tablenum} -trim'.format(infa=self.infa,outfa=self.outfa,regionsnum_num=regionsnum_num,tablenum=tablenum)
        subprocess.call(jobs,shell=True)

    def regions_clean_trans(self,regionsnum_num,tablenum):
        jobs = 'transeq -sequence {infa} -outseq {outfa} -regions {regionsnum_num} -table {tablenum} -clean '.format(infa=self.infa,outfa=self.outfa,regionsnum_num=regionsnum_num,tablenum=tablenum)
        subprocess.call(jobs,shell=True)

    def regions_trim_clean_trans(self,regionsnum_num,tablenum):
        jobs = 'transeq -sequence {infa} -outseq {outfa} -regions {regionsnum_num} -table {tablenum} -trim -clean '.format(infa=self.infa,outfa=self.outfa,regionsnum_num=regionsnum_num,tablenum=tablenum)
        print jobs
        subprocess.call(jobs,shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="qu can")
    parser.add_argument('--infa',required=True, help='the input file-fasta ')
    parser.add_argument('--oufa',required=True, help='the output file-fasta')
    parser.add_argument('--frame',required=False,default=1, help='Frame(s) to translate,value for 1,2,3,F(1+2+3),-1,-2,-3,R(-1and-2and-3),6(F+R)')
    parser.add_argument('--table',required=False, help='Code to use with species of your choose ,value for 0 (Standard); 1  (Standard (with alternative initiation codons)); 2 (Vertebrate Mitochondrial); 3  (Yeast Mitochondrial); 4 (Mold, Protozoan,   Coelenterate Mitochondrial and Mycoplasma/Spiroplasma); 5 (Invertebrate   Mitochondrial); 6 (Ciliate Macronuclear and  Dasycladacean); 9 (Echinoderm    Mitochondrial); 10 (Euplotid Nuclear); 11  (Bacterial); 12 (Alternative Yeast Nuclear); 13 (Ascidian Mitochondrial); 14 (Flatworm Mitochondrial); 15 (Blepharisma Macronuclear); 16 (Chlorophycean Mitochondrial); 21 (Trematode Mitochondrial); 22 (Scenedesmus obliquus); 23 (Thraustochytrium Mitochondrial))')
    parser.add_argument('--regions',required=False, help='sequence is translated for example 24-45,56-78')
    parser.add_argument('--trim',required=False, help='trim for the last X or *')
    parser.add_argument('--clean',required=False, help='This changes all STOP codon positions from the * character to X ')
    args = parser.parse_args()
###########################################
    if args.table:
        tablenum = args.table
    else:
        tablenum = 1

    myTranseq = TranSeq(args.infa,args.oufa)
    if args.frame==1 and args.regions:
        regionsnum_num = args.regions
        if args.trim and args.clean:
            myTranseq.regions_trim_clean_trans(regionsnum_num,tablenum)
        elif args.trim:
            myTranseq.regions_trim_trans(regionsnum_num,tablenum)
        elif args.clean:
            myTranseq.regions_clean_trans(regionsnum_num,tablenum)
        else:
            myTranseq.regions_trans(regionsnum_num,tablenum)
    else:
        framenum = args.frame
        if args.trim and args.clean:
            myTranseq.frame_trim_clean_trans(framenum,tablenum)
        elif args.trim:
            myTranseq.frame_trim_trans(framenum,tablenum)
        elif args.clean:
            myTranseq.frame_clean_trans(framenum,tablenum)
        else:
            myTranseq.frame_trans(framenum,tablenum)

