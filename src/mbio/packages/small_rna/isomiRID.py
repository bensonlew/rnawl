#!/usr/bin/python
# -*- coding: utf-8 -*-


__version__ = '0.53'
__LastModf__ = "2013, May 30"


import os
from os import makedirs
import datetime
import sys


def _cleanALLDIRS():
    top = "Results"
    try:
        for root, dirs, files in os.walk(top, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(top)
    except:
        next

def _cleanDIRS(x):
    listx = []
    if x == "r0":
        listx.append("Results/input")
    listx.append("Results/%s/SAM" % (x))
    listx.append("Results/%s/mapped" % (x))
    listx.append("Results/%s/unmapped"  % (x))
    try:
        for top in listx:
            for root, dirs, files in os.walk(top, topdown=False):
                for name in files:
                    os.remove(os.path.join(root, name))
                for name in dirs:
                    os.rmdir(os.path.join(root, name))
            os.rmdir(top)
    except:
        next

def _ParseFq(inputFile, outName):
    """
    :param inputFile:
    """
    output = open(outName, "w")
    count_seq = 0
    count = 1
    Seq = ""
    dic = dict()
    for line in open(inputFile):
        line = line.split("\n")[0]
        if count == 2:
            Seq = line.split()[0]
            if "N" in Seq:
                next;
            else:
                if Seq not in dic:
                    count_seq += 1
                    dic[Seq] = ">Seq" + str(count_seq) + ":1"
                elif Seq in dic:
                    dic[Seq] = ((dic[Seq]).split(":")[0]) + ":" + str(int((dic[Seq]).split(":")[-1]) + 1)
        if count == 4:
            count = 0
        count += 1
    for Seq in dic:
        tmp = dic[Seq] + "\n" + Seq + "\n"
        output.write(tmp)

def _ParseFa(inputFile, outName):
    """
    :param inputFile:
    """
    output = open(outName, "w")
    count_seq = 0
    dic = dict()
    for line in open(inputFile):
        line = line.split("\n")[0]
        if line.startswith(">"):
            next,
        else:
            Seq = line.split()[0]
            if "N" in Seq:
                next;
            else:
                if Seq not in dic:
                    count_seq += 1
                    dic[Seq] = ">Seq" + str(count_seq) + ":1"
                elif Seq in dic:
                    dic[Seq] = ((dic[Seq]).split(":")[0]) + ":" + str(int((dic[Seq]).split(":")[-1]) + 1)
    for Seq in dic:
        tmp = dic[Seq]+"\n"+Seq+"\n"
        output.write(tmp)

def _buildout(ref1, ref2, ref3, cab, output):
    def printlines(total, position, length, seq):
        first = position - 1
        end = total-first-length
        line = ("."*(position-1))+seq+("."*end)
        return line
    out = open(output, "w")
    for r in ref2:
        refName = r
        refSeq = ref1[r]["refSeq"]
        refLen = ref1[r]["refLen"]
        line1 = refName + ("-"*(refLen-len(refName))) +"\tRound\tLength\tVar\tPosition" + cab + "\n"
        line2 = ""
        line3 = refSeq + "\t" + refName + "\t" + str(refLen) + "\n"
        line4 = ""
        try:
            for m in ref3[r]:
                matName = m
                matLen = ref3[r][m]["matLen"]
                matSeq = ref3[r][m]["matSeq"]
                matPos = ref3[r][m]["matPos"]
                line2 += printlines(refLen, matPos, matLen, matSeq) + "\t" + matName + "\t" + str(matLen) +"\n"
        except:
            next;
        try:
            tmp1 = dict()
            tmp2 = dict()
            for t in ref2[r]:
                if ref2[r][t]["Position"] > 0:
                    tmp1[t] = ref2[r][t]["Position"]
                elif ref2[r][t]["Position"] < 0:
                    tmp2[t] = ref2[r][t]["Position"]
            list1 = sorted(tmp1.items(),lambda x,y: cmp(y[1],x[1]))
            list2 = sorted(tmp2.items(),lambda x,y: cmp(y[1],x[1]))
            for s1 in list1:
                s = s1[0]
                Modification = (ref2[r][s]["Modification"]).split(",")[0]
                if Modification != "no":
                    if ">" in Modification:
                        Position = ref2[r][s]["Position"]
                    else:
                        Position = ref2[r][s]["Position"]-len(Modification)
                else:
                    Position = ref2[r][s]["Position"]
                Freqs = ref2[r][s]["Freqs"]
                ModSeq = ref2[r][s]["ModSeq"]
                Mod = ref2[r][s]["Mod"]
                Rn = ref2[r][s]["Round"]
                Len = len(ModSeq)
                line1 += printlines(refLen, Position, Len, ModSeq) +"\t"+ Rn +"-"+ Mod + "\t" + str(Len) + "\t" + Modification + "\t" + str(Position) +"\t"+ Freqs + "\n"
            for s1 in list2:
                s = s1[0]
                Modification = (ref2[r][s]["Modification"]).split(",")[0]
                if Modification != "no":
                    if ">" in Modification:
                        Position = ref2[r][s]["Position"]*-1
                    else:
                        Position = (ref2[r][s]["Position"]-len(Modification))*-1
                else:
                    Position = ref2[r][s]["Position"]*-1
                Freqs = ref2[r][s]["Freqs"]
                ModSeq = _RevComp(ref2[r][s]["ModSeq"])
                Mod = ref2[r][s]["Mod"]
                Rn = ref2[r][s]["Round"]
                Len = len(ModSeq)
                line4 += printlines(refLen, Position, Len, ModSeq) +"\t"+ Rn + "-" + Mod + "\t" + str(Len) + "\t" + Modification + "\t" + str(Position*-1) + "\t"+ Freqs + "\n"

        except:
            next;
        lines = line1 +  line2 + line3 + line4 + "\n\n"
        out.write(lines)

def __mature(input1, sp, mir):
        cntrl = 0
        dic = dict()
        name = ""
        seq = ""
        rootname = len(sp)+len(mir)+1
        for line in open(input1):
            if line.startswith(">"):
                name = (line.split()[0])[1:]
                if name.split("-")[-1] == "3p" or name.split("-")[-1] == "5p":
                    Rname = sp + "-" + mir + name[rootname:-3]
                else:
                    Rname = sp + "-" + mir + name[rootname:]
                seq = ""
            else:
                seq += line.split()[0]
                if Rname not in dic:
                    dic[Rname] = dict()
                    dic[Rname][name] = dict()
                    dic[Rname][name]["matName"] = name
                    dic[Rname][name]["matSeq"] = seq
                    dic[Rname][name]["matLen"] = len(seq)
                if Rname in dic:
                    if name not in dic[Rname]:
                        dic[Rname][name] = dict()
                        dic[Rname][name]["matName"] = name
                        dic[Rname][name]["matSeq"] = seq
                        dic[Rname][name]["matLen"] = len(seq)
        return dic

def _mature():
        input1 = "Results/r0/SAM/mature-clean.sam"
        dic = dict()
        seq = ""
        name = ""
        for line1 in open(input1):
            if line1.startswith("@"):
                next
            else:
                tmpline = line1.split("\t")
                MIRname = tmpline[0]
                Sense = tmpline[1]
                Rname = tmpline[2]
                Position = int(tmpline[3])
                seq = tmpline[9]

                if Rname not in dic:
                    dic[Rname] = dict()
                    dic[Rname][MIRname] = dict()
                    dic[Rname][MIRname]["matName"] = MIRname
                    dic[Rname][MIRname]["matSeq"] = seq
                    dic[Rname][MIRname]["matLen"] = len(seq)
                    dic[Rname][MIRname]["matPos"] = Position
                if Rname in dic:
                    if name not in dic[Rname]:
                        dic[Rname][MIRname] = dict()
                        dic[Rname][MIRname]["matName"] = MIRname
                        dic[Rname][MIRname]["matSeq"] = seq
                        dic[Rname][MIRname]["matLen"] = len(seq)
                        dic[Rname][MIRname]["matPos"] = Position
        return dic

def _GetRef(input1):
    dic = dict()
    Rname = ""
    for line in open(input1):
            if line.startswith(">"):
                Rname = (line.split()[0])[1:]
            else:
                seq = line.split()[0]
                dic[Rname] = dict()
                dic[Rname]["refSeq"] = seq
                dic[Rname]["refLen"] = len(seq)
    return dic

def _ParseRn(inputFile, dic):
    count = 0
    cablibs = ""
    def _ModifySeq(Seq, mod):
        newSeq = ""
        if mod != "no":
            position = int(mod.split(",")[1])-1
            for i in range(0, len(Seq)):
                if i == position:
                    newSeq += Seq[i].lower()
                else:
                    newSeq += Seq[i]
        return newSeq
    for line in open(inputFile):
        refs = dict()
        freqs = ""
        col = len(line.split("\t"))
        if line.startswith("OriginalSeq"):
            lines = line.split("\t")
            start = 5
            while start < col:
                cablibs += "\t" + (lines[start]).split()[0]
                start += 1
        else:
            lines = line.split("\t")
            seq = lines[0]
            rn = lines[1]
            seqLen = lines[2]
            start = 4
            while start < col:
                freqs += "\t" + (lines[start]).split()[0]
                start += 1
            ref = lines[3][1:-1].split(", ")
            for i in ref:
                refs[(i[1:-1]).split(":")[0]] = dict()
                refs[(i[1:-1]).split(":")[0]]["Position"] = int((i[1:-1]).split(":")[1])
                refs[(i[1:-1]).split(":")[0]]["Modification"] = (i[1:-1]).split(":")[2]
                refs[(i[1:-1]).split(":")[0]]["Mod"] = (i[1:-1]).split(":")[3]

        for r in refs:
            count += 1
            SeqID = "Seq"+"-"+str(rn)+str(count)
            if r not in dic:
                dic[r] = dict()
                dic[r][SeqID] = dict()
                dic[r][SeqID]["Modification"] = refs[r]["Modification"]
                dic[r][SeqID]["Mod"] = refs[r]["Mod"]
                if dic[r][SeqID]["Modification"] != "no":
                    dic[r][SeqID]["ModSeq"] = _ModifySeq(seq, dic[r][SeqID]["Modification"])
                else:
                    dic[r][SeqID]["ModSeq"] = seq
                dic[r][SeqID]["Position"] = refs[r]["Position"]
                dic[r][SeqID]["Freqs"] = freqs
                dic[r][SeqID]["Round"] = rn


            if r in dic:
                dic[r][SeqID] = dict()
                dic[r][SeqID]["Modification"] = refs[r]["Modification"]
                dic[r][SeqID]["Mod"] = refs[r]["Mod"]
                if dic[r][SeqID]["Modification"] != "no":
                    dic[r][SeqID]["ModSeq"] = _ModifySeq(seq, dic[r][SeqID]["Modification"])
                else:
                    dic[r][SeqID]["ModSeq"] = seq
                dic[r][SeqID]["Position"] = refs[r]["Position"]
                dic[r][SeqID]["Freqs"] = freqs
                dic[r][SeqID]["Round"] = rn
    return dic, cablibs

class Reference:
    def __init__(self, name, refpath, buildindex, bowtiebuild):
        self.Name = name
        self.Refpath = refpath
        self.Buildindex = buildindex
        self.BowtieBuild = bowtiebuild

    def _BuidIndex(self):
        if self.Buildindex == "yes":
            print "Building indexes for %s ..." %(self.Name)
            print self.Refpath
            os.system("%s %s %s" %(self.BowtieBuild, self.Refpath, self.Refpath))
            print "... The indexes was built\n"
        else:
            print "... It was not necessary to build %s indexes\n" %(self.Name)

class Library():
    def __init__(self, libname, libpath, format):
        self.LibName = libname
        self.LibPath = libpath
        self.Format = format

    def _getLibInfo(self):
        print self.LibName, self.LibPath, self.Format

class SAMlist():
    def __init__(self, rounds, libs, add):
        self.rounds = rounds + 1
        self.libs = libs
        if add == 0:
            self.add = "r0"
        if add == 1:
            self.add = "r1"
        elif add > 1:
            self.add = "M" + str(add)
        self.samlist = dict()

    def _getSAMlist(self):
        for i in range(self.rounds):
            if i == 0:
                if self.add == "r0":
                    name = "r0"
                    self.samlist[name] = dict()
                    for l in self.libs:
                        self.samlist[name][l] = "Results/%s/SAM/r%i-%s-clean.sam" %(self.add, i, l)
                else:
                    next
            if i == 1:
                if self.add == "r1":
                    name = "r1"
                    self.samlist[name] = dict()
                    for l in self.libs:
                        self.samlist[name][l] = "Results/%s/SAM/r%i-%s-clean.sam" %(self.add, i, l)
                        print self.samlist[name][l]
                else:
                    next
            elif i > 1:
                name = "r"+str(i)
                self.samlist[name] = dict()
                for l in self.libs:
                    self.samlist[name][l] = "Results/%s/SAM/r%i-%s-clean.sam" %(self.add, i, l)
        return self.samlist

def _3sRNAaddfa(readFile, output, z):
    z = int(z)
    output1 = open(output, "w")
    yeseq = r'A|a|T|t|C|c|G|g'
    pegoFa = 0
    for line in open(readFile):
        line = line.split("\n")[0]
        if line.startswith(">"):
            if pegoFa == 0:
                seqName = line+"&"
                if (seqName.split("&")[1]).split("_")[0] == "SIZE":
                    name = seqName.split("&")[0]
                    size = (seqName.split("&")[1]).split("_")[1]
                    addx = (seqName.split("&")[2]).split("_")[1]
                else:
                    name = seqName.split("&")[0]
                    size = "x"
                    addx = ""
            elif pegoFa == 1:
                complete = name + "&SIZE_" + size + "&M3_" + add + addx + "\n" + seq + "\n"
                output1.write(complete)
                seqName = line+"&"
                if (seqName.split("&")[1]).split("_")[0] == "SIZE":
                    name = seqName.split("&")[0]
                    size = (seqName.split("&")[1]).split("_")[1]
                    addx = (seqName.split("&")[2]).split("_")[1]
                else:
                    name = seqName.split("&")[0]
                    size = "x"
                    addx = ""
                pegoFa = 0
        else:
            seq1 = line.split()[0]
            if z == 2:
                add = seq1[-2:]
                seq = seq1[:-2]
            elif z == 1:
                add = seq1[-1]
                seq = seq1[:-1]
            pegoFa = 1
    complete = name + "&SIZE_" + size + "&M3_" + add + addx + "\n" + seq + "\n"
    output1.write(complete)

def _5sRNAaddfa(readFile, output, z):
    output1 = open(output, "w")
    yeseq = r'A|a|T|t|C|c|G|g'
    pegoFa = 0
    for line in open(readFile):
        line = line.split("\n")[0]
        if line.startswith(">"):
            if pegoFa == 0:
                seqName = line+"&"
                if (seqName.split("&")[1]).split("_")[0] == "SIZE":
                    name = seqName.split("&")[0]
                    size = (seqName.split("&")[1]).split("_")[1]
                    addx = (seqName.split("&")[2]).split("_")[1]
                else:
                    name = seqName.split("&")[0]
                    size = "x"
                    addx = ""
            elif pegoFa == 1:
                complete = name + "&SIZE_" + size + "&M5_" + addx + add + "\n" + seq + "\n"
                output1.write(complete)
                seqName = line+"&"
                if (seqName.split("&")[1]).split("_")[0] == "SIZE":
                    name = seqName.split("&")[0]
                    size = (seqName.split("&")[1]).split("_")[1]
                    addx = (seqName.split("&")[2]).split("_")[1]
                else:
                    name = seqName.split("&")[0]
                    size = "x"
                    addx = ""
                pegoFa = 0
        else:
            seq1 = line.split()[0]
            add = seq1[:z]
            seq = seq1[z:]
            pegoFa = 1
    complete = name + "&SIZE_" + size + "&M5_" + addx + add + "\n" + seq + "\n"
    output1.write(complete)

def _appCutoff(input, name, cutoff, dir):
    inputA = open(input)
    cut = int(cutoff)
    out = dir+"Cutoff/"+name+"-Cutoff"+cutoff+".txt"
    output = open(out, "w")
    for line in inputA:
        if line.startswith("OriginalSeq"):
            output.write(line)
        else:
            sum = 0
            line1 = line.split("\t")[4:]
            n = len(line1)
            for i in line1:
                sum += int(i)
            if sum >= (cut):
                output.write(line)
    return out

def _RevComp(sequence):
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A','a':'t', 'c':'g', 'g':'c', 't':'a'}
    return "".join([complement.get(nt, 'N') for nt in sequence[::-1]])

def _Rev(sequence):
    return sequence[::-1]

def _Comp(sequence):
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A','a':'t', 'c':'g', 'g':'c', 't':'a'}
    return "".join([complement.get(nt, 'N') for nt in sequence])

def _freqSAM(input1, n, R, output, nadd):
    yesSeqNT = r'A|a|T|t|C|c|G|g'
    out = open(output, 'w')
    s = 5
    nlib = 4
    dic = dict()
    n = int(n)
    x = n+s
    cab = ""
    Seq = ""
    mod = ""
    startall = datetime.datetime.now()
    for line in input1:
        start = datetime.datetime.now()
        libName = line
        print "%s..."%(libName)
        lib2 = input1[line]
        nlib += 1
        cab += "\t" + libName
        samFile = open(lib2)
        for line1 in samFile.readlines():
            mod = "no"
            if line1.startswith("@"):
                next
            else:
                tmpName = line1.split("\t")[0]
                tmpName = tmpName+"&"
                if (tmpName.split("&")[1]).split("_")[0] == "SIZE":
                    rawname = (tmpName.split("&")[0]).split(":")[0]
                    Freq = int((tmpName.split("&")[0]).split(":")[1])
                    Size = (tmpName.split("&")[1]).split("_")[1]
                    Add = (tmpName.split("&")[2]).split("_")[1]
                    mod = Add
                    Add1 = (tmpName.split("&")[2]).split("_")[0]
                else:
                    rawname = (tmpName.split("&")[0]).split(":")[0]
                    Freq = int((tmpName.split("&")[0]).split(":")[1])
                    Size = ""
                    Add = ""
                    Add1 = "no"
                Test_Sense = line1.split("\t")[1]
                if Test_Sense == "0":
                    Sense = ""
                    Seqmap = line1.split("\t")[9]
                else:
                    Seqmap = _RevComp(line1.split("\t")[9])
                    Sense = "-"
                if nadd == 3:
                    Seq = Seqmap+Add
                    mod += ","+str(len(Seq)) + ":M3"
                if nadd == 5:
                    Seq = Add+Seqmap
                    mod += ",1" + ":M5"
                if nadd == 0:
                    Seq = Seqmap
                    mod += ":r0"
                if nadd == 1:
                    Seq = Seqmap
                    mod += ":r1"

                mm = (line1.split("\t")[12]).split(":")[-1]
                SNP = ""
                if mm == str(len(Seqmap)): # MD:Z:22
                    next,
                if mm != str(len(Seq)): # MD:Z:0C21 | MD:Z:21A0 | MD:Z:16A5 | MD:Z:4C17
                    try:
                        if mm[1]=="A" or mm[1]=="T" or mm[1]=="C" or mm[1]=="G":
                            SNP = Seq[int(mm[0])]
                            mod = mm[1] + ">" + SNP + "," + str(int(mm[0])+1)
                            if mm[0] == "0":
                                Add1 = "M5"
                                mod += ":" + Add1
                            else:
                                Add1 = "MM"
                                mod += ":" + Add1
                        elif mm[2]=="A" or mm[2]=="T" or mm[2]=="C" or mm[2]=="G":
                            SNP = Seq[int(mm[0:2])]
                            mod = mm[2] + ">" + SNP + "," + str(int(mm[0:2])+1)
                            if mm[0:2] == str(len(Seqmap)-1):
                                Add1 = "M3"
                                mod += ":" + Add1
                            else:
                                Add1 = "MM"
                                mod += ":" + Add1
                    except:
                        next;
                Size = str(len(Seq))
                RefSeq = line1.split("\t")[2]+":"+Sense+line1.split("\t")[3]+":"+mod
                if SNP != "N":
                    if Seq in dic:
                        dic[Seq][nlib] = Freq
                        if rawname not in dic[Seq][0]:
                            dic[Seq][0].append(rawname)
                            if RefSeq not in dic[Seq][4]:
                                dic[Seq][4].append(RefSeq)
                        elif rawname in dic[Seq][0]:
                            if RefSeq not in dic[Seq][4]:
                                dic[Seq][4].append(RefSeq)
                    elif Seq not in dic:
                        dic[Seq] = []
                        for i in range(s):
                            if i == 0 or i == 4:
                                dic[Seq].append([])
                            else:
                                dic[Seq].append("")
                        for i in range(n):
                            dic[Seq].append(0)
                        dic[Seq][0].append(rawname)
                        dic[Seq][1] = Seq
                        dic[Seq][2] = R
                        dic[Seq][3] = Size
                        dic[Seq][4].append(RefSeq)
                        dic[Seq][nlib] = Freq
        end = datetime.datetime.now() - start
        print "\t%s is done\t(time: %s)"%(libName,end)
    cab = "OriginalSeq"+"\t"+"Round"+"\t"+"Variation"+"\t"+"Length"+"\t"+"preMIRref"+ cab + "\n"
    out.write(cab)
    print "\tcreating Frequences Tab"
    for Seq in dic:
        tab = ""
        for i in range(x):
            if i != 0:
                if i < (x-1):
                    if i == 4:
                        dic[Seq][4].sort()
                        ",".join(dic[Seq][4])
                    tmp1 = dic[Seq][i]
                    tab += str(tmp1)+"\t"
                else:
                    tmp1 = dic[Seq][i]
                    tab += str(tmp1)
        tmp2 = tab + "\n"
        out.write(tmp2)
    total = datetime.datetime.now() - startall
    print "\ttotal time read libs: %s"%(total)

def _cleanSAM(inputA, output):
    out = open(output, "wb")
    cab = ""
    dic = dict()
    for line in open(inputA):
        if line.startswith("@"):
            out.write(line)
            if line.startswith("@SQ"):
                locus = (line.split("@SQ\tSN:")[1]).split("\t")[0]
                dic[locus] = [1]
        else:
            seq = line.split("\t")[2]
            if seq in dic:
                out.write(line)
            else:
                next

def _Round_mir(bowtie, refpath, mir_ref, root0):
    print "\n*** Mapping mature miRNA on precursor"
    inp = "f"
    os.system("%s --norc -n 0 -v 0 -a -l 6 -t %s -%s %s --sam %sSAM/mature.sam" % (bowtie, refpath, inp, mir_ref, root0))
    inputA = "%sSAM/mature.sam" %(root0)
    inputB = "%sSAM/mature-clean.sam" %(root0)
    _cleanSAM(inputA, inputB)
    os.remove("%sSAM/mature.sam" %(root0))

def _Round0(bowtie, refpath, libname, libpath, format, filter, filterpath, root0):

    """

    :param bowtie:
    :param refpath:
    :param libname:
    :param libpath:
    :param format:
    :param filter:
    :param filterpath:
    :param root0:
    """
    print "\n*** Mapping Round 0...%s" %(libname)
    if format == "fa":
        inp = "f"
    else:
        inp = "q"
    if filter == "yes":
        os.system("%s --norc -n 0 -v 0 -a -l 6 -t %s -%s %s --sam %sSAM/r0-%s.sam --un %sunmapped/unr0tmp-%s.%s --al %smapped/r0-%s.%s" %(bowtie, refpath, inp, libpath, root0, libname, root0, libname, format, root0, libname, format))
        os.system("%s --norc -n 0 -v 0 -a -l 6 -t %s -%s %sunmapped/unr0tmp-%s.%s --sam %sSAM/r0tmp-%s.sam --un %sunmapped/unr0-%s.%s --al %smapped/r0-%s.%s" %(bowtie, filterpath, inp, root0, libname, format, root0, libname, root0, libname, format, root0, libname, format))
        os.remove("%sunmapped/unr0tmp-%s.%s" %(root0, libname, format))
        os.remove("%sSAM/r0tmp-%s.sam" %(root0, libname))
    elif filter == "no":
        os.system("%s --norc -n 0 -v 0 -a -l 6 -t %s -%s %s --sam %sSAM/r0-%s.sam --un %sunmapped/unr0-%s.%s --al %smapped/r0-%s.%s" %(bowtie, refpath, inp, libpath, root0, libname, root0, libname, format, root0, libname, format))
    inputA = "%sSAM/r0-%s.sam" %(root0, libname)
    inputB = "%sSAM/r0-%s-clean.sam" %(root0, libname)
    _cleanSAM(inputA, inputB)
    os.remove("%sSAM/r0-%s.sam" %(root0, libname))

def _Round1(bowtie, refpath, libname, libpath, format, n, root1, root0):
    add = "M"
    inp = "f"
    print "*** Mapping %sadd Round %s...%s" %(add, n, libname)
    os.system("%s --norc -n 1 -v 1 -a -l 6 -t %s -%s %sunmapped/unr0-%s.%s --sam %sSAM/r1-%s.sam --un %sunmapped/unr1-%s.%s --al %smapped/r1-%s.%s" %(bowtie, refpath, inp, root0, libname, format, root1, libname, root1, libname, format, root1, libname, format))
    inputA = "%sSAM/r1-%s.sam" %(root1, libname)
    inputB = "%sSAM/r1-%s-clean.sam" %(root1, libname)
    _cleanSAM(inputA, inputB)
    os.remove("%sSAM/r1-%s.sam" %(root1, libname))

def _RoundN(add, bowtie, refpath, libname, libpath, format, n, rootN, root1):
    y=n-1
    if format == "fa":
        inp = "f"
    else:
        inp = "q"
    if n == 2:
        print "*** Mapping M%s Round %s...%s" %(add, n, libname)
        inputA = "%sunmapped/unr1-%s.%s" %(root1, libname, format)
        inputB = "%sunmapped/un%sr1-%s.%s" %(root1, add, libname, format)
        print "\t*** Creating un%sr1-%s.%s" %(add, libname, format)
        if add == "3":
            _3sRNAaddfa(inputA, inputB, 2)
        if add == "5":
            _5sRNAaddfa(inputA, inputB, 2)
        os.system("%s --norc -n 0 -v 0 -a -l 6 -t %s -%s %sunmapped/un%sr1-%s.%s --sam %sSAM/r2-%s.sam --un %sunmapped/unr2-%s.%s --al %smapped/r2-%s.%s" %(bowtie, refpath, inp, root1, add, libname, format, rootN, libname, rootN, libname, format, rootN, libname, format))
        inputA = "%sSAM/r%s-%s.sam" %(rootN, n, libname)
        inputB = "%sSAM/r%s-%s-clean.sam" %(rootN, n, libname)
        _cleanSAM(inputA, inputB)
        os.remove("%sSAM/r%s-%s.sam" %(rootN, n, libname))

    elif n > 2:
        print "\n*** Mapping M%s Round %s...%s" %(add, n, libname)
        inputA = "%sunmapped/unr%s-%s.%s" %(rootN, y, libname, format)
        inputB = "%sunmapped/un%sr%s-%s.%s" %(rootN, add, y, libname, format)
        print "\t*** Creating un%sr%s-%s.%s" %(add, y, libname, format)
        if add == "3":
            if format == "fq":
                _3sRNAaddfq(inputA, inputB, 1)
            elif format == "fa":
                _3sRNAaddfa(inputA, inputB, 1)
        if add == "5":
            if format == "fq":
                _5sRNAaddfq(inputA, inputB, 1)
            elif format == "fa":
                _5sRNAaddfa(inputA, inputB, 1)
        os.system("%s --norc -n 0 -v 0 -a -l 6 -t %s -%s %sunmapped/un%sr%s-%s.%s --sam %sSAM/r%s-%s.sam --un %sunmapped/unr%s-%s.%s --al %smapped/r%s-%s.%s" %(bowtie, refpath, inp, rootN, add, y, libname, format, rootN, n, libname, rootN, n, libname, format, rootN, n, libname, format))
        inputA = "%sSAM/r%s-%s.sam" %(rootN, n, libname)
        inputB = "%sSAM/r%s-%s-clean.sam" %(rootN, n, libname)
        _cleanSAM(inputA, inputB)
        os.remove("%sSAM/r%s-%s.sam" %(rootN, n, libname))

def _ReadConfig():
    config = open(sys.argv[1])
    bowtie = []
    libs = dict()
    mainref = []
    filterref = []
    mirnaref = []
    add3 = []
    add5 = []
    Size_cutoff = []
    cutoffs = []
    for line in config.readlines():
        line = line.strip()
        if line.startswith("#"):
            next;
        if line.startswith("bowtie_path:"): #bowtie
            bowtie_path = line.split()[1]
            bowtie.append(bowtie_path)
        if line.startswith("bowtie-build_path:"):
            bowtiebuild_path = line.split()[1]
            bowtie.append(bowtiebuild_path)
        if line.startswith("lib:"): #libs
            libName = line.split()[2]
            libPath = line.split()[1]
            format = line.split()[3]
            libs[libName]=[]
            libs[libName].append("")
            libs[libName].append("")
            libs[libName][0]= libPath
            libs[libName][1]= format
        if line.startswith("main_ref:"): #main ref
            mainref_path = line.split()[1]
            mainref_name = mainref_path.split("/")[-1]
            mainref_index = line.split()[2]
            mainref.append(mainref_name)
            mainref.append(mainref_path)
            mainref.append(mainref_index)
        if line.startswith("filter_ref:"): #filter ref
            filterref_path = line.split()[1]
            if filterref_path == "no":
                filterref = filterref_path
            else:
                filterref_name = os.path.split(filterref_path)[-1]
                #filterref_name = filterref_path.split("/")[-1]
                filterref_index = line.split()[2]
                filterref.append(filterref_name)
                filterref.append(filterref_path)
                filterref.append(filterref_index)
        if line.startswith("known_miRNAs:"): #known miRNA
            miRNAref_path = line.split()[1]
            if miRNAref_path == "no":
                mirnaref.append(miRNAref_path)
            else:
                miRNAref_name = os.path.split(miRNAref_path)[-1]
                #miRNAref_name = miRNAref_path.split("/")[-1]
                mirnaref.append(miRNAref_name)
                mirnaref.append(miRNAref_path)
        if line.startswith("M3:"): #3Add and Rounds
            add3.append((line.split()[1]))
            add3.append((line.split()[2]))
        if line.startswith("M5:"): #5Add and Rounds
            add5.append((line.split()[1]))
            add5.append((line.split()[2]))
        if line.startswith("RangeSize:"): #Size cutoff
            min = int(line.split()[1])
            max = int(line.split()[2])
            Size_cutoff.append(min)
            Size_cutoff.append(max)
        if line.startswith("cutoff:"): #Frequence cutoff
            cutoff = line.split()[1]
            cutoffs.append(cutoff)
        else:
            next;
    return bowtie, libs, mainref, filterref, mirnaref, add3, add5, Size_cutoff, cutoffs

def _Complete_analysis():
    bowtie, libs, mainref, filterref, mirnaref, add3, add5, Size_cutoff, cutoffs = _ReadConfig()

    ### Creating Folders
    #For r0
    try:
        makedirs("Results/input")
        makedirs("Results/r0/SAM")
        makedirs("Results/r0/mapped")
        makedirs("Results/r0/unmapped")
        makedirs("Results/r0/RawResults")
        makedirs("Results/r0/Cutoff")

        makedirs("Results/MapResults")

        makedirs("Results/r1/SAM")
        makedirs("Results/r1/mapped")
        makedirs("Results/r1/unmapped")
        makedirs("Results/r1/RawResults")
        makedirs("Results/r1/Cutoff")

        # For add3
        add3_ = "no"
        if add3[0] == "yes":
            round3 = int(add3[1])+1
            add3_ = add3[0]
            makedirs("Results/M3/SAM")
            makedirs("Results/M3/mapped")
            makedirs("Results/M3/unmapped")
            makedirs("Results/M3/RawResults")
            makedirs("Results/M3/Cutoff")

        # For M5
        add5_ = "no"
        if add5[0] == "yes":
            round5 = int(add5[1])+1
            add5_ = add5[0]
            makedirs("Results/M5/SAM")
            makedirs("Results/M5/mapped")
            makedirs("Results/M5/unmapped")
            makedirs("Results/M5/RawResults")
            makedirs("Results/M5/Cutoff")

    except:
        print "\nWARNING!\n\tThe Results directory already exists\n\n\tPlease, rename this old directory or include  -f or --force parameter on terminal in order to erase this directory.\n\n\tEx: python sRNAadd.py --force\n\n"
        exit(1)

    # set bowtie
    bowtie_path = bowtie[0]
    bowtiebuild_path = bowtie[1]
    # set main ref
    mainref_name = mainref[0]
    mainref_path = mainref[1]
    mainref_index = mainref[2]
    main_ref = Reference(mainref_name, mainref_path, mainref_index, bowtiebuild_path)
    main_ref._BuidIndex()
    # set filter ref
    if filterref != "no":
        filter = "yes"
        filterref_name = filterref[0]
        filterref_path = filterref[1]
        filterref_index = filterref[2]
        filter_ref = Reference(filterref_name, filterref_path, filterref_index, bowtiebuild_path)
        filter_ref._BuidIndex()
    else:
        filter = "no"
        filterref_path = ""

        next;
        # set Size Cutoff
    min = Size_cutoff[0]
    max = Size_cutoff[1]

    ### Set other variables:
    mature_dic = ""

    ### Perform Parse and mapping
    for l in libs:
        libname = l
        libpath_tmp = libs[l][0]
        format = libs[l][1]

        #Parse input
        libpath = "Results/input/"+ libname +"-nonRed.fa"
        print "\ncreating %s..." % (libpath)
        if libs[l][1] == "fq":
            _ParseFq(libpath_tmp, libpath)
            format = "fa"
        if libs[l][1] == "fa":
            _ParseFa(libpath_tmp, libpath)
            format = "fa"

        # For R0
        root0 = "Results/r0/"
        print "\n*** Performing mapping R0"
        _Round0(bowtie_path, mainref_path, libname, libpath, format, filter, filterref_path, root0)

        # For R1
        root1 = "Results/r1/"
        print "\n*** Performing mapping R1"
        _Round1(bowtie_path, mainref_path, libname, libpath, format, 1, root1, root0)

        # For ADD3
        if add3_ == "yes":
            root3 = "Results/M3/"
            add3 = "3"
            print "\n*** Performing mapping M3"
            for i in range(round3):
                n = i + 1
                if n == 0:
                    next,
                if n == 1:
                    next,
                else:
                    _RoundN(add3, bowtie_path, mainref_path, libname, libpath, format, n, root3, root1)

        # For ADD5
        if add5_ == "yes":
            root5 = "Results/M5/"
            add5 = "5"
            print "\n*** Performing mapping M5"
            for i in range(round5):
                n = i+1
                if n == 0:
                    next
                if n == 1:
                    next
                else:
                    _RoundN(add5, bowtie_path, mainref_path, libname, libpath, format, n, root5, root1)

    #Perform freq tab r0:
    print "\n*** Getting frequences for r0"
    outListr0 = []
    rnList0 = []
    roundr0 = 0
    SAMr0 = SAMlist(roundr0, libs, 0)
    SAMlistr0 = SAMr0._getSAMlist()
    for rn in SAMlistr0:
        rnList0.append(rn)
        samlistr0 = SAMlistr0[rn]
        nlibs = len(SAMlistr0[rn])
        out = "Results/r0/RawResults/"+rn+".txt"
        outListr0.append(out)
        header0 = _freqSAM(samlistr0, nlibs, rn, out, 0)

        #Applying cutoffs R0
        cutoffList = []
        cutoffd = dict()
        for co in cutoffs:
            if co == "no":
                print "\tThere is no cutoff to applying for R0"
            else:
                print "Applying frequencing cutoff of %s reads:" % (co)
                cutoffd[co] = dict()
                r = "r0"
                for r in rnList0:
                    cutoffd[co][r] = dict()
                    print "\tFor %s..." % (r)
                    input = "Results/r0/RawResults/"+r+".txt"
                    cutofftemp = _appCutoff(input, r, co, root0)
                    cutoffd[co][r]["r0"] = cutofftemp
                    cutoffList.append(cutofftemp)
    #Perform freq tab r1:
    print "\n*** Getting frequences for r1"
    outListr1 = []
    rnList1 = []
    roundr1 = 1
    SAMr1 = SAMlist(roundr1, libs, 1)
    SAMlistr1 = SAMr1._getSAMlist()
    Tablist = []
    for rn in SAMlistr1:
        rnList1.append(rn)
        samlistr1 = SAMlistr1[rn]
        nlibs = len(SAMlistr1[rn])
        out = "Results/r1/RawResults/"+rn+".txt"
        outListr1.append(out)
        header0 = _freqSAM(samlistr1, nlibs, rn, out, 1)

        #Applying cutoffs R0
        cutoffList = []
        for co in cutoffs:
            if co == "no":
                print "\tThere is no cutoff to applying for R0"
            else:
                print "Applying frequencing cutoff of %s reads:" % (co)
                r = "r0"
                for r in rnList1:
                    cutoffd[co][r] = dict()
                    print "\tFor %s..." % (r)
                    input = "Results/r1/RawResults/"+r+".txt"
                    cutofftemp = _appCutoff(input, r, co, root1)
                    cutoffd[co][r]["r1"] = cutofftemp
                    cutoffList.append(cutofftemp)

    #Perform freq tab add3:
    if add3_ == "yes":
        print "\n*** Getting frequences for M3"
        outListM3= []
        rnList3 = []
        SAM3 = SAMlist(round3, libs, 3)
        SAMlist3 = SAM3._getSAMlist()
        for rn in SAMlist3:
            rnList3.append(rn)
            samlist3 = SAMlist3[rn]
            nlibs = len(SAMlist3[rn])
            out = "Results/M3/RawResults/"+rn+".txt"
            outListM3.append(out)
            header3add = _freqSAM(samlist3, nlibs, rn, out, 3)

        #Applying cutoffs add3
        cutoffList = []
        for co in cutoffs:
            if co == "no":
                print "\tThere is no cutoff to applying for M3"
            else:
                print "Applying frequencing cutoff of %s reads:" % (co)
                r = "r0"
                for r in rnList3:
                    if r in cutoffd[co]:
                        next,
                    else:
                        cutoffd[co][r] = dict()
                    print "\tFor %s..." % (r)
                    input = "Results/M3/RawResults/"+r+".txt"
                    cutofftemp = _appCutoff(input, r, co, root3)
                    cutoffd[co][r]["M3"] = cutofftemp
                    cutoffList.append(cutofftemp)
    #Perform freq tab add5:
    if add5_ == "yes":
        print "\n*** Getting frequences for M5"
        outListM5= []
        rnList5 = []
        SAM5 = SAMlist(round5, libs, 5)
        SAMlist5 = SAM5._getSAMlist()
        for rn in SAMlist5:
            rnList5.append(rn)
            samlist5 = SAMlist5[rn]
            nlibs = len(SAMlist5[rn])
            out = "Results/M5/RawResults/"+rn+".txt"
            outListM5.append(out)
            header5add = _freqSAM(samlist5, nlibs, rn, out, 5)

        #Applying cutoffs add5
        cutoffList = []
        for co in cutoffs:
            if co == "no":
                print "\tThere is no cutoff to applying for M5"
            else:
                print "Applying frequencing cutoff of %s reads:" % (co)
                r = "r0"
                for r in rnList5:
                    if r in cutoffd[co]:
                        next,
                    else:
                        cutoffd[co][r] = dict()
                    print "\tFor %s..." % (r)
                    input = "Results/M5/RawResults/"+r+".txt"
                    cutofftemp = _appCutoff(input, r, co, root5)
                    cutoffd[co][r]["M5"] = cutofftemp
                    cutoffList.append(cutofftemp)

    #Applying Filter
    if mirnaref[0] == "no":
        ref_dic = _GetRef(mainref_path)
    else:
        ref_dic = _GetRef(mainref_path)
        mirna_path = mirnaref[1]
        _Round_mir(bowtie_path, mainref_path, mirna_path, root0)
        mature_dic = _mature()

    ##Creating MapResults
    print "\nCreating MapResults"
    lib_dic = dict()
    for r in outListr0:
        lib_dic, cab = _ParseRn(r, lib_dic)
        _cleanDIRS("r0")
    for r in outListr1:
        lib_dic, cab = _ParseRn(r, lib_dic)
        _cleanDIRS("r1")
    if add3_ == "yes":
        for r in outListM3:
            lib_dic, cab = _ParseRn(r, lib_dic)
            _cleanDIRS("M3")
    if add5_ == "yes":
        for r in outListM5:
            lib_dic, cab = _ParseRn(r, lib_dic)
            _cleanDIRS("M5")
    out = "Results/MapResults/All_MapResults.txt"
    _buildout(ref_dic, lib_dic, mature_dic, cab, out)
    for c in cutoffd:
        print "\nCreating MapResults with %s reads Cutoff" % (c)
        lib_dic2 = dict()
        out = "Results/MapResults/cutoff_%s_MapResults.txt" % (c)
        for r in cutoffd[c]:
            for i in cutoffd[c][r]:
                #print "\tGetting frequences for %s" % (cutoffd[c][r][i])
                lib_dic2, cab = _ParseRn(cutoffd[c][r][i], lib_dic2)
        _buildout(ref_dic, lib_dic2, mature_dic, cab, out)


if __name__=="__main__":

    print "\n Start Running isomiRID Version 0.53\n"

    if "-f" in sys.argv or "--force" in sys.argv:
        _cleanALLDIRS()
        _Complete_analysis()
    else:
        _Complete_analysis()