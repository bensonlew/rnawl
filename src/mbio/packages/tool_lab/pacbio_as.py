# -*- coding: utf-8 -*-
import os
import argparse
from Bio.Blast import NCBIXML
import pandas as pd


def step1(AS_blast_out, outpath):
    output_file1 = os.path.join(outpath, '1.txt')
    with open(AS_blast_out, 'r') as a, open(output_file1, 'w') as o:
        Blast_records = NCBIXML.parse(a)
        header = ['QueryName', 'QueryLength', 'SubjectName', 'SubjectLength', 'QueryIdentityPercent', 'HspLength',
                  'QhStart',
                  'QhEnd', 'ShStart', 'ShEnd', 'HspNum']
        o.write("\t".join(header) + "\n")
        for blast_record in Blast_records:
            queryName = blast_record.query
            queryLength = int(blast_record.query_length)
            subjectName = ""
            subjectLength = ""
            queryIdentityPercent = 0
            for alignment, description in zip(blast_record.alignments, blast_record.descriptions):
                hspNum = int(len(alignment.hsps))
                subjectName = alignment.title
                subjectLength = int(alignment.length)
                this_matches = 0
                for p in alignment.hsps:
                    if p.expect <= 1000000:
                        this_query = blast_record.query_length
                        this_matches += p.identities
                        if this_matches > 0:
                            this_percent = (this_matches * 100.) / this_query
                            queryIdentityPercent = this_percent
                if queryName != subjectName and hspNum == 2:
                    if 1000 <= queryLength <= subjectLength:
                        for hsp in alignment.hsps:
                            hspLength = hsp.align_length
                            qhStart = hsp.query_start
                            qhEnd = hsp.query_end
                            shStart = hsp.sbjct_start
                            shEnd = hsp.sbjct_end
                            outList = [queryName, str(queryLength), subjectName, str(subjectLength),
                                       str(round(queryIdentityPercent, 2)), str(hspLength), str(qhStart), str(qhEnd),
                                       str(shStart), str(shEnd), str(hspNum)]
                            o.write("\t".join(outList) + "\n")
    return output_file1
def step2_1(output_file1, outpath):
    WriteOutFile = True
    output_file2_1 = os.path.join(outpath, '2_1.txt')
    with open(output_file1, 'r') as a, open(output_file2_1, 'w') as o:
        HeaderLine = 'QName\tQLength\tSName\tSLength\tQueryIdentityPercent\tHspLength\tQhspStart\tQhspEnd\tShspStart\tShspEnd\tDiffSLen'
        o.write(HeaderLine + '\n')
        LineNumber = 0

        for Line in a:
            if LineNumber > 0:
                Line = Line.strip('\n')
                ElementList = Line.split('\t')
                queryName = ElementList[0]
                queryLength = ElementList[1]
                subjectName = ElementList[2]
                subjectLength = ElementList[3]
                queryIDPercent = ElementList[4]
                hspLength = ElementList[5]
                Qhstart = ElementList[6]
                QhEnd = ElementList[7]
                ShStart = ElementList[8]
                ShEnd = ElementList[9]
                Delta_S_Len = int(ElementList[9]) - int(ElementList[8])
                outList = [str(queryName), str(queryLength), str(subjectName), str(subjectLength), str(queryIDPercent),
                           str(hspLength), str(Qhstart), str(QhEnd), str(ShStart), str(ShEnd), str(Delta_S_Len)]
                if WriteOutFile:
                    o.write("\t".join(outList) + "\n")
            LineNumber += 1
    return output_file2_1

def step2_2(output_file2_1, outpath):
    output_file2_2 = os.path.join(outpath, '2_2.txt')
    with open(output_file2_1, 'r') as a, open(output_file2_2, 'w') as o:
        HeaderLine = 'QName\tQLength\tSName\tSLength\tQueryIdentityPercent\tHspLength1\tQhspStart1\tQhspEnd1\tShspStart1\tShspEnd1\tDiffSLen1\tHspLength2\tQhspStart2\tQhspEnd2\tShspStart2\tShspEnd2\tDiffSLen2'
        o.write(HeaderLine + '\n')

        LineNumber = 0
        dataDict = dict()
        for Line in a:
            if LineNumber > 0:
                Line = Line.strip('\n')
                ElementList = Line.split('\t')

                keyName = ElementList[0] + "\t" + ElementList[1] + "\t" + ElementList[2] + "\t" + ElementList[3] + "\t" + \
                          ElementList[4]

                if dataDict.__contains__(keyName):
                    dataDict[keyName] = dataDict[keyName] + "\t" + ElementList[5] + "\t" + ElementList[6] + "\t" + \
                                        ElementList[
                                            7] + "\t" + ElementList[8] + "\t" + ElementList[9] + "\t" + ElementList[10]

                else:
                    dataDict[keyName] = ElementList[5] + "\t" + ElementList[6] + "\t" + ElementList[7] + "\t" + ElementList[
                        8] + "\t" + ElementList[9] + "\t" + ElementList[10]
            # print dataDict[keyName]+"\n"
            LineNumber += 1

        for KeyName in dataDict:
            values = dataDict[KeyName]
            outList = [str(KeyName), str(values)]
            o.write("\t".join(outList) + "\n")
    return output_file2_2

def step3(output_file2_2, outpath):
    output_file3 = os.path.join(outpath, '3.txt')
    with open(output_file2_2, 'r') as a, open(output_file3, 'w') as o:
        HeaderLine = 'QueryName\tQueryLength\tSubjectName\tSubjectLength\tQueryIdentity\tHspLength1\tQhspStart1\tQhspEnd1\tShspStart1\tShspEnd1\tDiffSLen1\tHspLength2\tQhspStart2\tQhspEnd2\tShspStart2\tShspEnd2\tDiffSLen2'
        o.write(HeaderLine + '\n')
        LineNumber = 0

        for Line in a:
            if LineNumber > 0:
                Line = Line.strip('\n')
                ElementList = Line.split('\t')
                DiffLen1 = int(ElementList[10])
                DiffLen2 = int(ElementList[16])
                if (DiffLen1 < 0 and DiffLen2 < 0) or (DiffLen1 > 0 and DiffLen2 > 0):
                    o.write(Line + '\n')
            LineNumber += 1
    return output_file3

def step4(output_file3, outpath):
    output_file4 = os.path.join(outpath, '4.txt')
    with open(output_file3, 'r') as a, open(output_file4, 'w') as o:
        HeaderLine = 'QueryName\tSubjectName\tQIdentity\tQhspStart1\tQhspEnd1\tQhspStart2\tQhspEnd2\tShspStart1\tShspEnd1\tShspStart2\tShspEnd2'
        o.write(HeaderLine + '\n')
        LineNumber = 0

        for Line in a:
            if LineNumber > 0:
                Line = Line.strip('\n')
                ElementList = Line.split('\t')
                QName = ElementList[0]
                QLength = ElementList[1]
                SName = ElementList[2]
                Slength = ElementList[3]
                QIdentity = float(ElementList[4])
                Hsp_Length_I = int(ElementList[5])
                Qhsp_Start_I = ElementList[6]
                Qhsp_End_I = ElementList[7]
                Shsp_Start_I = ElementList[8]
                Shsp_End_I = ElementList[9]
                Diff_SLen_I = ElementList[10]
                Hsp_Length_II = int(ElementList[11])
                Qhsp_Start_II = ElementList[12]
                Qhsp_End_II = ElementList[13]
                Shsp_Start_II = ElementList[14]
                Shsp_End_II = ElementList[15]
                Diff_SLen_II = ElementList[16]
                if 95 <= QIdentity <= 100.5:
                    if Hsp_Length_I >= 100 and Hsp_Length_II >= 100:
                        outList = [str(QName), str(SName), str(QIdentity), str(Qhsp_Start_I), str(Qhsp_End_I),
                                   str(Qhsp_Start_II), str(Qhsp_End_II), str(Shsp_Start_I), str(Shsp_End_I),
                                   str(Shsp_Start_II),
                                   str(Shsp_End_II)]
                        o.write("\t".join(outList) + "\n")
            LineNumber += 1
    return output_file4

def step5(output_file4, outpath):
    output_file5 = os.path.join(outpath, '5.txt')
    with open(output_file4, 'r') as a, open(output_file5, 'w') as o:
        HeaderLine = 'QueryName\tSubjectName\tQIdentity\tQhspStart1\tQhspEnd1\tQhspStart2\tQhspEnd2\tShspStart1\tShspEnd1\tShspStart2\tShspEnd2'
        o.write(HeaderLine + '\n')
        LineNumber = 0

        for Line in a:
            if LineNumber > 0:
                Line = Line.strip('\n')
                ElementList = Line.split('\t')
                QName = ElementList[0]
                SName = ElementList[1]
                Qhsp_Start_I = ElementList[3]
                Qhsp_End_I = ElementList[4]
                Qhsp_Start_II = ElementList[5]
                Qhsp_End_II = ElementList[6]
                Shsp_Start_I = ElementList[7]
                Shsp_End_I = ElementList[8]
                Shsp_Start_II = ElementList[9]
                Shsp_End_II = ElementList[10]

                if not (Qhsp_Start_II > Qhsp_Start_I and Qhsp_End_II < Qhsp_End_I) and not (
                        Qhsp_Start_II < Qhsp_Start_I and Qhsp_End_II > Qhsp_End_I):
                    o.write(Line + '\n')
            LineNumber += 1
    return output_file5

def step6_1(output_file5, outpath):
    output_file6_1 = os.path.join(outpath, '6_1.txt')
    with open(output_file5, 'r') as a, open(output_file6_1, 'w') as o:
        HeaderLine = 'QueryName\tSubjectName\tQueryIdentityPercent\tQhspStart1\tQhspEnd1\tQhspStart2\tQhspEnd2\tShspStart1\tShspEnd1\tShspStart2\tShspEnd2'
        o.write(HeaderLine + '\n')
        LineNumber = 0

        for Line in a:
            if LineNumber > 0:
                Line = Line.strip('\n')
                ElementList = Line.split('\t')
                QName = ElementList[0]
                SName = ElementList[1]
                Qhsp_Start_I = ElementList[3]
                Qhsp_End_I = ElementList[4]
                Qhsp_Start_II = ElementList[5]
                Qhsp_End_II = ElementList[6]
                Shsp_Start_I = ElementList[7]
                Shsp_End_I = ElementList[8]
                Shsp_Start_II = ElementList[9]
                Shsp_End_II = ElementList[10]
                GapQueryI = abs(int(Qhsp_End_I) - int(Qhsp_Start_II))
                GapQueryII = abs(int(Qhsp_Start_I) - int(Qhsp_End_II))
                GapSbjctI = abs(int(Shsp_End_I) - int(Shsp_Start_II))
                GapSbjctII = abs(int(Shsp_End_II) - int(Shsp_Start_I))
                if (GapQueryI <= 5 and GapSbjctI >= 100) or (GapQueryII <= 5 and GapSbjctII >= 100):
                    o.write(Line + '\n')
            LineNumber += 1
    return output_file6_1

def step6_2(output_file6_1, output):
    output_path = os.path.join(output, 'output/AS_result.txt')
    out = pd.read_table(output_file6_1, header=0, sep='\t')
    out.sort_values(by=['QueryIdentityPercent'], ascending=False, inplace=True)
    out.drop_duplicates(subset=['SubjectName'], keep='first', inplace=True)
    out.drop_duplicates(subset=['QueryName'], keep='first', inplace=True)
    SubjectName = out['SubjectName'].tolist()
    out['SubjectName'] = [i.rstrip('No definition line') for i in SubjectName]
    out.to_csv(output_path, header=True, index=False, sep='\t')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="请输入blast_result", type=str, required=True)
    parser.add_argument("-o", help="请输入输出路径", type=str, required=True)
    args = parser.parse_args()
    output_file1 = step1(args.i, args.o)
    output_file2_1 = step2_1(output_file1, args.o)
    output_file2_2 = step2_2(output_file2_1, args.o)
    output_file3 = step3(output_file2_2, args.o)
    output_file4 = step4(output_file3, args.o)
    output_file5 = step5(output_file4, args.o)
    output_file6_1 = step6_1(output_file5, args.o)
    step6_2(output_file6_1, args.o)