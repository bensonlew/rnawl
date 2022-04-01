#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__ = 'Wanjin.Hu'

import argparse
import os
import re
import datetime
from pathlib import Path

parser = argparse.ArgumentParser(description="get sample raw data path of bacterialGenome and microbial longReads")
parser.add_argument("-i", "--input", dest="inList", required=True,
        type=str, help="input projectPathTable, such as:/mnt/ilustre/centos7users/wanjin.hu/script_test/pacbioDataTable.xls")
parser.add_argument("-o", "--output", dest="outList", required=True,
        type=str, help="output samplePathTable, such as:/mnt/ilustre/centos7users/wanjin.hu/script_test/result.xls")
parser.add_argument("-t", "--tempdir",dest="tempDir", required=False, default="/mnt/ilustre/centos7users/wanjin.hu/script_test/testData",
        type=str, help="Temporarily stored folder for trim sequences, such as:/mnt/ilustre/centos7users/wanjin.hu/script_test/testData")
args = parser.parse_args()

# trim_fqSeq="/mnt/ilustre/centos7users/wanjin.hu/script/trim_fqSeq_pacbio.pl"
# fastq_rename = "/mnt/ilustre/centos7users/wanjin.hu/script/fastqRename.py"

def new_split_path(split_path):
    if os.path.exists(split_path):
        return split_path
    if "ilustre" in split_path:
        split_path1 = split_path.replace("ilustre", "clustre")
        if os.path.exists(split_path1):
            return split_path1
    if "sglustre" in split_path:
        split_path1 = split_path.replace("sglustre", "ilustre")
        if os.path.exists(split_path1):
            return split_path1
    return split_path

start = datetime.datetime.now()
dp = open(args.outList,'w')
with open(args.inList, 'r') as list:
    #next(list)
    r = list.readline().strip('\n')
    dp.write(r+'\trawdata path\ttrimdata path\tmerge\ttrim\trename\n')
    for line in list:
        dataPath = ''
        trimDataPath = ''
        merge, trim, rename = 'false', 'false', 'false'
        wholeLine = line.strip("\n")
        projectType = wholeLine.split('\t')[0]
        projectNum = wholeLine.split('\t')[2]
        mjNum = wholeLine.split('\t')[4]
        sampleName = wholeLine.split('\t')[5]
        samplePrimer = wholeLine.split('\t')[6]
        #sequencePlatform = wholeLine.split('\t')[10]
        #sequenceCompany = wholeLine.split('\t')[14]
        samplePath = wholeLine.split('\t')[-2]
        # if not os.path.exists(samplePath):
        #     samplePath = samplePath.replace("ilustre", "clustre")
        samplePath = new_split_path(samplePath)
        samplePath1 = Path(samplePath)
        if projectType in ['细菌基因组完成图', '真菌基因组精细图', "微生物基因组纯测序"]:
            trimDataPath = '-'
            if samplePath == '':
                print('样本{}的路径为空，请查看表格{}中的项目路径是否完整\n'.format(mjNum, args.inList))
                dataPath = '-'
            else:
                if samplePath1.is_dir():
                    for subFile in os.listdir(samplePath):
                        fileName = samplePath.split("/")[-1]
                        if subFile == '1.rawdata':
                            fileList=[]
                            samplePath2 = os.path.join(samplePath,subFile)
                            for nsubFile in os.listdir(samplePath2):
                                if nsubFile.endswith ('.bam'):
                                    tempDataPath = os.path.join(samplePath2,nsubFile)
                                    fileList.append(tempDataPath)
                            if len(fileList) == 1:
                                dataPath = fileList[0]
                            else:
                                print('样本{}存在补测的情况，补测数据路径为{}'.format(mjNum, samplePath))
                                print('\n------样本{}正在合并中------\n'.format(mjNum))
                                # os.system('mkdir -p '+args.tempDir+'/'+fileName)
                                # os.system('pbmerge -o '+args.tempDir+'/'+fileName+'/merge.bam '+samplePath2+'/*.bam')
                                # dataPath = os.path.join(args.tempDir,fileName,'merge.bam')
                                dataPath = ";".join(fileList)
                                merge = fileName+".merge.bam"
                        if subFile.endswith ('.bam'):
                            dataPath = os.path.join(samplePath,subFile)
                        if 'pass' in subFile:
                            dataPath = os.path.join(samplePath,subFile)
                        else:
                            dataPath == '-'
        if projectType == '全长微生物多样性分析':
            # valuedSampleName = sampleName+'_valued.fastq'
            valuedSampleName = projectNum+'-'+mjNum+'-'+sampleName+'-'+samplePrimer+'_valued.fastq'
            if samplePath == '':
                print('样本{}的路径为空，请查看表格{}中的项目路径是否完整\n'.format(mjNum, args.inList))
                dataPath = '-'
                trimDataPath = '-'
            else:
                if samplePath1.is_dir():
                    for subFile in os.listdir(samplePath):
                        # if subFile.startswith(mjNum) and subFile.endswith('ccs.fastq'):
                        if subFile.startswith(mjNum) and (subFile.endswith('.fastq') or subFile.endswith('.fq')):
                            dataPath = os.path.join(samplePath,subFile)
                            if samplePrimer == '27F_1492R' or '20F_1492R':
                                print('样本{}正在进行长度质控，请稍等'.format(subFile))
                                # os.system('perl {} -i {} -o {}/trim.{} -m {} -x {}'.format(trim_fqSeq,dataPath,args.tempDir,subFile,1000,1800))
                                # trimDataName = 'trim.'+subFile
                                # preTrimDataPath = os.path.join(args.tempDir,trimDataName)
                                trim = 'trim.'+subFile+":1000-1800"
                                print('样本{}正在进行序列名修改，请稍等'.format(subFile))
                                # os.system('python3 {} -i {} -o {}/{} -n {}'.format(fastq_rename,preTrimDataPath,args.tempDir,valuedSampleName,sampleName))
                                # trimDataPath = os.path.join(args.tempDir,valuedSampleName)
                                rename = valuedSampleName + ":" + sampleName
                            if samplePrimer == 'ITS1F_ITS4R':
                                print('样本{}正在进行长度质控，请稍等'.format(subFile))
                                # os.system('perl {} -i {} -o {}/trim.{} -m {} -x {}'.format(trim_fqSeq,dataPath,args.tempDir,subFile,300,900))
                                # trimDataName = 'trim.'+subFile
                                # preTrimDataPath = os.path.join(args.tempDir,trimDataName)
                                trim = 'trim.'+subFile+":300-900"
                                print('样本{}正在进行序列名修改，请稍等'.format(subFile))
                                # os.system('python3 {} -i {} -o {}/{} -n {}'.format(fastq_rename,preTrimDataPath,args.tempDir,valuedSampleName,sampleName))
                                # trimDataPath = os.path.join(args.tempDir,valuedSampleName)
                                rename = valuedSampleName + ":" + sampleName
                        # if subFile.startswith(sampleName) and subFile.endswith('ccs.fastq'):
                        # if subFile.startswith(sampleName+'.') and (subFile.endswith('.fastq') or subFile.endswith('.fq')):
                        #     dataPath = os.path.join(samplePath,subFile)
                        #     if samplePrimer == '27F_1492R' or '20F_1492R':
                        #         print('样本{}正在进行长度质控，请稍等'.format(subFile))
                        #         # os.system('perl {} -i {} -o {}/trim.{} -m {} -x {}'.format(trim_fqSeq,dataPath,args.tempDir,subFile,1000,1800))
                        #         # trimDataName = 'trim.'+subFile
                        #         # preTrimDataPath = os.path.join(args.tempDir,trimDataName)
                        #         trim = 'trim.'+subFile+":1000-1800"
                        #         print('样本{}正在进行序列名修改，请稍等'.format(subFile))
                        #         # os.system('python3 {} -i {} -o {}/{} -n {}'.format(fastq_rename,preTrimDataPath,args.tempDir,valuedSampleName,sampleName))
                        #         # trimDataPath = os.path.join(args.tempDir,valuedSampleName)
                        #         rename = valuedSampleName + ":" + sampleName
                        #     if samplePrimer == 'ITS1F_ITS4R':
                        #         print('样本{}正在进行长度质控，请稍等'.format(subFile))
                        #         # os.system('perl {} -i {} -o {}/trim.{} -m {} -x {}'.format(trim_fqSeq,dataPath,args.tempDir,subFile,300,900))
                        #         # trimDataName = 'trim.'+subFile
                        #         # preTrimDataPath = os.path.join(args.tempDir,trimDataName)
                        #         trim = 'trim.'+subFile+":300-900"
                        #         print('样本{}正在进行序列名修改，请稍等'.format(subFile))
                        #         # os.system('python3 {} -i {} -o {}/{} -n {}'.format(fastq_rename,preTrimDataPath,args.tempDir,valuedSampleName,sampleName))
                        #         # trimDataPath = os.path.join(args.tempDir,valuedSampleName)
                        #         rename = valuedSampleName + ":" + sampleName
        if "PacBio全长转录组测序" in projectType:  # 上传tar压缩文件
            if samplePath1.is_dir():
                for subFile in os.listdir(samplePath):
                    if subFile.endswith(".tar"):
                        dataPath = os.path.join(samplePath,subFile)
        dp.write(wholeLine+'\t'+dataPath+'\t'+trimDataPath+'\t'+merge+'\t'+trim+'\t'+rename+'\n')
    dp.close()
    print('\n结果文件查看{}'.format(args.outList))
end = datetime.datetime.now()
runningTime = end -start
print('totally run time is {}:'.format(runningTime))
