# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2021.01.20


import re, Bio, urllib2, regex, os
import subprocess
from biocluster.iofile import Directory
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError
import re


'''
菌种鉴定小工具参数gff文件夹有list.txt检查
'''


class GffDirFile(Directory):
    """
    定义gff_dir文件夹
    """
    def __init__(self):
        super(GffDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(GffDirFile, self).get_info()

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        list_txt = os.path.join(self.prop['path'], "list.txt")
        sample = []
        files = []
        if os.path.exists(list_txt):
            with open(list_txt, 'r') as f:
                lines = f.readlines()
                for i in lines[1:]:
                    lin = i.strip().split("\t")
                    if len(lin) != 2:
                        raise FileError('list.txt必须是两列！')
                    else:
                        files.append(lin[0])
                        sample.append(lin[1])
        filelist = os.listdir(self.prop['path'])
        if len(sample) != len(set(sample)):
            raise FileError('样品名称有一样的！')
        if len(files) != len(filelist) - 1:
            raise FileError('list文件记录个数和文件夹下的个数不一致！')
        else:
            for i in files:
                if i not in filelist:
                    raise FileError('文件夹下{}不在list.txt文件中！'.format(i))
        if not len(filelist):
            raise FileError('gff文件夹为空，请检查确认')
        self.check_file2()

    def check_file2(self):
        """
        检查gff文件是否有必须的信息
        :return:
        """
        dig_pat = re.compile('^\d+$')
        files  = os.listdir(self.prop['path'])
        has_cds = False
        normal_gene_num = 0
        for f in files:
            if f.endswith('.gff'):
                gff = self.prop['path'] + '/' + f
                with open(gff) as fr:
                    line_n = 0
                    for line in fr:
                        line_n += 1
                        if line.startswith('#'):
                            continue
                        new_line = line.strip()
                        sp_line = new_line.split('\t')
                        if len(sp_line) == 0: #允许空行存在，不报错
                            continue
                        if self.is_whitespace(line):
                            raise FileError('文件尾部存在空行，请检查！')
                        type = sp_line[2]
                        if type not in ['CDS','tRNA','rRNA']:
                            continue
                        if type == 'CDS':
                            has_cds = True

                        if len(sp_line) != 9:
                            raise FileError('%s 文件第%s 行少于9列'%(gff, line_n))
                        if sp_line[0] == "" :
                            raise FileError('%s 文件第%s 行 第1列不能为空'%(gff,line_n))

                        if not (dig_pat.match(sp_line[3]) and dig_pat.match(sp_line[4])):
                            raise  FileError('%s 文件第%s行 第4和第5列必须都是数字'%(gff,line_n))

                        #检查第7列
                        if type == 'CDS' or type == 'rRNA':
                            if sp_line[6] not in ['+','-','.']:
                                raise  FileError('%s 文件第%s行 第7列必须都是+或-'%(gff,line_n))

                        ##检查第9列
                        if type == 'CDS':
                            if ';Parent=' not in  sp_line[8]:
                                raise FileError('%s 文件第%s行 第9列必须有Parent信息'%(gff,line_n))
                                #print('%s 文件第%s行 第9列必须有Parent信息'%(gff,line_n))
                            else:
                                normal_gene_num +=1
                        elif type == 'tRNA':
                            if ';Parent=' not in  sp_line[8]:
                                #raise FileError('%s 文件第%s行 第9列必须有Parent'%(gff,line_n))
                                print('%s 文件第%s行 第9列必须有Parent'%(gff,line_n))
                        elif type == 'rRNA':
                            if ';Parent=' not in  sp_line[8]:
                                #raise FileError('%s 文件第%s行 第9列必须有Parent'%(gff,line_n))
                                print('%s 文件第%s行 第9列必须有Parent'%(gff,line_n))

                if has_cds == False:
                    raise FileError('%s没有CDS 信息'%(gff))
                if normal_gene_num ==0:
                    raise FileError('%sCDS 都没有Parent'%(gff))

        print('check_file2 检查通过')

    def is_whitespace(self,line):
        if line != ' ' and line != '\n':
            return False
        else:
            return True


if __name__ == "__main__":
    raw_dir = GffDirFile()
    raw_dir.get_info()
    raw_dir.check()
    raw_dir.check_file1()
    raw_dir.check_file2()