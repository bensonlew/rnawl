# -*- coding: utf-8 -*-
# __author__ = zouguanqing
# last_modify：2018.03.16


import re, Bio, urllib2, regex, os
import subprocess
from biocluster.iofile import Directory
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError
import re


'''
微生物基因组参数gff文件夹检查
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
        #seqinfo = self.get_raw_info()
        #self.set_property("sample_list", seqinfo[0])

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        filelist = os.listdir(self.prop['path'])
        gff_list_path = os.path.join(self.prop['path'], "list.txt")

        if not len(filelist):
            raise FileError('gff文件夹为空，请检查确认')
        if not os.path.exists(gff_list_path):
            raise FileError('gff文件夹的list.txt为不存在，请检查确认')
        os.system('dos2unix %s'%gff_list_path)  #转unix格式

        self.check_file2()

        print('check 检查通过')


    def check_file1(self):
        """
        检查gff文件的命名
        :return:
        """
        gff_list_path = os.path.join(self.prop['path'], "list.txt")
        pat = re.compile('\s+')
        with open(gff_list_path, "rb") as l:
            lines = l.readlines()
            for line in lines[1:]:
                line2 = pat.split(line.strip())
                if len(line2) == 2:
                    if not os.path.exists(self.prop['path'] + '/'+line2[0]+'.gff'):
                        print (self.prop['path'] + '/'+line2[1]+'.gff')
                        raise FileError('gff文件夹的gff文件名称必须是样本名.gff。样本名.gff 不存在')
                else:
                    raise FileError('list.txt文件格式有误')
        print('check_file1 检查通过')

    def check_file2(self):
        """
        检查gff文件是否有必须的信息
        :return:
        """
        dig_pat = re.compile('^\d+$')
        files  = os.listdir(self.prop['path'])
        has_cds = False
        normal_gene_num = 0
        for  f in files:
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
                            if sp_line[6] not in ['+','-']:
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

    def get_raw_info(self):
        """
        获取dir文件夹的样品信息
        :return: (sample_list)
        """
        sample_list ={}
        raw_list_path = os.path.join(self.prop['path'], "list.txt")
        with open(raw_list_path, "rb") as l:
            lines = l.readlines()
            for line in lines[1:]:
                line2 = line.strip().split("\t")
                sample_list[line2[0]] = line2[0]
        return sample_list



if __name__ == "__main__":
    raw_dir = GffDirFile()
    raw_dir.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/bac/bac_update/test_anno/gff_dir/")
    raw_dir.get_info()
    raw_dir.check()
    raw_dir.check_file1()
    raw_dir.check_file2()
