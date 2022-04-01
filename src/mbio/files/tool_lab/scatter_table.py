# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'  
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from biocluster.config import Config
import os
import codecs
import chardet
import random
import re
import subprocess


class ScatterTableFile(File):
    """
    定义小工具二维表格
    """
    def __init__(self):
        super(ScatterTableFile, self).__init__()
        self.unicode = False
        self.col_number = 0
        self.col_sample = []
        self.row_sample = []

    def get_info(self):
        if 'path' in self.prop.keys():
            # file_name = os.path.basename(self.prop['path'])
            # if file_name.split(".")[-1] != 'txt':
            #     raise FileError("文件类型不对,应该为TXT格式,其后缀应该为.txt")  # changed by wzy,20170926  注释掉
            f = open(self.prop['path'], 'r')
            code_dic = chardet.detect(f.read())
            if code_dic['encoding'] != 'ascii' and code_dic['encoding'] != 'UTF-16LE':
                raise FileError('文件编码格式不符合要求')
            self.new_path = self.get_newtable(code_dic['encoding'])
            self.set_property("new_table", self.new_path)
            self.check_info(self.new_path)  # 判断是否符合数据表格的要求
            self.set_property('sample_num', self.col_number)
            self.set_property('col_sample', self.col_sample)
            self.set_property('row_sample', self.row_sample)
        else:
            raise FileError("文件路径不正确，请设置正确的文件路径!")

    def get_newtable(self, encoding):
        dos2unix_path = Config().SOFTWARE_DIR + '/bioinfo/hd2u-1.0.0/bin/dos2unix'
        dir_path = Config().WORK_DIR + '/tmp/convert_table'
        if os.path.exists(dir_path):  # 刚开始还不存在covert_table这个文件夹
            pass
        else:
            os.mkdir(dir_path)

        # 定义新文件的路径
        new_file_path = os.path.join(dir_path, "{}_{}.txt".format(os.path.basename(self.prop['path']),
                                                                  random.randint(1000, 10000)))
        if encoding == 'ascii':  # 如果是ASCII码的时候，根据行数判断是否可能是mac的文件，采取不同的转换方式
            with open(self.prop['path'], 'r') as old:
                line = old.readlines()
            if len(line) == 1:
                subprocess.check_output(dos2unix_path + ' -C ' + self.prop['path'], shell=True)
            else:
                subprocess.check_output(dos2unix_path + ' -U ' + self.prop['path'], shell=True)
            with open(self.prop['path'], 'r') as old:  # 转换完以后判断分隔符是否为\t(tab分隔)，如果不是进行转换
                first_line = old.readline().strip("\n").split("\t")
            if len(first_line) > 1:
                return self.prop['path']
            else:
                with open(self.prop['path'], 'r') as old, open(new_file_path, 'w') as new:
                    for line in old:
                        line = line.strip('\n').split(' ')
                        new_line = []
                        for i in line:
                            if i != '':
                                new_line.append(i)
                            else:
                                continue
                        new.write(('\t').join(new_line) + '\n')
                return new_file_path
        else:  # 编码格式为utf-16-le时，直接读取另写入新的文件路径中
            fp2 = codecs.open(self.prop['path'], 'r', 'utf-16-le')
            lineList = fp2.readlines()
            with open(new_file_path, 'w') as new:
                n = 0
                for line in lineList:
                    n += 1
                    line = line.strip('\r\n').split('\t')
                    line_c = []
                    if n == 1:
                        line = line[1:]
                    for content in line:
                        line_c.append(str(content))
                    if n == 1:
                        new.write("#sample\t" + ('\t').join(line_c) + '\n')
                    else:
                        new.write(('\t').join(line_c) + '\n')
            return new_file_path

    def check_info(self, file_path):  # 判断二维表是否符合要求
        with open(file_path, 'r') as f:
            first_line = f.readline().strip('\r\n').split('\t')
            f1 = set(first_line)
            if len(f1) != len(first_line):
                raise FileError('列名不能重复_{}'.format(first_line))
            if first_line[0] == "":
                raise FileError('矩阵第一行第一列不能为空')  #add by liulinmeng , 20180416
            for m in first_line[1:]:
                if re.match('^[a-zA-Z0-9_]+$', m):
                    continue
                else:
                    raise FileError('列名中只能含数字/字母/下划线_{}'.format(m))
            self.col_number = len(first_line)
            self.col_sample = first_line[1:]
            for i in first_line:
                if i.isdigit():
                    raise FileError('列名中不能存在数字_{}'.format(i))
                else:
                    continue
            row_name = []
            for line in f:
                content = line.strip('\r\n').split('\t')
                if content[0] in row_name:
                    raise FileError('行名不能重复_{}'.format(content[0]))
                else:
                    row_name.append(content[0])
                    if re.match('^[a-zA-Z0-9_]+$', content[0]):
                        pass
                    else:
                        raise FileError('行名中只能含数字/字母/下划线_{}'.format(content[0]))
                if len(content) != self.col_number:
                    raise FileError('该表格行列信息不全——{}'.format(content))
                if content[0].isdigit():
                    raise FileError('行名中不能存在数字')
                # for i in content[1:]:
                #     if float(i) or i == '0':
                #         continue
                #     else:
                #         raise FileError('二维数据表格内容必须为数字_{}'.format(i))
                self.row_sample = row_name

    @staticmethod
    def get_table_of_main_table(origin_table, target_path, group):
        r_sample = []
        final_lines =[]
        with open(group, 'r') as g, open(origin_table)as fr, open(target_path, "w")as fw:
            for line in g:
                sample = line.strip().split("\t")[0]
                r_sample.append(sample)   # 获取分组文件中存在的样本
            lines = fr.readlines()
            table_list = [i.rstrip().split('\t') for i in lines]
            T_lines = map(lambda *a: '\t'.join(a), *table_list)
            # T_lines = map(lambda *a: '\t'.join(a) + '\n', *lines)
            final_lines.append(T_lines[0])
            for line in T_lines[1:]:
                line_split = line.strip().split("\t")
                if line_split[0] not in r_sample:
                    pass
                else:
                    final_lines.append(line)  # 转置挑选需要的样本
            table_list = [i.rstrip().split('\t') for i in final_lines]
            T_lines = map(lambda *a: '\t'.join(a) + '\n', *table_list)  # 再次转置回来
            for line in T_lines:
                fw.write(line)
            # 此方法不能投递运行，故换写法
            # import pandas as pd
            # data = pd.read_table(origin_table, header=0, sep="\t")
            # for i in data.columns[1:]:
            #     if i not in r_sample:
            #         exp.append(i)  # 获取要去除的样本
            # if len(exp) >= 1:
            #     data = data.drop(exp, axis=1)
            # else:
            #     pass
            # data.to_csv(target_path, sep="\t", index=False)  # 数据框写入文件

    def check(self):
        if super(ScatterTableFile, self).check():
            self.get_info()
            return True


if __name__ == "__main__":
    a = ScatterTableFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/wangzhaoyue/toolapps/single_table_input/matrix_column.txt")
    group_file = '/mnt/ilustre/users/sanger-dev/workspace/20170913/Single_test_bar_log_group1/SimpleBar/group1/group_1'
    a.check()
    out = '/mnt/ilustre/users/sanger-dev/workspace/20170913/Single_test_bar_log_group1/SimpleBar/group1/input_1'
    a.get_table_of_main_table('/mnt/ilustre/users/sanger-dev/sg-users/wangzhaoyue/toolapps/single_table_input/matrix_column.txt', out, group_file)
