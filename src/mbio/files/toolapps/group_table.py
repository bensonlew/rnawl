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
from mbio.files.toolapps.table import TableFile


class GroupTableFile(TableFile):
    """
    定义小工具二维分组文件不存在数字
    """
    def __init__(self):
        super(GroupTableFile, self).__init__()
        self.new_file = ''

    def get_info(self):
        if 'path' in self.prop.keys():
            file_name = os.path.basename(self.prop['path'])
            if file_name.split(".")[-1] != 'txt':
                raise FileError("文件类型不对,应该为TXT格式,其后缀应该为.txt")
            f = open(self.prop['path'], 'r')
            code_dic = chardet.detect(f.read())
            if code_dic['encoding'] != 'ascii' and code_dic['encoding'] != 'UTF-16LE':
                raise FileError('文件编码格式不符合要求')
            self.set_property("new_table", self.get_newtable(code_dic['encoding']))
            self.new_file = self.get_newtable(code_dic['encoding'])
            self.check_info(self.new_file)  # 判断是否符合数据表格的要求
            info = self.get_file_info()
            self.set_property("group_scheme", info[1])  # 分组方案中涉及的分组名称
            self.set_property('sample_name', info[0])  # 分组方案中涉及的样本
        else:
            raise FileError("文件路径不正确，请设置正确的文件路径!")

    def get_file_info(self):
        """
		获取group_table文件的信息
		"""
        row = 0
        with open(self.new_file, 'r') as f:
            sample = list()
            line = f.readline().rstrip()  # 将rstrip("\r\n") 全部替换为rstrip()
            line = re.split("\t", line)
            if line[1] == "##empty_group##":
                is_empty = True
            else:
                is_empty = False
            header = list()
            len_ = len(line)  # 标明group文件有几列
            for i in range(1, len_):  # 记录列名,也就是记录分组方案的名称
                header.append(line[i])
            for line in f:
                line = line.rstrip()
                line = re.split("\t", line)
                row += 1
                if line[0] not in sample:
                    sample.append(line[0])  # 记录所有分组方案中的样本名称
            return (sample, header, is_empty)

    def sub_group(self, target_path, header):
        """
		:param target_path:  生成的子group表的位置
		:param header: 需要提取的子分组方案名，列表
		"""
        if not isinstance(header, list):
            raise Exception("第二个参数的格式错误， 应该是一个python的列表")
        my_index = list()
        for h in header:
            if h not in self.prop['group_scheme']:
                raise Exception("{}不存在该表的分组方案中".format(h))
        len_ = len(self.prop['group_scheme'])
        for i in range(0, len_):
            if self.prop['group_scheme'][i] in header:
                my_index.append(i + 1)
        with open(self.new_file, 'r') as f, open(target_path, 'w') as w:
            line = f.readline().rstrip()
            line = re.split("\t", line)
            new_header = list()
            for i in my_index:
                new_header.append(line[i])
            w.write("#sample\t{}\n".format("\t".join(new_header)))
            for line in f:
                sub_line = list()
                line = line.rstrip()
                line = re.split("\t", line)
                sub_line.append(line[0])
                for i in my_index:
                    if line[i] == '' or line[i] == ' ':  # add by wzy 20170909
                        pass
                    else:
                        sub_line.append(line[i])
                if len(sub_line) > 1:  # add by wzy 20170909
                    new_line = "\t".join(sub_line)
                    w.write(new_line + "\n")
                else:
                    pass

    def get_group_detail(self):
        """
        根据分组文件得到具体的分组方案
        """
        group_samples = {}  # 分组对应中新样本对应的旧样本
        with open(self.prop['path'], "r") as f:
            line = f.readline().rstrip()
            line = re.split("\t", line)
            if line[1] == "##empty_group##":
                is_empty = True
            else:
                is_empty = False
            for i in range(1, len(line)):
                group_samples[line[i]] = {}
            for item in f:
                item = item.rstrip().split("\t")
                for i in range(1, len(line)):
                    try:
                        if item[i] and item[i] not in group_samples[line[i]]:
                            group_samples[line[i]][item[i]] = []
                            group_samples[line[i]][item[i]].append(item[0])
                        elif item[i]:
                            group_samples[line[i]][item[i]].append(item[0])
                        else:
                            print "{}样本不在分组方案{}内".format(item[0], line[i])
                    except:
                        print "{}样本不在分组方案{}内".format(item[0], line[i])
        return group_samples

    @staticmethod
    def check_info(file_path):
        with open(file_path, 'r') as f:
            first_line = f.readline().strip('\n').split('\t')
            f1 = set(first_line)
            if len(f1) != len(first_line):
                raise FileError('列名不能重复_{}'.format(first_line))
            col_number = len(first_line)
            for i in first_line:
                if i.isdigit():
                    raise FileError('列名中不能存在数字——{}'.format(i))
                else:
                    continue
            row_name = []
            for line in f:
                content = line.strip('\n').split('\t')
                if content[0] in row_name:
                    raise FileError('行名不能重复_{}'.format(content[0]))
                else:
                    row_name.append(content[0])
                if len(content) != col_number:
                    raise FileError('该表格行列信息不全——{}'.format(content))
                for i in content:
                    if i.isdigit():
                        raise FileError('行名中不能存在数字——{}'.format(content))
                    else:
                        pass
                    if i == '' or i == ' ':
                        raise FileError('分组中样本不能为空——{}'.format(content))

    def group_num(self, name):
        """
        传入分组方案name,判断该分组方案下分组类别的数量
        :param name: 某一分组方案的名字,string
        :return group:该分组方案下分组类别的数量
        """
        with open(self.prop['path'], 'r') as f:
            head = f.readline().split()
            site = 0
            for i in head:
                if i == name:
                    site = head.index(i)
            gnames = []
            while True:
                line = f.readline()
                if not line:
                    break
                gnames.append(line.split()[site])
            gnum = {}
            for n in gnames:
                gnum[n] = gnames.count(n)
            group = len(gnum)
            return group

    def check(self):
        if super(GroupTableFile, self).check():
            self.get_info()
            return True

if __name__ == '__main__':
    a = GroupTableFile()
    # a.set_path('/mnt/ilustre/users/sanger-dev/sg-users/wangzhaoyue/toolapps/single_table_input/group2.txt')
    a.set_path('/mnt/ilustre/users/sanger-dev/sg-users/zengjing/toolapps/hc_heatmap/more_group.txt')
    a.get_info()
    a.sub_group('/mnt/ilustre/users/sanger-dev/sg-users/wangzhaoyue/toolapps/group1.txt', ['BBB'])
    group_detail = a.get_group_detail()
    print group_detail
