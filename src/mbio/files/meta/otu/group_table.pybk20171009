# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import defaultdict


class GroupTableFile(File):
    """
    定义group_table格式文件
    """
    def __init__(self):
        super(GroupTableFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(GroupTableFile, self).get_info()
        info = self.get_file_info()
        self.set_property("sample_number", len(info[0]))
        self.set_property("sample", info[0])
        self.set_property("group_scheme", info[1])
        self.set_property("is_empty", info[2])

    def get_file_info(self):
        """
        获取group_table文件的信息
        """
        row = 0
        self.format_check()
        is_empty = False
        with open(self.prop['path'], 'r') as f:
            sample = list()
            line = f.readline().rstrip()  # 将rstrip("\r\n") 全部替换为rstrip()
            line = re.split("\t", line)
            if line[1] == "##empty_group##":
                is_empty = True
            else:
                is_empty = False
            header = list()
            len_ = len(line)
            for i in range(1, len_):
                header.append(line[i])
            for line in f:
                line = line.rstrip()
                line = re.split("\t", line)
                row += 1
                if line[0] not in sample:
                    sample.append(line[0])
                else:
                    raise FileError("存在相同的样本名")
            return (sample, header, is_empty)

    def format_check(self):
        with open(self.prop['path'], 'r') as f:
            line = f.readline().rstrip()
            if not re.search("^#", line[0]):
                raise FileError("该group（或环境因子）文件不含表头，group（或环境因子）表第一列应该以#号开头")
            line = line.split("\t")
            length = len(line)
            if length < 2:
                raise FileError('分组（或环境因子）文件至少应该有两列')
            for i in line[1:]:
                if re.search("\s", i):
                    raise FileError('分组（或环境因子）文件表头必须以字母或数字开头，不能含有空格，支持字母，数字，下划线')
                if re.search("\.", i):
                    raise FileError("分组（或环境因子）文件表头必须以字母或数字开头，不能包含'.'，支持字母，数字，下划线")
        with open(self.prop['path'], 'r') as f:
            for line in f:
                if "#" in line:
                    continue
                line = line.rstrip()
                line = re.split("\t", line)
                for l in line:
                    if re.search("\s", l):
                        raise FileError('分组名里不可以包含空格')
                    s = ''
                    for i in range(48,58) + [95] + range(97,123) + range(65,91):
                        s += chr(i)
                    check_str = l.rstrip(s)
                    if check_str != '':
                        print('错误分组名：' + l)
                        raise FileError('分组名里不可以有符号、数字、下划线以外的字符' + check_str)
                    if l[0] in s[0:11]:
                        print('分组名建议不要以数字或下划线开头：' + l)
                len_ = len(line)
                if len_ != length:
                    raise FileError("文件的列数不相等")

    def check(self):
        if super(GroupTableFile, self).check():
            self.get_info()
            if not self.prop['is_empty']:
                if self.prop['sample_number'] == 0:
                    raise FileError('应该至少包含一个样本')

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
        with open(self.prop['path'], 'r') as f, open(target_path, 'w') as w:
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
                    sub_line.append(line[i])
                new_line = "\t".join(sub_line)
                w.write(new_line + "\n")

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

    def get_edger_group(self, name, edger_path):
        """
        将分组文件转换成edger分组文件的格式
        """
        self.sub_group('./group_file', name)
        with open('group_file', 'rb') as r, open(edger_path, 'wb') as w:
            lines = r.readlines()
            for line in lines[1:]:
                info = line.strip('\n').split('\t')
                w.write('%s\t%s\n' % (info[1], info[0]))

    def get_group_name(self, name):
        """
        传入分组方案name,获取该分组方案下分组类别的列表
        :param name: 某一分组方案的名字,string
        :return gnames:该分组方案下分组类别的列表
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
            return list(set(gnames))

    def get_group_spname(self, name=None):
        """
        获取某一分组方案下分组类别对应的样本信息详细
        :param name: 某一分组方案的名字,string
        :return group_spname:该分组方案下分组类别对应的样本信息详细的字典，eg：{'A': [1,2,3], 'B': [4,5,6]}
        add by qiuping, last_modify:20161027
        """
        with open(self.prop['path'], 'r') as f:
            head = f.readline().split()
            group_spname = defaultdict(list)
            if not name:
                gname = head[1]
            else:
                gname = name
            _index = head.index(gname)
            for line in f:
                line = line.strip('\n').split()
                sample = line[0]
                group = line[_index]
                group_spname[group].append(sample)
        return group_spname

if __name__ == "__main__":
    g = GroupTableFile()
    g.set_path("example.group")
    g.get_info()
    g.sub_group("example.group.sub", ["g1", "g3", "g4"])
