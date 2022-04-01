# -*- coding: utf-8 -*-
# __author__ = 'gdq'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import OrderedDict

class GroupTableFile(File):
    """
    group_info: path of group info, file with at least two columns. if no replicate
               exist, just use sample name as group name. Header line starts with '#'.
               --------------------
               #sample group_name  group_name
               s1   group1
               s2   group1  group3
               s3   group2
               s4   group2  group3
               s5   s5
               s6   s6
               --------------------
    """
    def __init__(self):
        super(GroupTableFile, self).__init__()

    def get_info(self):
        super(GroupTableFile, self).get_info()
        header, group_dict, samples = self.parse_file()
        self.set_property("sample_number", len(samples))
        self.set_property("sample", samples)
        self.set_property("group_scheme", header)
        self.set_property("group_dict", group_dict)

    def parse_file(self):
        # group_info -> dict, group_name as key, list of sample names as values. {group:[s1,s2,]}
        sample_list = set()
        group_dict = OrderedDict()
        header_list = list()
        with open(self.prop['path']) as f:
            header_line = f.readline()
            if not header_line.startswith('#'):
                raise FileError("该group文件不含表头，group表第一列应该以#号开头", code = "42001201")
            header_list = header_line.lstrip('#').rstrip().split()[1:]
            for line in f:
                if not line.strip():
                    continue
                tmp_list = line.strip().split()
                sample_list.add(tmp_list[0])
                for g in tmp_list[1:]:
                    group_dict.setdefault(g, list())
                    group_dict[g].append(tmp_list[0])
            for g in group_dict.keys():
                group_dict[g] = sorted(list(set(group_dict[g])))
        if not group_dict:
            raise FileError("分组信息内容为空", code = "42001202")
        return header_list, group_dict, sorted(sample_list)

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

    def check(self):
        super(GroupTableFile, self).check()
        self.get_info()
