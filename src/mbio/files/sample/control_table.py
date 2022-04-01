# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import Counter


class ControlTableFile(File):
    """
    """

    def __init__(self):
        super(ControlTableFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(ControlTableFile, self).get_info()
        num, vs_list = self.get_control_info()
        if num == 0:
            raise FileError('对照方案至少为1', code = "43900101")
        self.set_property('vs_list', vs_list)
        cmp_list = self.parse_file()
        self.set_property('cmp_list', cmp_list)

    def get_control_info(self):
        """
        :return:对照样本（组）数目：num；包含两两比较的样本（分组）元组的列表:[(对照，实验), (对照，实验)]
        """
        with open(self.prop['path'], 'rb') as r:
            lines = r.readlines()
            if not lines:
                raise FileError('对照组文件为空', code = "43900102")
            if not re.match(r'^#', lines[0]):
                raise FileError('对照文件格式有误，表头应为#开头', code = "43900103")
            num = len(lines) - 1
            vs_list = []
            sam_list = []
            for line in lines[1:]:
                if not line.strip():
                    continue
                control = line.strip('\n').split()[0]
                other = line.strip('\n').split()[1]
                vs_list.append((control, other))
                sam_list.append((other, control))
                if control == other:
                    raise FileError('对照组：%s与实验组：%s名字相同！', variables = (control, other), code = "43900104")
            sam_list += vs_list
            count = Counter(sam_list).values()
            for i in count:
                if i != 1:
                    raise FileError("同一个两两比较分组中出现不同的对照组，请检查", code = "43900105")
            return num, vs_list

    def parse_file(self):
        # comparison info -> list. [(ctrl, test), ...]
        with open(self.prop['path']) as f:
            cmp_list = list()
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                tmp_ctrl, tmp_test = line.strip().split()
                cmp_list.append((tmp_ctrl, tmp_test))
        cmp_list = sorted(list(set(cmp_list)))
        if not cmp_list:
            raise FileError('比较信息的内容为空', code = "43900106")
        return cmp_list

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(ControlTableFile, self).check():
            # 父类check方法检查文件路径是否设置，文件是否存在，文件是否为空
            self.get_info()
            return True
