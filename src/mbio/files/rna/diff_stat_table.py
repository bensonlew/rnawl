# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

import re
import os
from collections import defaultdict
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class DiffStatTableFile(File):
    """
    """

    def __init__(self):
        super(DiffStatTableFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(DiffStatTableFile, self).get_info()
        diff_gene, regulate_dict = self.get_table_info()
        self.set_property('diff_genes', diff_gene)
        self.set_property('regulate_dict', regulate_dict)

    def get_table_info(self):
        """
        获取并返回table的信息
        :return:
        """
        with open(self.prop['path']) as f:
            header = f.readline().strip()
            if not header.lower().endswith('pvalue\tpadjust\tsignificant\tregulate\tncbi'):
                raise Exception('错误的表头:{}，表头应与基因差异表达分析结果格式一致'.format(header))
            diff_gene = []
            regulate_dict = defaultdict(list)
            for line in f:
                line = line.strip('\n').split('\t')
                if line[-3] == 'yes':
                    diff_gene.append(line[0])
                if line[-2] == 'up':
                    regulate_dict['up'].append(line[0])
                if line[-2] == 'down':
                    regulate_dict['down'].append(line[0])
        return diff_gene, regulate_dict

    def get_stat_file(self, output, kegg_file):
        f = open(kegg_file, "r")
        ko_list = []
        f.readline()
        for line in f:
            line = line.strip().split("\t")
            ko_list.append(line[0])
        f.close()
        new_file_name = output + "/" + os.path.basename(self.prop["path"]).split("_edgr_stat.xls")[0] + ".DE.list"
        with open(self.prop['path'], "r") as f, open(new_file_name, "w") as w:
            header = f.readline().strip()
            for line in f:
                line = line.strip('\n').split('\t')
                if line[0] in ko_list:
                    w.write(line[0] + "\t" + line[-2] + "\n")

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(DiffStatTableFile, self).check():
            # 父类check方法检查文件路径是否设置，文件是否存在，文件是否为空
            self.get_info()

if __name__ == "__main__":
    cls = DiffStatTableFile()
    cls.set_path("/mnt/ilustre/users/sanger-dev/workspace/20170613/Refrna_ore_test_for_api/Express/output/diff/trans_diff/diff_stat_dir/A_vs_B_edgr_stat.xls")
    cls.check()
    cls.get_stat_file()