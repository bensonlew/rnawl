# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class BlastTableRelaxedFile(File):
    """
    定义blast+比对输出类型为6结果文件的table格式 注意：文件实际并非来之blast，而是通过xml转化，参见mibo.packages.align.blast.xml2table
    """

    def __init__(self):
        super(BlastTableRelaxedFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(BlastTableRelaxedFile, self).get_info()
        table_info = self.get_table_info()
        self.set_property('count', table_info)

    def get_table_info(self):
        """
        获取blast结果table的信息

        :return
        """
        with open(self.path) as f:
            header = f.readline().strip().split()
            header_len = len(header)
            count = 1
            for line in f:
                count += 1
                line_sp = line.split('\t')
                if len(line_sp) != header_len:
                    raise FileError('文件中存在格式不整齐的行:{}'.format(line))
        return count


    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        if super(BlastTableRelaxedFile, self).check():
            # 父类check方法检查文件路径是否设置，文件是否存在，文件是否为空
            self.get_info()
            return True


if __name__ == '__main__':  # for test
    a = BlastTableRelaxedFile()
    # a.set_path('C:\\Users\\sheng.he.MAJORBIO\\Desktop\\annotation\\annotation\\NR\\transcript.fa_vs_nr.blasttable.xls')
    a.set_path('C:\\Users\\sheng.he.MAJORBIO\\Desktop\\blast_test_1.xls')
    a.check()
    a.get_info()
    print(a.prop['header'])
    print(a.prop['count'])
    print(a.prop['query_num'])
    a.retain_highest_score('C:\\Users\\sheng.he.MAJORBIO\\Desktop\\test_1.xls')
    a.fliter_evalue(0.00001, 'C:\\Users\\sheng.he.MAJORBIO\\Desktop\\test_2.xls')
