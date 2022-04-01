# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from mbio.packages.align.blast.xml2table import default_header


class BlastTableFile(File):
    """
    定义blast+比对输出类型为6结果文件的table格式 注意：文件实际并非来之blast，而是通过xml转化，参见mibo.packages.align.blast.xml2table
    """

    def __init__(self):
        super(BlastTableFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(BlastTableFile, self).get_info()
        table_info = self.get_table_info()
        self.set_property('count', table_info[0])
        self.set_property('header', table_info[1])
        self.set_property('query_num', table_info[2])
        self.set_property('query_list', table_info[3])

    def get_table_info(self):
        """
        获取blast结果table的信息

        :return
        """
        with open(self.path) as f:
            header = f.readline().strip().split()
            if header != default_header:
                raise FileError("文件的表头与默认表头不一致：\n默认表头：%s，\n文件错误表头：%s",
                                variables = ('\t'.join(default_header), '\t'.join(header)),
                                code = "45000201")
            count = 0
            query_count = 0
            query_list = []
            flag = None
            for line in f:
                count += 1
                line_sp = line.split('\t')
                if len(line_sp) != 16:
                    raise FileError("文件中存在格式不整齐的行：%s", variables = (line), code = "45000202")
                if flag == line_sp[5]:
                    pass
                else:
                    query_list.append(line_sp[5])
                    query_count += 1
                    flag = line_sp[5]
        return count, default_header, query_count, query_list



    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        if super(BlastTableFile, self).check():
            # 父类check方法检查文件路径是否设置，文件是否存在，文件是否为空
            self.get_info()
            return True

    def retain_highest_score(self, outputfile):
        """保留每个query的比对的最高得分的结果"""
        query = None
        with open(self.path) as f, open(outputfile, 'w') as w:
            w.write(f.readline())
            for line in f:
                new_query = line.split()[5]
                if new_query == query:
                    continue
                else:
                    query = new_query
                    w.write(line)
        table = BlastTableFile()
        table.set_path(outputfile)
        return table

    def fliter_evalue(self, evalue, outputfile):
        """按照设定evalue值进行筛选"""
        if not 1 > evalue >= 0:
            raise ValueError('设定的evalue值不再限制范围内[0,1)')
        with open(self.path) as f, open(outputfile, 'w') as w:
            w.write(f.readline())
            for line in f:
                e = float(line.split()[1])
                if e >= evalue:
                    continue
                else:
                    w.write(line)
        table = BlastTableFile()
        table.set_path(outputfile)
        return table

    def sub_blast_table(self, genes, new_fp):
        """根据提供的基因列表，查找table中的查询序列，生成新的table"""
        genes = dict(zip(genes, xrange(len(genes))))
        with open(self.path, "rb") as f, open("tmp.xls", "wb") as w:
            lines = f.readlines()
            w.write(lines[0])
            for line in lines[1:]:
                item = line.strip().split("\t")
                query_id = item[5]
                if query_id in genes:
                    w.write(line)
        with open("tmp.xls", "rb") as f, open(new_fp, "wb") as w:
            lines = f.readlines()
            w.writelines(lines)


if __name__ == '__main__':  # for test
    a = BlastTableFile()
    # a.set_path('C:\\Users\\sheng.he.MAJORBIO\\Desktop\\annotation\\annotation\\NR\\transcript.fa_vs_nr.blasttable.xls')
    a.set_path('C:\\Users\\sheng.he.MAJORBIO\\Desktop\\blast_test_1.xls')
    a.check()
    a.get_info()
    print(a.prop['header'])
    print(a.prop['count'])
    print(a.prop['query_num'])
    a.retain_highest_score('C:\\Users\\sheng.he.MAJORBIO\\Desktop\\test_1.xls')
    a.fliter_evalue(0.00001, 'C:\\Users\\sheng.he.MAJORBIO\\Desktop\\test_2.xls')
