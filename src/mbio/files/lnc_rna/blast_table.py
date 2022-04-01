# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

from biocluster.iofile import File
from mbio.packages.align.blast.xml2table import default_header
from biocluster.core.exceptions import FileError

class BlastTableFile(File):
    '''
    Define a table format for blast+ result of outfmt 6
    The file actually come not from blast, but through XML transformation
    '''

    def __init__(self):
        super(BlastTableFile, self).__init__()

    def get_info(self):
        '''
        Get file properties
        '''
        super(BlastTableFile, self).get_info()
        table_info = self.get_table_info()
        self.set_property('count', table_info[0])
        self.set_property('header', table_info[1])
        self.set_property('query_num', table_info[2])
        self.set_property('query_list', table_info[3])

    def get_table_info(self):
        '''
        Get information for the blast result table
        '''
        with open(self.path) as f:
            header = f.readline().strip().split()
            if header != default_header:
                raise FileError(
                    '默认表头与文件表头不一致：\n默认表头：%s\n文件表头：%s', variables = (
                        '\t'.join(default_header), '\t'.join(header)), code = '43700501')
            count = 0
            query_count = 0
            query_list = []
            flag = None
            for line in f:
                count += 1
                line_sp = line.split('\t')
                if len(line_sp) != 16:
                    raise FileError('文件中存在格式不整齐的行：\n%s', variables = (line), code = '43700502')
                if flag == line_sp[5]:
                    pass
                else:
                    query_list.append(line_sp[5])
                    query_count += 1
                    flag = line_sp[5]
        return count, default_header, query_count, query_list

    def check(self):
        '''
        Check whether the documents meet the requirements
        '''
        if super(BlastTableFile, self).check():
            self.get_info()
            return True

    def retain_highest_score(self, outputfile):
        '''
        Retain the result of the highest score for each query comparison
        '''
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
        '''
        Filter by setting evalue
        '''
        if not 1 > evalue >= 0:
            raise FileError('设定的evalue值不再限制范围内[0,1)')
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
        '''
        According to the gene list provided, search the query sequence in the table to generate a new table
        '''
        genes = dict(zip(genes, xrange(len(genes))))
        with open(self.path, 'rb') as f, open('tmp.xls', 'wb') as w:
            lines = f.readlines()
            w.write(lines[0])
            for line in lines[1:]:
                item = line.strip().split('\t')
                query_id = item[5]
                if query_id in genes:
                    w.write(line)
        with open('tmp.xls', 'rb') as f, open(new_fp, 'wb') as w:
            lines = f.readlines()
            w.writelines(lines)
