# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError

class GoLevel2File(File):
    def __init__(self):
        super(GoLevel2File, self).__init__()

    def get_info(self):
        super(GoLevel2File, self).get_info()
        go2level_info = self.get_table_info()
        self.set_property('header', go2level_info[1])
        self.set_property('count', go2level_info[0])

    def get_table_info(self):
        with open(self.path) as f:
            go_info = f.read().split('\n')
            header = go_info[0]
            if len(header.split('\t')) != 6:
                raise FileError('Error in file header, default header:\nterm_type\tterm\tGO\tnumber\tpercent\tsequence')
            count = 0
            for record in go_info[1:]:
                if record != '':
                    count += 1
                    record_info = record.split('\t')
                    if len(record_info) != 6:
                        raise FileError('There is a disordered line in the file')
                    if record_info[0] not in ['biological_process', 'cellular_component', 'molecular_function']:
                        raise FileError(
                            'There is an incorrect first level classification in the file -> {}'.format(record_info[0])
                        )
        return count, header

    def check(self):
        if super(GoLevel2File, self).check():
            self.get_info()
            return True

    def get_gene(self):
        with open(self.path) as f:
            f.readline()
            gene_list = []
            for line in f:
                line = line.strip('\n').split('\t')
                genes = line[-1].split(';')
                for i in genes:
                    gene_list.append(i.split('(')[0])
        return gene_list
