# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/25 13:33

import re
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fasta import FastaFile
import os


class MapspliceFaDirFile(Directory):
    """
    fasta 单序列/文件格式
    """

    def __init__(self):
        super(MapspliceFaDirFile, self).__init__()

    def check(self):
        if self.check_content():
            return True

    def check_content(self):
        super(MapspliceFaDirFile, self).check()
        f_list = [os.path.join(self.path, f) for f in os.listdir(self.path) if re.match(r'^(\S+)\.fa$', f.strip())]
        for fa in f_list:
            file_name = re.match(r'^(\S+)\.fa$', os.path.basename(fa.strip())).group(1)
            fasta_obj = FastaFile()
            fasta_obj.set_path(fa)
            if fasta_obj.check():
                seq_number = int(fasta_obj.prop['seq_number'].strip())
                if seq_number == 1:
                    fr = open(fa)
                    n = 0
                    while n < 1:
                        line = fr.readline()
                        if not line.strip():
                            continue
                        m = re.match(r'^>(\S+)\s*$', line.strip())
                        if m:
                            n = n + 1
                            seq_name = m.group(1)
                            if file_name == seq_name:
                                print '{}的文件名与其包含的序列名一致'.format(fa, file_name, seq_name)
                            else:
                                raise FileError('{}的文件名{} 与其包含的序列名{}不一致'.format(fa, file_name, seq_name))
                else:
                    raise FileError('{}含有的序列数目不是1条，为{}条'.format(fa, seq_number))
        return True

# if __name__ == '__main__':
#     test_dir = '/mnt/ilustre/users/sanger-dev/workspace/20170125/Single_mapsplice_split_linfang/MapspliceSplitFa/output'
#     obj = MapspliceFaDirFile()
#     obj.set_path(test_dir)
#     obj.check()
