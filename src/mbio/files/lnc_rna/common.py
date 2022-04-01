# -*- coding: utf-8 -*-
# __author__ = "liubinxu"
import os

from biocluster.iofile import File

class CommonFile(File):
    def __init__(self):
        super(CommonFile, self).__init__()

    def check(self):
        super(CommonFile, self).check()

    def hard_link(self, new_path):
        if os.path.isdir(os.path.dirname(new_path)):
            cmd = 'ln {old} {new}'.format(old=self.path, new=new_path)
            return_code = os.system(cmd)
            if return_code not in (0, '0'):
                cmd = 'cp {old} {new}'.format(old=self.path, new=new_path)
                os.system(cmd)
        else:
            raise Exception('文件创建link失败，请确认文件路径是否正确')

    def basename(self):
        return os.path.basename(self.path)
