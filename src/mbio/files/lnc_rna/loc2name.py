# -*_ coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.iofile import File
import os
from biocluster.core.exceptions import FileError

class Loc2nameFile(File):
    def __init__(self):
        super(Loc2nameFile, self).__init__()
        self.names = list()
        self.locations = dict()

    def check(self):
        super(Loc2nameFile, self).check()
        if not os.path.exists(self.prop['path']):
            raise FileError('{} does not exist'.format(self.prop['path']))
        self.obtain()
        self.set_property('names', self.names)
        self.set_property('locations', self.locations)

    def obtain(self):
        for n, line in enumerate(open(self.prop['path'])):
            if line.strip() == '':
                continue
            items = line.strip().split('\t')
            if len(items) != 2:
                raise FileError('find format error in {} at line {}'.format(self.prop['path'], n))
            location = items[0]
            name = items[1]
            if name in self.names:
                raise FileError('find duplicate {} in {} at line {}'.format(name, self.prop['path']), n)
            else:
                self.names.append(name)
                self.locations[name] = location
