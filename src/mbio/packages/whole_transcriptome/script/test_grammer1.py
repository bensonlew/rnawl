# -*- coding: utf-8 -*-

from __future__ import print_function

if __name__ == '__main__':
    x1 = []
    y1 = []
    b1 = id(x1) == id(y1)
    x2 = {}
    y2 = x2
    b2 = id(x2) == id(y2)
    x3 = ()
    y3 = ()
    b3 = id(x3) == id(y3)
    print('b1: {}, b2: {}, b3: {}'.format(b1, b2, b3))
