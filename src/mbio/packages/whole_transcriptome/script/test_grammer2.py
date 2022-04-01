# -*- coding: utf-8 -*-

from __future__ import print_function

global oid


def m(x):
    print('b2: {}'.format(id(x) == oid))
    x += 1
    print('b3: {}'.format(id(x) == oid))
    return x - 1


if __name__ == '__main__':
    a = 1
    oid = id(a)
    print('b1: {}'.format(id(a) == oid))
    b = m(a)
    print('b4: {}'.format(id(b) == oid))
