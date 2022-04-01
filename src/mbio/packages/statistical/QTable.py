# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
import os
import math


class QTable:

    def __init__(self, preferences):
        CVs = ['0.001',
         '0.01',
         '0.02',
         '0.05',
         '0.1']
        tables = ['tukeyQ_001.txt',
         'tukeyQ_01.txt',
         'tukeyQ_02.txt',
         'tukeyQ_05.txt',
         'tukeyQ_10.txt']
        self.qTables = {}
        file_path = os.path.dirname(os.path.realpath(__file__))
        for i in xrange(0, len(tables)):
            dataPath = '%s/tukey_table/%s' % (file_path, tables[i])
            fin = open(dataPath, 'U')
            data = fin.readlines()
            fin.close()
            qTable = {}
            for j in xrange(1, len(data)):
                lineSplit = data[j].split('\t')
                values = [float(x) for x in lineSplit[1:] ]
                qTable[j + 1] = values

            self.qTables[CVs[i]] = qTable

    def cv(self, alpha, k, df):
        if k < 2 or df < 2:
            return 1000.0
        if k > 29:
            k = 29
        if df > 200:
            df = 200
        alphaStr = '%.3g' % alpha
        return self.qTables[alphaStr][df][k - 2]

    def cvInterpolate(self, alpha, k, df):
        dfLower = int(math.floor(df))
        dfUpper = int(math.ceil(df))
        cvLower = self.cv(alpha, k, dfLower)
        cvUpper = self.cv(alpha, k, dfUpper)
        cv = cvLower + (df - dfLower) * (cvUpper - cvLower)
        return cv
