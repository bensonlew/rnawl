'''
random forest tool constants
'''
import os
import numpy as np


out_num = 4
data_num = 8

n_jobs = 4

type_num = [np.dtype('int'), np.dtype('float')]
type_o = [np.dtype('O')]

data_sep = '\t'

mode, class_col, filter = 'median', 0, .5

n_cpu = 10

db_public = 'table'
db_bar = 'bar'
db_curve = 'curve'
details = 'detail'
db_name_sep = '_'

table_id = '_id'

table_status = 'table_status'
table_set = u'$set'

status_start = 'failed'
status_end = 'end'

sh = 'library/rftool/bin/rf.sh'

step_rf = 'step_rf'
input = 'input'

variable = 'variable'

tables = [
    dict(name='test_probability_predict',
         columns=None,  # ['id', 'predict', 'features'],
         data_path='test_probability_predict_rf.txt',
         db=db_public, ),
    dict(name='feature_importances',
         columns=['feature', 'importance', ],
         data_path='feature_importances_rf.txt',
         db=db_public, ),
    dict(name='feature_importances',
         columns=['feature', 'importance', ],
         data_path='feature_importances_rf.txt',
         db=db_bar, ),
    dict(name='best_feature_numbers',
         columns=['feature_number', 'error', ],
         data_path='rfe_scores.txt',
         db=db_curve, ), ]


def get_header(file):
    with open(file) as f:
        return f.readline().split(data_sep)


def read0(file, header_only=False):
    """
    get data header or data body

    Parameters
    ----------
    file : str
        file name

    Returns
    -------
    list of list of str, or None if read if only has data header

    """
    data = []
    with open(file) as f:
        if header_only:
            header = f.readline().split(data_sep)
            data.append[header]
            return data
        
        _ = f.readline()
        for line in f:
            line_split = line.split(data_sep)
            data.append(line_split)
    if data  == []:
        return None
    return data
    
    
def read(file):
    """
    get data header or data body

    Parameters
    ----------
    file : str
        file name

    Returns
    -------
    list of list of str, or None if read if only has data header

    """
    data = []
    with open(file) as f:        
        _ = f.readline()  # has header default
        for line in f:
            line_split = line.split(data_sep)
            data.append(line_split)
    if data  == []:
        return None
    return data
