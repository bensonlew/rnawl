#!/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "zhangyitong"


from mako.template import Template
from biocluster.config import Config
import os


this_file_dir = os.path.dirname(os.path.realpath(__file__))


def preprocess(exp_file, group_file, na_ratio, method, specific, group_threshold, fill_type, out_prefix,
               interactive=None):
    '''
    生成并运行R脚本，进行数据预处理，包括定量值去除，组特异性判断，缺失值填充
    :param exp_file: 原始表达谱
    :param group_file: 分组文件
    :param na_ratio: 定量值去除threshold
    :param method: 缺失值填充方法
    :param specific: 是否进行组特异性判断
    :param group_threshold: 组特异性判断threshold
    :param fill_type: 按组或全部样本填充
    :param out_prefix: 输出文件
    :param interactive: 是否为交互分析
    :return:
    '''
    f = Template(filename=this_file_dir + '/dia_preprocess.r', input_encoding='utf-8', output_encoding='utf-8')
    exp_preprocess = f.render(exp=exp_file, group=group_file, na_ratio=na_ratio, method=method, specific=specific,
                              specific_threshold=group_threshold, fill_type=fill_type, outprefix=out_prefix,
                              interactive=interactive)
    with open("preprocess_%s.r" % method, 'w') as rfile:
        rfile.write("%s" % exp_preprocess)

