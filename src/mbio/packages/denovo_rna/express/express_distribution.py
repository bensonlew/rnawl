#!/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "qiuping"


from mako.template import Template
from biocluster.config import Config
import os


def distribution(rfile, input_matrix, outputfile, filename):
    """
    生成表达量分布图的数据
    :param input_matrix: 表达量矩阵
    :param outputfile: 输出文件路径
    :param filename: gene或transcript
    """
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    f = Template(filename=this_file_dir + '/express_distribution.r')
    distribution = f.render(input_matrix=input_matrix, outputfile=outputfile, filename = filename)
    with open(rfile, 'w') as rf:
        rf.write("%s" % distribution)

if __name__ == "__main__":
    distribution("express_distribution.r","/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tofiles/gene_fpkm","/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tofiles/",filename = "gene")
    print 'end'
