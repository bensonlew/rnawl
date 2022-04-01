#!/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "zhujuan"
from mako.template import Template
import os


this_file_dir = os.path.dirname(os.path.realpath(__file__))


def env_vif(abund_table, env_table, viflim, method, output_dir):
    """
    生成并运行R脚本，进行vif方差膨胀因子分析
    :param abund_table: 输入丰度文件
    :param env_table: 输入环境因子文件
    :param viflim: vif limit to filter[10-20 integer,default:10]
    :param method: rda|cca; if no argument it will choose by DCA result:DCA1>=3.5,CCA;DCA1<3.5,RDA
    :param output_dir: 输出文件存放目录
    """
    f = Template(filename=this_file_dir + '/env_vif.r', strict_undefined=True)
    mul_test = f.render(abund_table=abund_table, env_table=env_table, viflim=viflim, method=method, output_dir=output_dir)
    with open("run_env_vif.r", 'w') as rfile:
        rfile.write("%s" % mul_test)
