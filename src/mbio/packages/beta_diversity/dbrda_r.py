# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
# __version__ = 'v1.0'
# __last_modified__ = '20151117'
"""

"""

import os
import platform
from mako.template import Template
from biocluster.config import Config


def create_r_new(otu_file, env_file, output_dir, distance_algorithm, sqrt_dist='FALSE', add_dist='FALSE'):
    """
    生成可以运行的R脚本

    :rtype : object
    :param otu_file: 输出文件夹
    :param env_file: 环境因子文件
    :param output_dir: 输出文件夹
    :param distance_algorithm: 距离计算方法（此处计算方法是R中的vegdist方法）
    """
    output_dir = output_dir.rstrip('\\')
    output_dir = output_dir.rstrip('/')
    with open(env_file) as env:
        env_formula = ' + '.join(env.readline().strip().split('\t')[1:])
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    f = Template(filename=this_file_dir + '/db_rda.r')
    content_r = f.render(otu_file=otu_file, env_file=env_file, output_dir=output_dir,
                         distance_algorithm=distance_algorithm, env_formula=env_formula,
                         sqrt_dist=sqrt_dist, add=add_dist)
    tempr = open(output_dir + '/temp_r.R', 'w')
    tempr.writelines([i.rstrip() + '\n' for i in content_r.split('\r\n')])
    tempr.close()


def create_r(dist_matrix, env_file, species_table, output_dir, sqrt_dist='FALSE', add_dist='FALSE'):
    """
    生成可以运行的R脚本,输入文件为距离矩阵

    :param dist_matrix: 输入矩阵
    :param env_file: 输入环境因子文件
    :param species_table 输入物种丰度表格  # add 'species_table' by zhujuan for add db_rda_species.xls 20171010
    :param output_dir: 输出文件夹
    """
    output_dir = output_dir.rstrip('\\')
    output_dir = output_dir.rstrip('/')
    with open(env_file) as env:
        env_formula = ' + '.join(env.readline().strip().split('\t')[1:])
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    f = Template(filename=this_file_dir + '/db_rda_dist.r')
    print f
    print "this_file_dir:", this_file_dir
    content_r = f.render(dis_matrix=dist_matrix, env_file=env_file, species_table=species_table, output_dir=output_dir,
                         env_formula=env_formula, sqrt_dist=sqrt_dist, add=add_dist)
    tempr = open(output_dir + '/temp_r_dist.R', 'w')
    tempr.writelines([i.rstrip() + '\n' for i in content_r.split('\r\n')])
    tempr.close()


def run_r_script(script, delscript=True):
    """
    分平台运行R脚本，运行完成脚本会被删除
    :param script:R脚本路径（路径使用斜杠，不要使用反斜杠）
    """
    if platform.system() == 'Windows':
        os.system('R CMD BATCH --vanilla --slave %s ' % script)
    elif platform.system() == 'Linux':
        os.system('%s/program/R-3.3.1/bin/Rscript %s' % (Config().SOFTWARE_DIR, script))
    else:
        pass
    if delscript:
        os.remove(script)
        if os.path.exists(os.path.dirname(script) + '/temp_r.Rout'):
            os.remove(os.path.dirname(script) + '/temp_r.Rout')


def db_rda_dist(dis_matrix, env, species_table, output_dir, sqrt_dist='FALSE', add_dist='FALSE'):
    """
    输入距离矩阵，物种丰度文件，分组信息，输出文件夹，进行db_rda分析
    :pararm return: 成功完成返回‘0’
    """
    create_r(dis_matrix, env, species_table, output_dir, sqrt_dist, add_dist)
    script = output_dir + '/temp_r_dist.R'
    run_r_script(script, delscript=False)
    return 0


def db_rda_new(otu_file, env_file, output_dir, distance_algorithm='bray', sqrt_dist='FALSE', add_dist='FASLE'):
    """
    """
    create_r_new(otu_file, env_file, output_dir, distance_algorithm, sqrt_dist, add_dist)
    script = output_dir + '/temp_r.R'
    run_r_script(script, delscript=False)
    return 0
