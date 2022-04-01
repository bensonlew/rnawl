# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
# __version__ = 'v1.0'
# __last_modified__ = '20161012'
"""
adonis分析脚本， 提供adonis分析需要的 距离矩阵，group分组文件，选定的分组文件列名， permutations置换次数值
"""

import os
import platform
from mako.template import Template
from biocluster.config import Config
import subprocess


def create_r(dis_matrix_file, group_file, output_dir, one_group, permutations=999):
    """
    生成可以运行的R脚本

    :param otu_file: 输出文件夹
    :param group_file: 分组文件
    :param output_dir: 输出文件夹
    """
    output_dir = output_dir.rstrip('\\')
    output_dir = output_dir.rstrip('/')
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    f = Template(filename=this_file_dir + '/adonis.r')
    lab, group = check_group_lab(one_group, group_file)
    content_r = f.render(dist_matrix=dis_matrix_file, group_file=group, output_dir=output_dir,
                         flied=lab, permutations=permutations)
    tempr = open(output_dir + '/temp_r.R', 'w')
    tempr.writelines([i.rstrip() + '\n' for i in content_r.split('\r\n')])
    tempr.close()
    return output_dir + '/temp_r.R'


def check_group_lab(lab, group, new_name='Group', new_group=None):
    if lab[0].isdigit():
        if not new_group:
            new_group = group + '.temp'

        def rplace_group(matched):
            return '\t' + matched.groups()[0] + '\n'
        with open(group) as f:
            data = f.readlines()
            names = data[0].rstrip().split('\t')
            new_headers = []
            for i in names:
                if i == lab:
                    new_headers.append(new_name)
                else:
                    new_headers.append(i)
            data[0] = '\t'.join(new_headers) + '\n'
        with open(new_group, 'w') as w:
            w.writelines(data)
        return new_name, new_group
    else:
        return lab, group


def run_r_script(script, delscript=True):
    """
    分平台运行R脚本，运行完成脚本会被删除
    :param script:R脚本路径（路径使用斜杠，不要使用反斜杠）
    """
    if platform.system() == 'Windows':
        cmd = 'R CMD BATCH --vanilla --slave %s ' % (script)
        # os.system('R CMD BATCH --vanilla --slave %s ' % (script))
    elif platform.system() == 'Linux':
        cmd = '%s/program/R-3.3.1/bin/Rscript %s' % (Config().SOFTWARE_DIR, script)
        # os.system('%s/R-3.2.2/bin/Rscript %s' % (Config().SOFTWARE_DIR, script))
    else:
        pass
    try:
        subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        return "error, info:{}".format(e)
    if delscript:
        os.remove(script)
        if os.path.exists(os.path.dirname(script) + '/temp_r.Rout'):
            os.remove(os.path.dirname(script) + '/temp_r.Rout')
    return 0


def adonis(dist_matrix, group_file, output_dir, one_group, permutations):
    """
    运行主函数
    adonis分析产生两个结果文件在结果输出目录 分别为adonis.txt(结果直接print)和adonis_format.txt(结果格式整理文件)

    :param dis_matrix:距离矩阵
    :param group_file:分组文件
    :param output_dir:输出文件夹
    :param one_group:分组文件选定一个分组方案名
    :param permutations:分析置换次数
    """
    R_file = create_r(dis_matrix_file=dist_matrix, group_file=group_file, output_dir=output_dir, one_group=one_group, permutations=permutations)
    return run_r_script(R_file, delscript=False)

# adonis("C:\\\\Users\\\\sheng.he.MAJORBIO\\\\Desktop\\\\bray_curtis_otu_table.txt", "C:\\\\Users\\\\sheng.he.MAJORBIO\\\\Desktop\\\\env_table_group.txt",  "C:\\\\Users\\\\sheng.he.MAJORBIO\\\\Desktop", 'GROUP1')
