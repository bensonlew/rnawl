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
import subprocess


def create_r(otu_file, group_file, output_dir, one_group):
    """
    生成可以运行的R脚本

    :param otu_file: 输出文件夹
    :param group_file: 分组文件
    :param output_dir: 输出文件夹
    """
    output_dir = output_dir.rstrip('\\')
    output_dir = output_dir.rstrip('/')
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    lab, group = check_group_lab(one_group, group_file)
    f = Template(filename=this_file_dir + '/plsda.r')
    content_r = f.render(otu_file=otu_file, env_file=group, output_dir=output_dir,
                         group_name=lab)
    tempr = open(output_dir + '/temp_r.R', 'w')
    tempr.writelines([i.rstrip() + '\n' for i in content_r.split('\r\n')])
    tempr.close()


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


def plsda(otu_file, group_file, output_dir, one_group):
    create_r(otu_file=otu_file, group_file=group_file, output_dir=output_dir, one_group=one_group)
    return run_r_script(output_dir + '/temp_r.R', delscript=False)

# plsda("C:\\Users\\sheng.he.MAJORBIO\\Desktop\\Plsda\\otu.new.xls", "C:\\Users\\sheng.he.MAJORBIO\\Desktop\\Plsda\\map.txt",  "C:\\Users\\sheng.he.MAJORBIO\\Desktop", 'group')
