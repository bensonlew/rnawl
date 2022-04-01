#!/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "guhaidong"
import os, re


def mask_taxon(abund, new_abund):
    """
    将注释信息（丰度表第一列信息转换为"name+数字"形式）
    :param abund: 输入的原始丰度文件
    :param new_abund: 输出掩蔽掉原有注释信息的丰度文件
    :return:taxon_to_name[new_name] = name，并返回文件taxon_to_name.xls
    """
    taxon_to_name = {}
    out_path = os.path.dirname("new_abund")
    taxon_to_name_file = os.path.join(out_path, "taxon_to_name.xls")
    with open(abund, "r") as f, open(new_abund, "w") as w, open(taxon_to_name_file, "w") as nf:
        first_line = f.readline()
        w.write(first_line)
        n = 1
        for line in f:
            line = line.split("\t")
            name = line[0]
            new_name = "name" + str(n)
            nf.write(new_name + "\t" + name + "\n")
            taxon_to_name[new_name] = name
            n += 1
            new_line = new_name + "\t" + "\t".join(line[1:])
            w.write(new_line)
    return taxon_to_name


def mask_env(env, new_env):
    """
    将环境因子信息（env表第一行信息转换为"name+数字"形式）
    :param env: 输入的环境因子文件
    :param new_env: 输出掩蔽掉原有环境因子
    :return:env_to_name[new_name] = name，并返回文件env_to_name.xls
    """
    env_to_name = {}
    out_path = os.path.dirname("new_env")
    env_to_name_file = os.path.join(out_path, "env_to_name.xls")
    with open(env, "r") as f, open(new_env, "w") as w, open(env_to_name_file, "w") as nf:
        first_line = f.readline().rstrip("\n").split("\t")
        w.write(first_line[0])
        n = 1
        for i in first_line[1:]:
            name = i
            new_name = "name" + str(n)
            nf.write(new_name + "\t" + name + "\n")
            env_to_name[new_name] = name
            n += 1
            new = "\t" + new_name
            w.write(new)
        w.write("\n")
        for line in f:
            w.write(line)
    return env_to_name
