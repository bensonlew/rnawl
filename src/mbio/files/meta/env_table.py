# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.iofile import File
import os
import copy
from biocluster.core.exceptions import FileError


class EnvTableFile(File):
    """
    定义envtable环境因子表文件
    """

    def __init__(self):
        super(EnvTableFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(EnvTableFile, self).get_info()
        envinfo = self.get_env_info()
        self.set_property('env', envinfo[0])
        self.set_property('sample', envinfo[1])

    def get_env_info(self):
        """
        获取并返回环境因子信息

        :return :
        """
        tempfile = open(self.prop['path'])
        lines = tempfile.readlines()
        env = lines[0].rstrip('\n').split('\t')[1:]
        sample = [i.split('\t')[0] for i in lines[1:]]
        tempfile.close()
        return env, sample

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        if super(EnvTableFile, self).check():
            # 父类check方法检查文件路径是否设置，文件是否存在，文件是否为空
            self.get_info()
            if len(self.prop['sample']) != len(set(self.prop['sample'])):
                raise FileError('存在重复的样本名')
            tempfile = open(self.prop['path'])
            env_num = 0
            for line in tempfile:
                if env_num == 0:
                    env_num = len(line.split('\t'))
                else:
                    if env_num == len(line.split('\t')):
                        pass
                    else:
                        raise FileError('存在数据冗余或者数据缺失')
        return True


    def choose(self, env_list=[], Except=False, path='unknown'):
        """
        选择部分环境因子生成新的环境因子表

        :param env_list:  一个样品名称列表
        :param Except:  除去env_list列表的样本,默认是保留
        :return newenv:  返回一个新的EnvTableFile对象
        """
        if Except:
            keep_env = copy.deepcopy(self.prop['env'])
            for i in env_list:  #
                keep_env.remove(i)
        else:
            keep_env = env_list
        if path == 'unknown':
            filename = os.path.basename(self.prop['path']).split('.')
            dirname = os.path.dirname(self.prop['path'])
            if len(filename) > 1:
                path = dirname + '/' + filename[0] + '_new.' + filename[-1]
            else:
                path = dirname + '/' + filename[0] + '_new'
        self.create_new(keep_env, path)
        newenv = EnvTableFile()
        newenv.set_path(path)
        return newenv

    def create_new(self, env_list, path):
        """
        根据样品名和路径，创建一个新的环境因子表

        :param env_list:  样品名列表
        :param path:  新文件路径
        """
        newfile = open(path, 'w')
        oldfile = open(self.prop['path'])
        oldlines = oldfile.readlines()
        newfile.write('\t' + '\t'.join(env_list) + '\n')
        column = []
        for m in env_list:
            column.append(self.prop['env'].index(m) + 1)
        column.insert(0, 0)

        for line in oldlines[1:]:
            newline = '\t'.join([line.strip().split('\t')[i] for i in column])
            newfile.write(newline + '\n')
        oldfile.close()
        newfile.close()
