# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.iofile import File
import os
import copy
from biocluster.core.exceptions import FileError


class DistanceMatrixFile(File):
    """
    定义DistanceMatrix文件
    """
    METHOD = ['abund_jaccard', 'binary_chisq', 'binary_chord',
              'binary_euclidean', 'binary_hamming', 'binary_jaccard',
              'inary_lennon', 'binary_ochiai',
              'binary_pearson', 'binary_sorensen_dice',
              'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran',
              'canberra', 'chisq', 'chord', 'euclidean', 'gower',
              'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
              'pearson', 'soergel', 'spearman_approx', 'specprof',
              'unifrac',
              'unweighted_unifrac', 'unweighted_unifrac_full_tree',
              'weighted_normalized_unifrac', 'weighted_unifrac']
    # 现有的矩阵计算类型

    def __init__(self):
        super(DistanceMatrixFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(DistanceMatrixFile, self).get_info()
        distancematrixinfo = self.get_matrix_info()
        self.set_property('method', distancematrixinfo[0])
        self.set_property('samp_numb', len(distancematrixinfo[1]))
        self.set_property('samp_list', distancematrixinfo[1])

    def get_matrix_info(self):
        """
        获取并返回矩阵信息

        :return method,sample_list:  返回矩阵计算方法和样品列表
        """
        method = self.get_method()
        tempfile = open(self.prop['path'])
        sample_list = tempfile.readline().strip('\n').split('\t')[1:]
        tempfile.close()
        return method, sample_list

    def get_method(self):
        """
        根据提供文件名获取开头的矩阵计算方法

        :return method:  返回矩阵计算方法
        """
        method = os.path.basename(self.prop['path'])
        method = method.split('_')
        if len(method) < 2:
            method = 'unknown_method'
            return method
        if method[0] in DistanceMatrixFile.METHOD:  # 方法为一个单词组成
            method = method[0]
        elif '_'.join([method[0], method[1]]) in DistanceMatrixFile.METHOD:  # 方法为两个单词组成
            method = '_'.join([method[0], method[1]])
        else:
            method = 'unknown_method'
        return method

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        if super(DistanceMatrixFile, self).check():
            # 父类check方法检查文件路径是否设置，文件是否存在，文件是否为空
            self.get_info()
            dist_dict = dict()
            all_values = []
            with open(self.path, 'r') as f:
                head = f.readline().rstrip().split('\t')
                head_len = len(head)
                head = head[1:]
                for line in f:
                    all_nums = line.rstrip().split('\t')
                    if len(all_nums) == head_len:
                        pass
                    else:
                        raise FileError('矩阵每行数据量格式不正确', code="42700201")
                    values = dict(zip(head, all_nums[1:]))
                    all_values.extend(all_nums[1:])
                    dist_dict[all_nums[0]] = values
                for samp1 in head:
                    for samp2 in head:
                        # print dist_dict[samp2][samp1], dist_dict[samp1][samp2]
                        if dist_dict[samp1][samp2] != dist_dict[samp2][samp1]:
                            raise FileError('矩阵数据不对称', code="42700202")
                all_values = [float(i) for i in all_values]
                all_plus = sum(all_values)
                if all_plus == 0 and len(all_values) > 1:  # 只有一个样本时距离就为零，不做处理，但是后续不能做任何分析
                    raise FileError('所有距离矩阵值全部为零', code="42700203")
            if len(self.prop['samp_list']) != len(set(self.prop['samp_list'])):
                raise FileError('存在重复的样本名', code="42700204")
            return True

    def choose(self, sample_list=[], Except=False, path='unknown'):
        """
        选择部分样本生成新的矩阵

        :param samp_list:  一个样品名称列表
        :param except:  除去samp_list列表的样本
        :return newmatrix:  返回一个新的矩阵DistanceMatrixFile对象
        """
        if Except:
            keep_samp = copy.deepcopy(self.prop['samp_list'])
            for i in sample_list:  #
                keep_samp.remove(i)
        else:
            keep_samp = sample_list
        if path == 'unknown':
            filename = os.path.basename(self.prop['path']).split('.')
            dirname = os.path.dirname(self.prop['path'])
            if len(filename) > 1:
                path = dirname + '/' + filename[0] + '_new.' + filename[-1]
            else:
                path = dirname + '/' + filename[0] + '_new'
        self.create_new(keep_samp, path)
        newmatrix = DistanceMatrixFile()
        newmatrix.set_path(path)
        return newmatrix

    def create_new(self, samp_list, path):
        """
        根据样品名和路径，创建一个新的矩阵文件

        :param samp_list:  样品名列表
        :param path:  新文件路径
        """
        dist_dict = dict()
        head = list()
        with open(self.path, 'r') as f:
            head = f.readline().rstrip().split('\t')
            head = head[1:]
            for line in f:
                all_nums = line.rstrip().split('\t')
                values = dict(zip(head, all_nums[1:]))
                dist_dict[all_nums[0]] = values
        newfile = open(path, 'wb')
        newfile.write('\t' + '\t'.join(samp_list) + '\n')
        for samp in samp_list:
            line_list = [samp]
            for samp_col in samp_list:
                line_list.append(dist_dict[samp][samp_col])
            newfile.write('\t'.join(line_list) + '\n')
        newfile.close()

if __name__ == '__main__':
    path = "C:\\Users\\sheng.he.MAJORBIO\\Desktop\\dis.txt"
    a = DistanceMatrixFile()
    a.set_path(path)
    a.check()
