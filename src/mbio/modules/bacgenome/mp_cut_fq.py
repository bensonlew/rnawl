# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
class MpCutFqModule(Module):
    """
    微生物基因组组装MP文库截取
    author: gaohao
    last_modify: 2018.03.25
    """
    def __init__(self, work_id):
        super(MpCutFqModule, self).__init__(work_id)
        options = [
            {"name": "fq_dir", "type": "infile", "format": "sequence.fastq_dir"},#clean data 的list文件
            {'name': 'sample_name', "type": "string"},  # 样本名
            {'name': "sample_info", "type": "infile", "format": "sequence.barcode_info"},  # 样本信息表
            {"name": "mp_list", "type": "outfile", "format": "meta.otu.otu_table"}  # 输出目录
        ]
        self.sample_path ={}
        self.sample ={}
        self.modules = []
        self.sample_info = {}
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fq_dir').is_set:
            raise OptionError('必须输入clean data文件夹', code="21401701")
        if not self.option('sample_name'):
            raise OptionError('必须输入sample_name样品名称！', code="21401702")
        if not self.option('sample_info').is_set:
            raise OptionError('必须输入sample_info文件！', code="21401703")

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run_mp_cut(self):
        self.sample_path=self.get_info(self.option('fq_dir').prop['path'] + '/list.txt')
        base_dir = self.option("fq_dir").prop['path']
        n = 0
        for sample in self.sample_path:
            if re.search(r'MP',sample):
                self.mp_cut = self.add_tool('sequence.fastq_mp_cut')
                self.step.add_steps('mp_cut{}'.format(n))
                opts = {
                    "fastq1": base_dir + '/' + self.sample_path[sample]['l'],
                    "fastq2": base_dir + '/' + self.sample_path[sample]['r'],
                    "sample_name": sample,
                }
                self.mp_cut.set_options(opts)
                step = getattr(self.step, 'mp_cut{}'.format(n))
                step.start()
                self.step.update()
                self.mp_cut.on('end', self.finish_update, 'mp_cut{}'.format(n))
                self.modules.append(self.mp_cut)
                n += 1
        self.logger.info(self.modules)
        self.on_rely(self.modules, self.set_output)
        self.step.update()
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(MpCutFqModule, self).run()
        self.get_info(self.option('fq_dir').prop['path'] + '/list.txt')
        time.sleep(2)
        self.run_mp_cut()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        allfiles = os.listdir(dirpath)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for module in self.modules:
            self.linkdir(module.output_dir,self.output_dir)
        time.sleep(2)
        path = self.output_dir + "/list.txt"
        with open(path, 'w') as f:
            files = os.listdir(self.output_dir)
            for file in files:
                if re.search(r'.fq', file):
                    tmp = file.split('.')
                    if tmp[-2] == '1':
                        f.write(file + '\t' + tmp[0] + '\t' + 'l' + '\n')
                    elif tmp[-2] == '2':
                        f.write(file + '\t' + tmp[0] + '\t' + 'r' + '\n')
        f.close()
        time.sleep(2)
        self.get_list()
        file = self.output_dir + '/MP_list.txt'
        self.option('mp_list').set_path(file)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MpCutFqModule, self).end()

    def get_info(self,file):
        """
        获得样本对应的路径信息，以及样本的
        :return:
        """
        sample ={}
        with open(file,'r')as fr:
            for line in fr:
                tmp = line.strip('\r\n').split('\t')
                if tmp[1] in sample.keys():
                    sample[tmp[1]][tmp[2]]=tmp[0]
                else:
                    sample[tmp[1]]={tmp[2]:tmp[0]}
        return sample

    def get_list(self):
        """
        得到生成config索要的输入文件list
        :return:
        """
        with open(self.option("sample_info").prop['path'])as f:
            lines = f.readlines()
            for line in lines[0:]:
                tmp = line.strip().split('\t')
                self.sample_info[tmp[0]] = tmp[1] + '\t' + tmp[2]# 样本的插入片段长度
        file = self.output_dir + '/MP_list.txt'
        with open(file, 'w') as f:
            for sample in self.sample_info:
                if re.search(r'MP', sample):
                    self.sample_mp = self.get_info(self.output_dir + "/list.txt")
                    if sample in self.sample_mp.keys():
                        f.write(sample + '\t' +  self.output_dir + '/' + self.sample_mp[sample]['l'] + '\t' + self.output_dir + '/' + self.sample_mp[sample]['r'] + '\t' + self.sample_info[sample] + '\n')
        f.close()