#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == zhouxuan

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import shutil

class FungiGenomeModule(Module):
    """
    对数据进行解压统计碱基质量
    """
    def __init__(self, work_id):
        super(FungiGenomeModule, self).__init__(work_id)
        options = [
            {"name": "raw_dir", "type": "infile", "format": "bacgenome.raw_dir"},
            {"name": "sample_name", "type": "string"},
            {"name": "sequence_type", "type": "string"},
            {"name": "dir", "type": "outfile", "format": "sequence.fastq_dir"},
            {"name": "pacbio_dir", "type": "outfile", "format": "bacgenome.pacbio_dir"},
        ]
        self.add_option(options)
        self.samplesi = ""
        self.tools = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option("raw_dir").is_set:
            raise OptionError("必须输入样本文件夹！", code="22100801")
        elif not self.option("sample_name"):
            raise OptionError("必须输入样本名！", code="22100802")
        elif not self.option("sequence_type"):
            raise OptionError("必须read的文库类型！", code="22100803")
        else:
            return True

    def run(self):
        super(FungiGenomeModule, self).run()
        self.samples = self.get_list()
        if self.option("raw_dir").is_set and self.option("sequence_type") in ['PE','PE,MP','MP,PE']:
            self.run_ungiz()
        elif self.option("sequence_type") in ['pacbio','Pacbio','PACBIO']:
            self.set_output()
        elif self.option("sequence_type") in ['PE,pacbio', 'pacbio,PE', 'PE,Pacbio', 'Pacbio,PE']:
            self.run_ungiz()


    def get_list(self):
        sample = {}
        if self.option("sequence_type") in ['PE','PE,MP','MP,PE']:
            raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
            output2 = self.work_dir + "/" + 'sample_info'
            with open(raw_list_path, "rb") as l, open(output2, "w") as file2:
                lines = l.readlines()
                lib_list = []
                raw_list = []
                for line in lines[1:]:
                    line2 = line.strip().split()
                    if len(line2) == 6:
                        lib_type = line2[5] + line2[2]
                        lib_list.append(lib_type)
                        raw_list.append(line2[1])
                        file2.write(lib_type + '\t' + line2[2] + '\t' + line2[3] + '\t' + line2[4] +'\n')
                        if re.search(';', line2[1]):
                            raw_path = line2[1].split(';')
                            for raw in raw_path:
                                files = raw.split(',')
                                for i in range(len(files)):
                                    if i % 2:
                                        type = 'l'
                                    else:
                                        type = 'r'
                                    sample_path = self.option("raw_dir").prop['path'] + "/" + files[i]
                                    if lib_type not in sample:
                                        sample[lib_type] = {type: sample_path}
                                    elif type in sample[lib_type].keys():
                                        sample[lib_type][type] = sample[lib_type][type] + " " + sample_path
                                    else:
                                        sample[lib_type][type] = sample_path
                        else:
                            files = line2[1].split(',')
                            for i in range(len(files)):
                                if i % 2:
                                    type = 'l'
                                else:
                                    type = 'r'
                                sample_path = self.option("raw_dir").prop['path'] + "/" + files[i]
                                if lib_type not in sample:
                                    sample[lib_type] = {type: sample_path}
                                elif type in sample[lib_type].keys():
                                    sample[lib_type][type] = sample[lib_type][type] + " " + sample_path
                                else:
                                    sample[lib_type][type] = sample_path
                    else:
                        self.set_error('list.txt文件格式有误', code="22100801")
        elif self.option("sequence_type") in ['pacbio','Pacbio','PACBIO']:
            raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
            output = self.work_dir + "/" + 'pacbio.rawdata.list'
            with open(raw_list_path, "rb") as l, open(output, "w") as file:
                lines = l.readlines()
                for line in lines[1:]:
                    line2 = line.strip('\r\n').split()
                    if len(line2) == 6:
                        if re.search(',', line2[1]):
                            pacs = line2[1].split(',')
                            new_path = ''
                            for pac in pacs:
                                sample_path = self.option("raw_dir").prop['path'] + "/" + pac
                                new_path = self.work_dir + "/all.pacbio.fq"
                                os.system('cat %s >>%s' % (sample_path, new_path))
                            file.write(
                                line2[0] + '\t' + 'all.pacbio.fq' + '\t-\t' + '-\t' + line2[4] + '\t' + line2[5] + '\n')
                        else:
                            if os.path.exists( self.work_dir + "/all.pacbio.fq"):
                                os.remove( self.work_dir + "/all.pacbio.fq")
                            os.link(self.option("raw_dir").prop['path'] + "/" + line2[1],
                                    self.work_dir + "/all.pacbio.fq")
                            file.write(
                                line2[0] + '\t' + 'all.pacbio.fq' + '\t-\t' + '-\t' + line2[4] + '\t' + line2[5] + '\n')
                    else:
                        self.set_error('list.txt文件格式有误', code="22100802")
        elif self.option("sequence_type") in ['PE,pacbio','pacbio,PE','PE,Pacbio','Pacbio,PE']:
            raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
            output = self.work_dir + "/" + 'pacbio.rawdata.list'
            output2 = self.work_dir + "/" + 'sample_info'
            with open(raw_list_path, "rb") as l, open(output, "w") as file,open(output2,'w') as file2:
                lines = l.readlines()
                for line in lines[1:]:
                    line2 = line.strip('\r\n').split()
                    if len(line2) == 6:
                        if line2[5] in ["Pacbio", "pacbio", "PACBIO"]:
                            if re.search(',', line2[1]):
                                pacs = line2[1].split(',')
                                new_path = ''
                                for pac in pacs:
                                    sample_path = self.option("raw_dir").prop['path'] + "/" + pac
                                    new_path = self.work_dir + "/all.pacbio.fq"
                                    os.system('cat %s >>%s' % (sample_path, new_path))
                                file.write(
                                    line2[0] + '\t' + 'all.pacbio.fq' + '\t-\t' + '-\t' + line2[4] + '\t' + line2[5] + '\n')
                            else:
                                if os.path.exists(self.work_dir + "/all.pacbio.fq"):
                                    os.remove(self.work_dir + "/all.pacbio.fq")
                                os.link(self.option("raw_dir").prop['path'] + "/" + line2[1],self.work_dir + "/all.pacbio.fq")
                                file.write(line2[0] + '\t' + 'all.pacbio.fq' + '\t-\t' + '-\t' + line2[4] + '\t' + line2[5] + '\n')
                        else:
                            lib_type = line2[5] + line2[2]
                            file2.write(lib_type + '\t' + line2[2] + '\t' + line2[3] + '\t' + line2[4] +'\n')
                            if re.search(';', line2[1]):
                                raw_path = line2[1].split(';')
                                for raw in raw_path:
                                    files = raw.split(',')
                                    for i in range(len(files)):
                                        if i % 2:
                                            type = 'l'
                                        else:
                                            type = 'r'
                                        sample_path = self.option("raw_dir").prop['path'] + "/" + files[i]
                                        if lib_type not in sample:
                                            sample[lib_type] = {type: sample_path}
                                        elif type in sample[lib_type].keys():
                                            sample[lib_type][type] = sample[lib_type][type] + " " + sample_path
                                        else:
                                            sample[lib_type][type] = sample_path
                            else:
                                files = line2[1].split(',')
                                for i in range(len(files)):
                                    if i % 2:
                                        type = 'l'
                                    else:
                                        type = 'r'
                                    sample_path = self.option("raw_dir").prop['path'] + "/" + files[i]
                                    if lib_type not in sample:
                                        sample[lib_type] = {type: sample_path}
                                    elif type in sample[lib_type].keys():
                                        sample[lib_type][type] = sample[lib_type][type] + " " + sample_path
                                    else:
                                        sample[lib_type][type] = sample_path
                    else:
                        self.set_error('list.txt文件格式有误', code="22100803")
        return sample

    def run_ungiz(self):
        samples = self.samples
        reslut_path = os.path.join(self.work_dir, "ungiz_dir")
        if not os.path.exists(reslut_path):
            os.mkdir(reslut_path)
        for lib in samples:
            for d in samples[lib]:
                direct = ""
                if d == "r":
                    direct = "2"
                elif d == "l":
                    direct = "1"
                else:
                    self.set_error("序列的方向不对，必须为：l/r", code="22100804")
                gunzip_fastq = self.add_tool('bacgenome.fastq_ungz')
                gunzip_fastq.set_options({
                    "fastq": samples[lib][d],
                    "sample_name": self.option("sample_name"),
                    "direction": direct,
                    "lib_type": lib,
                    "result_path": reslut_path
                })
                self.tools.append(gunzip_fastq)
        if len(self.tools) > 1:
            self.on_rely(self.tools, self.set_output)
        else:
            self.tools[0].on('end', self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        if self.option("sequence_type") in ['PE','PE,MP','MP,PE']:
            if os.path.exists(self.work_dir + "/ungiz_dir"):
                try:
                    self.linkdir(self.work_dir + "/ungiz_dir", "data")
                except Exception, e:
                    self.set_error('解压的结果linkdir时出错%s', variables=(e), code="22100805")
            samples = self.samples
            list = os.path.join(self.output_dir, "data/list.txt")
            with open(list, "wb") as w:
                for sample in samples:
                    for d in samples[sample]:
                        direct = ""
                        if d == "r":
                            direct = "2"
                        if d == "l":
                            direct = "1"
                        fq_name = self.option("sample_name") + "_" + sample + "." + direct + ".fq"
                        w.write(fq_name + '\t' + sample + '\t' + d + '\n')
            self.option('dir', self.output_dir + "/data")
            if os.path.exists(self.work_dir + "/" + 'sample_info'):
                try:
                    if os.path.exists(self.output_dir + "/data/" + 'sample_info'):
                        os.remove(self.output_dir + "/data/" + 'sample_info')
                    os.link(self.work_dir + "/" + 'sample_info', self.output_dir + "/data/" + 'sample_info')
                except Exception, e:
                    self.set_error('sample_info的文件link时出错%s', variables=(e), code="22100806")
            self.end()
        elif self.option("sequence_type") in ['pacbio','Pacbio','PACBIO']:
            if not os.path.exists(self.output_dir + '/pacbio_data'):
                os.mkdir(self.output_dir + '/pacbio_data')
            if os.path.exists(self.output_dir + '/pacbio_data/all.pacbio.fq'):
                os.remove(self.output_dir + '/pacbio_data/all.pacbio.fq')
            os.link(self.work_dir + "/all.pacbio.fq",self.output_dir + '/pacbio_data/all.pacbio.fq')
            if os.path.exists(self.work_dir + "/" + 'pacbio.rawdata.list'):
                try:
                    if os.path.exists(self.output_dir + "/pacbio_data/" + 'pacbio.rawdata.list'):
                        os.remove(self.output_dir + "/pacbio_data/" + 'pacbio.rawdata.list')
                    os.link(self.work_dir + "/" + 'pacbio.rawdata.list',self.output_dir + "/pacbio_data/" + '/pacbio.rawdata.list')
                except Exception, e:
                    self.set_error('pacbio.rawdata.list的文件link时出错%s', variables=(e), code="22100807")
            self.option('pacbio_dir', self.output_dir + "/pacbio_data")
            self.end()
        elif self.option("sequence_type") in ['PE,pacbio','pacbio,PE','PE,Pacbio','Pacbio,PE']:
            if os.path.exists(self.work_dir + "/ungiz_dir"):
                try:
                    self.linkdir(self.work_dir + "/ungiz_dir", "data")
                except Exception, e:
                    self.set_error('解压的结果linkdir时出错%s', variables=(e), code="22100808")
            samples = self.samples
            list = os.path.join(self.output_dir, "data/list.txt")
            with open(list, "wb") as w:
                for sample in samples:
                    for d in samples[sample]:
                        direct = ""
                        if d == "r":
                            direct = "2"
                        if d == "l":
                            direct = "1"
                        fq_name = self.option("sample_name") + "_" + sample + "." + direct + ".fq"
                        w.write(fq_name + '\t' + sample + '\t' + d + '\n')
            self.option('dir', self.output_dir + "/data")
            if os.path.exists(self.work_dir + "/" + 'sample_info'):
                try:
                    if os.path.exists(self.output_dir + "/data/" + 'sample_info'):
                        os.remove(self.output_dir + "/data/" + 'sample_info')
                    os.link(self.work_dir + "/" + 'sample_info', self.output_dir + "/data/" + 'sample_info')
                except Exception, e:
                    self.set_error('sample_info的文件link时出错%s', variables=(e), code="22100809")
            if os.path.exists(self.output_dir + '/pacbio_data'):
                shutil.rmtree(self.output_dir + '/pacbio_data')
            os.mkdir(self.output_dir + '/pacbio_data')
            os.link(self.work_dir + "/all.pacbio.fq", self.output_dir + '/pacbio_data/all.pacbio.fq')
            if os.path.exists(self.work_dir + "/" + 'pacbio.rawdata.list'):
                try:
                    if os.path.exists(self.output_dir + "/pacbio_data/" + 'pacbio.rawdata.list'):
                        os.remove(self.output_dir + "/pacbio_data/" + 'pacbio.rawdata.list')
                    os.link(self.work_dir + "/" + 'pacbio.rawdata.list',
                            self.output_dir + "/pacbio_data/" + '/pacbio.rawdata.list')
                except Exception, e:
                    self.set_error('pacbio.rawdata.list的文件link时出错%s', variables=(e), code="22100810")
            self.option('pacbio_dir', self.output_dir + "/pacbio_data")
            self.end()

    def end(self):
        super(FungiGenomeModule, self).end()

    def linkdir(self, dirpath, dirname):
        """
		link一个文件夹下的所有文件到本module的output目录
		:param dirpath: 传入文件夹路径
		:param dirname: 新的文件夹名称
		:return:
		"""
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
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
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))