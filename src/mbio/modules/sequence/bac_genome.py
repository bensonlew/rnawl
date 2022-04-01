#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == zhouxuan

import os
import re,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class BacGenomeModule(Module):
    """
    对数据进行解压统计碱基质量
    """
    def __init__(self, work_id):
        super(BacGenomeModule, self).__init__(work_id)
        options = [
            {"name": "analysis", "type": "string", "default": "uncomplete"}, ##流程分析模式complete，uncomplete
            {"name": "raw_dir", "type": "infile", "format": "bacgenome.raw_dir"},
            {"name": "sample_name", "type": "string"},
            {"name": "sequence_type", "type": "string"},
            {"name": "dir", "type": "outfile", "format": "sequence.fastq_dir"},
            {"name": "pacbio_dir", "type": "outfile", "format": "bacgenome.pacbio_dir1"},  #guanqing.zou 20180905
            {"name": "nanopore_fq", "type": "outfile", "format": "sequence.fastq"}  # gaohao 20210519
        ]
        self.add_option(options)
        self.samplesi = ""
        self.tools = []
        self.h5_tools = []
        self.nano_tools = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option("raw_dir").is_set:
            raise OptionError("必须输入样本文件夹！", code="24000101")
        elif not self.option("analysis"):
            raise OptionError("必须输入流程运行类型！", code="24000102")
        elif not self.option("sample_name"):
            raise OptionError("必须输入样本名！", code="24000103")
        elif not self.option("sequence_type"):
            raise OptionError("必须read的文库类型！", code="24000104")
        else:
            return True

    def run(self):
        super(BacGenomeModule, self).run()
        self.samples = self.get_list()
        if self.option("analysis") in ["uncomplete"]:
            self.run_ungiz()
            if len(self.tools) > 1:
                self.on_rely(self.tools, self.set_output)
            else:
                self.tools[0].on("end",self.set_output)
            for t in self.tools:
                t.run()

        elif self.option("analysis") in ["complete"]:
            self.logger.info(self.nanopore_paths)
            self.logger.info(self.nanopore_paths2)
            self.run_h5_to_fq()
            self.run_bam_to_fq() ## add bam转fq文件
            self.run_nanopore()
            if not re.search(r'PE',self.option("sequence_type")):
                if len(self.h5_tools) == 0 and len(self.nano_tools) >0:
                    if len(self.nano_tools) ==1:
                        self.nano_tools[0].on('end', self.set_output)
                        self.nano_tools[0].run()
                    elif len(self.nano_tools) >1:
                        self.on_rely(self.nano_tools, self.set_output)
                        for t in self.nano_tools:
                            t.run()
                elif len(self.h5_tools) > 0 and len(self.nano_tools) ==0:
                    if len(self.h5_tools) == 1:
                        self.h5_tools[0].on('end', self.set_output)
                        self.h5_tools[0].run()
                    else:
                        self.on_rely(self.h5_tools, self.set_output)
                        for t in self.h5_tools:
                            t.run()

            elif re.search(r'PE', self.option("sequence_type")):
                self.run_ungiz()
                all_tools = []
                self.logger.info(self.tools)
                if len(self.h5_tools) == 0 and len(self.nano_tools) > 0:
                    all_tools = self.tools + self.nano_tools
                elif len(self.h5_tools) > 0 and len(self.nano_tools) == 0:
                    all_tools = self.tools + self.h5_tools
                elif len(self.h5_tools) == 0 and len(self.nano_tools) == 0:
                    all_tools = self.tools
                self.logger.info(all_tools)
                if len(all_tools) > 1:
                    self.on_rely(all_tools, self.set_output)
                    for t in all_tools:
                        t.run()
                else:
                    all_tools[0].on('end', self.set_output)
                    all_tools[0].run()


    def get_list(self):
        sample = {}
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                output2 = self.work_dir + "/" + 'sample_info'
                with open(raw_list_path, "rb") as l,open(output2,"w") as file2:
                    lines = l.readlines()
                    lib_list =[]
                    raw_list =[]
                    for line in lines[1:]:
                        line2 = line.strip().split()
                        if len(line2) == 6 or len(line2) == 7:
                            lib_type = line2[5] + line2[2]
                            lib_list.append(lib_type)
                            raw_list.append(line2[1])
                            file2.write(lib_type + '\t' + line2[2] + '\t' + line2[3] + '\n')
                            if re.search(';', line2[1]):
                                raw_path = line2[1].split(';')
                                for raw in raw_path:
                                    files = raw.split(',')
                                    for i in range(len(files)):
                                        if i == 0:
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
                                    if i == 0:
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
                            self.set_error('list.txt文件格式有误', code="24000101")

        elif self.option("analysis") in ["complete"]:
            #guanqing.zou 20180905>>>
            self.h5_paths = {}
            self.bam_paths = {}
            self.nanopore_paths = {}
            self.nanopore_paths2 = {}
            if not os.path.exists(self.output_dir + "/pacbio_data"):
                os.mkdir(self.output_dir + "/pacbio_data")
            new_path = self.output_dir + "/pacbio_data/all.pacbio.fq"   #guanqing.zou 20180905
            if os.path.exists(new_path):
                os.remove(new_path)
            if not re.search(r'PE', self.option("sequence_type")) and  self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                output = self.work_dir + "/" + 'pacbio.rawdata.fofn'
                with open(raw_list_path, "rb") as l, open(output, "w") as file:
                    lines = l.readlines()
                    #new_path = self.work_dir + "/all.pacbio.fq"   #guanqing.zou 20180905
                    for line in lines[1:]:
                        line2 = line.strip().split()
                        sample_name = line2[0]
                        if sample_name not in self.h5_paths.keys():
                            self.h5_paths[sample_name] = []
                        if sample_name not in self.bam_paths.keys():
                            self.bam_paths[sample_name] = []
                        if len(line2) == 6:
                            if re.search(';', line2[1]):
                                pacs = line2[1].split(';')
                                for pac in pacs:
                                    file_path = pac.split(',')
                                    for k in file_path:
                                        if not k.endswith('bax.h5') and  not k.endswith('bas.h5') :  #zouguanqing 201904
                                            if not k.endswith('.bam'): # qingchen.zhang 20201229 用于限制bam文件不参与后面的质控
                                                sample_path = self.option("raw_dir").prop['path'] + "/" + k
                                                os.system('cat %s >>%s' % (sample_path, new_path))  #guanqing.zou 20180905
                                                file.write(sample_path + '\n')
                                            else:
                                                if k.endswith('.bam'): # qingchen.zhang 20201229 用于限制bam文件不参与后面的质控
                                                    self.bam_paths[sample_name].append(self.option("raw_dir").prop['path'] + "/" + k)
                                        else:
                                            if k.endswith('bax.h5'):
                                                    self.h5_paths[sample_name].append(self.option("raw_dir").prop['path'] + "/" + k)
                            else:
                                file_path = line2[1].split(',')
                                for k in file_path:
                                    if not k.endswith('bax.h5') and  not k.endswith('bas.h5') :
                                        if not k.endswith('.bam'): # qingchen.zhang 20201229 用于限制bam文件不参与后面的质控
                                            sample_path = self.option("raw_dir").prop['path'] + "/" + k
                                            os.system('cat %s >>%s' % (sample_path, new_path))  #guanqing.zou 20180905
                                            file.write(sample_path + '\n')
                                        else:
                                            if k.endswith('.bam'): # qingchen.zhang 20201229 用于限制bam文件不参与后面的质控
                                                self.bam_paths[sample_name].append(self.option("raw_dir").prop['path'] + "/" + k)
                                    else:
                                        if k.endswith('bax.h5'):
                                            self.h5_paths[sample_name].append(self.option("raw_dir").prop['path'] + "/" + k)
                        else:
                            self.set_error('list.txt文件格式有误', code="24000102")

            elif re.search(r'PE', self.option("sequence_type")) and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                output = self.work_dir + "/" + 'pacbio.rawdata.fofn'
                output2 = self.work_dir + "/" + 'sample_info'
                with open(raw_list_path, "rb") as l, open(output, "w") as file,open(output2,'w') as file2:
                    lines = l.readlines()
                    for line in lines[1:]:
                        line2 = line.strip().split()
                        sample_name = line2[0]
                        if sample_name not in self.h5_paths.keys():
                            self.h5_paths[sample_name] = []
                        if sample_name not in self.bam_paths.keys():
                            self.bam_paths[sample_name] = []
                        if len(line2) == 6 or len(line2) == 7:
                            if line2[5] in ["Pacbio", "pacbio", "PACBIO"]:
                                if re.search(';', line2[1]):
                                    pacs = line2[1].split(';')
                                    for pac in pacs:
                                        file_path = pac.split(',')
                                        for k in file_path:
                                            if not k.endswith('bax.h5') and not k.endswith('bas.h5'):
                                                if not k.endswith('.bam'): # qingchen.zhang 20201229 用于限制bam文件不参与后面的质控
                                                    sample_path = self.option("raw_dir").prop['path'] + "/" + k
                                                    os.system('cat %s >>%s' % (sample_path, new_path))  #guanqing.zou 20180905
                                                    file.write(sample_path + '\n')
                                                else:
                                                    if k.endswith('.bam'): # qingchen.zhang 20201229 用于限制bam文件不参与后面的质控
                                                        self.bam_paths[sample_name].append(self.option("raw_dir").prop['path'] + "/" + k)
                                            else:
                                                if k.endswith('bax.h5'):
                                                    self.h5_paths[sample_name].append(self.option("raw_dir").prop['path'] + "/" + k)
                                else:
                                    file_path = line2[1].split(',')
                                    for k in file_path:
                                        if not k.endswith('bax.h5') and  not k.endswith('bas.h5') :
                                            if not k.endswith('.bam'): # qingchen.zhang 20201229 用于限制bam文件不参与后面的质控
                                                sample_path = self.option("raw_dir").prop['path'] + "/" + k
                                                os.system('cat %s >>%s' % (sample_path, new_path))  #guanqing.zou 20180905
                                                file.write(sample_path + '\n')
                                            else:
                                                if k.endswith('.bam'): # qingchen.zhang 20201229 用于限制bam文件不参与后面的质控
                                                    self.bam_paths[sample_name].append(self.option("raw_dir").prop['path'] + "/" + k)
                                        else:
                                            if k.endswith('bax.h5'):
                                                self.h5_paths[sample_name].append(self.option("raw_dir").prop['path'] + "/" + k)
                            elif line2[5] in ["nanopore", "Nanopore", "NANOPORE"]:##gaohao 2021.05.19
                                if line2[1].endswith(".gz") or line2[1].endswith(".tar.gz"):
                                    self.nanopore_paths[sample_name] = self.option("raw_dir").prop['path'] + "/" + line2[1]
                                elif line2[1].endswith(".fq") or line2[1].endswith(".fastq"):
                                    self.nanopore_paths2[sample_name] = self.option("raw_dir").prop['path'] + "/" + line2[1]
                            elif line2[5] in ["PE", "pe", "Pe"]:
                                lib_type = line2[5] + line2[2]
                                file2.write(lib_type + '\t' + line2[2] + '\t' + line2[3] + '\n')
                                if re.search(';', line2[1]):
                                    raw_path = line2[1].split(';')
                                    for raw in raw_path:
                                        files = raw.split(',')
                                        for i in range(len(files)):
                                            if i == 0:
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
                                        if i == 0:
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
                            self.set_error('list.txt文件格式有误', code="24000103")
        return sample

    # guanqing 20190627
    def run_h5_to_fq(self):
        for  s in self.h5_paths.keys():
            if self.h5_paths[s] == []:
                continue
            h5fq_tool = self.add_tool("bacgenome.pacbio_convert2")
            opt = {
                "files" : ' '.join(self.h5_paths[s]),
                "sample_name" : s,
                "bam_result_path": h5fq_tool.output_dir,
                "fq_result_path" : self.output_dir   #h5fq_tool.output_dir
            }
            h5fq_tool.set_options(opt)
            self.h5_tools.append(h5fq_tool)

    def run_bam_to_fq(self):
        """
        如果三代数据只上传bam文件不上传fastq，则需要根据bam文件转fastq文件
        :return:
        """
        for s in self.bam_paths.keys():
            if self.bam_paths[s] == []:
                continue
            h5fq_tool = self.add_tool("bacgenome.pacbio_convert")
            opt = {
                "files" : ' '.join(self.bam_paths[s]),
                "sample_name" : s,
            }
            h5fq_tool.set_options(opt)
            self.h5_tools.append(h5fq_tool)

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
                    self.set_error("序列的方向不对，必须为：l/r", code="24000104")
                gunzip_fastq = self.add_tool('bacgenome.fastq_ungz')
                gunzip_fastq.set_options({
                    "fastq": samples[lib][d],
                    "sample_name": self.option("sample_name"),
                    "direction": direct,
                    "lib_type": lib,
                    "result_path": reslut_path
                })
                self.tools.append(gunzip_fastq)

    def run_nanopore(self):
        """
        解压nanopore的数据
        :return:
        """
        if len(self.nanopore_paths) > 0:
            reslut_path = os.path.join(self.work_dir, "nanopore_dir")
            if not os.path.exists(reslut_path):
                os.mkdir(reslut_path)
            self.logger.info(self.nanopore_paths[self.option("sample_name")])
            gunzip_nanopore = self.add_tool('bacgenome.fastq_ungz')
            gunzip_nanopore.set_options({
                "fastq": self.nanopore_paths[self.option("sample_name")],
                "sample_name": self.option("sample_name"),
                "direction": 'nanopore',
                "result_path": reslut_path
            })
            self.nano_tools.append(gunzip_nanopore)

    def set_output(self):
        ## guanqing 20190627
        if self.h5_tools != []:
            new_path = self.output_dir + "/pacbio_data/all.pacbio.fq"
            output = self.work_dir + "/" + 'pacbio.rawdata.fofn'
            t1 = self.h5_tools[0]
            fq = t1.option('fq_result_path').prop['path']
            if os.path.exists(new_path):
                os.system("cat %s >> %s" %(fq, new_path))
            else:
                os.link(fq,new_path)
            with open(output,'a') as file:
                file.write(new_path+'\n')

        if self.option("analysis") in ["uncomplete"] and self.option("raw_dir").is_set:
            if os.path.exists(self.work_dir + "/ungiz_dir"):
                try:
                    self.linkdir(self.work_dir + "/ungiz_dir", "data")
                except Exception, e:
                    self.set_error('解压的结果linkdir时出错%s', variables=(e), code="24000105")
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
                    self.set_error('sample_info的文件link时出错%s', variables=(e), code="24000106")
            self.end()
        elif self.option("analysis") in ["complete"] and self.option("raw_dir").is_set:
            if not re.search(r'PE',self.option("sequence_type")):
                if os.path.exists(self.work_dir + "/" + 'pacbio.rawdata.fofn') and os.path.getsize(self.work_dir + "/" + 'pacbio.rawdata.fofn'):  #zouguanqing
                    try:
                        if os.path.exists(self.output_dir + "/" + 'pacbio.rawdata.fofn'):
                            os.remove(self.output_dir + "/" + 'pacbio.rawdata.fofn')
                        if os.path.exists(self.work_dir + "/" + 'pacbio.rawdata.fofn'):
                            os.link(self.work_dir + "/" + 'pacbio.rawdata.fofn',
                                self.output_dir + "/" + '/pacbio.rawdata.fofn')
                        if len(os.listdir(self.output_dir + "/pacbio_data")) >=1:
                            self.option('pacbio_dir',self.output_dir + "/pacbio_data")
                    except Exception, e:
                        self.set_error('pacbio.rawdata.list的文件link时出错%s', variables=(e), code="24000107")
            elif re.search(r'PE',self.option("sequence_type")):
                if os.path.exists(self.work_dir + "/ungiz_dir"):
                    try:
                         self.linkdir(self.work_dir + "/ungiz_dir", "data")
                    except Exception, e:
                         self.set_error('解压的结果linkdir时出错%s', variables=(e), code="24000108")
                samples = self.samples
                list2 = os.path.join(self.output_dir, "data/list.txt")
                with open(list2, "wb") as w:
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
                if os.path.getsize(self.work_dir + "/" + 'pacbio.rawdata.fofn') >0:
                    try:
                        #self.option('pacbio_list',self.work_dir + "/" + 'pacbio.rawdata.fofn')  #guanqing.zou 20180905
                        self.option('pacbio_dir',self.output_dir + "/pacbio_data" ) #guanqing.zou 20180905
                        if os.path.exists(self.output_dir + "/" + 'pacbio.rawdata.fofn'):
                            os.remove(self.output_dir + "/" + 'pacbio.rawdata.fofn')
                        if os.path.exists(self.work_dir + "/" + 'pacbio.rawdata.fofn'):
                            os.link(self.work_dir + "/" + 'pacbio.rawdata.fofn',
                                self.output_dir + "/" + 'pacbio.rawdata.fofn')
                        if len(os.listdir(self.output_dir + "/pacbio_data")) >=1:
                            self.option('pacbio_dir',self.output_dir + "/pacbio_data")
                    except Exception, e:
                        self.set_error('pacbio.rawdata.list的文件link时出错%s', variables=(e), code="24000109")
                if len(self.nanopore_paths) >=1:
                    if os.path.exists(self.output_dir + "/" + self.option("sample_name") + ".nanopore.fq"):
                        os.remove(self.output_dir + "/" + self.option("sample_name") + ".nanopore.fq")
                    os.link(self.work_dir+"/nanopore_dir/"+ self.option("sample_name") + ".nanopore.fq",self.output_dir + "/" + self.option("sample_name") + ".nanopore.fq")
                    self.option("nanopore_fq", self.output_dir + "/" + self.option("sample_name") + ".nanopore.fq")
                if len(self.nanopore_paths2) >=1:
                    if os.path.exists(self.output_dir + "/" + self.option("sample_name") + ".nanopore.fq"):
                        os.remove(self.output_dir + "/" + self.option("sample_name") + ".nanopore.fq")
                    os.link(self.nanopore_paths2[self.option("sample_name")],self.output_dir + "/" + self.option("sample_name") + ".nanopore.fq")
                    self.option("nanopore_fq", self.output_dir + "/" + self.option("sample_name") + ".nanopore.fq")
                if os.path.exists(self.work_dir + "/" + 'sample_info'):
                    path = os.path.join(self.output_dir, "data")
                    if not os.path.exists(path):
                        os.mkdir(path)
                    try:
                        if os.path.exists(self.output_dir + "/data/" + 'sample_info'):
                            os.remove(self.output_dir + "/data/" + 'sample_info')
                        os.link(self.work_dir + "/" + 'sample_info', self.output_dir + "/data/" + 'sample_info')
                    except Exception, e:
                        self.set_error('sample_info的文件link时出错%s', variables=(e), code="24000110")
            self.end()

    def end(self):
        super(BacGenomeModule, self).end()

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
