# -*- coding: utf-8 -*-
# __author__ = 'sj'
# last modified by guhaidong 20171109

from __future__ import division
import os
import re
import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool


class SampleRenameAgent(Agent):
    """
    从fastq或者fastq文件夹里提取样本的信息
    """
    def __init__(self, parent):
        super(SampleRenameAgent, self).__init__(parent)
        options = [
            {"name": "file_list", "type": "string", "default": "null"},
            {"name": "info_txt", "type": "infile", "format": "sequence.info_txt"},
            {"name": "file_sample_list", "type": "outfile", "format": "sequence.info_txt"},
        ]
        self.add_option(options)
        self.step.add_steps("sample_rename")
        self.on('start', self.start_sample_rename)
        self.on("end", self.end_sample_rename)

    def start_sample_rename(self):
        self.step.sample_rename.start()
        self.step.update()             

    def end_sample_rename(self):
        self.step.sample_rename.finish()
        self.step.update()

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = 4
        self._memory = "4G" 


class SampleRenameTool(Tool):
    def __init__(self, config):
        super(SampleRenameTool, self).__init__(config)
        old_info_path = self.option("info_txt").prop["path"]
        if self.option("file_list") == "null":
            dict = {}
            with open(old_info_path, "r") as r:
                r.readline()
                for line in r:
                    line = line.strip()
                    lst = line.split("\t")
                    file_name = os.path.basename(lst[0])
                    sample_name = lst[1]
                    key = file_name + "::" + sample_name
                    dict[key] = [sample_name, sample_name]
            self.file_list = dict
        else:
            self.file_list = eval(self.option("file_list"))

    
    def run(self):
        super(SampleRenameTool, self).run()
        old_info_path = self.option("info_txt").prop["path"]
        new_info_path = os.path.join(self.work_dir, "info_tmp.txt")
        self.info_rename(old_info_path, new_info_path)

    
    def info_rename(self, old_info_path, new_info_path):
        self.logger.info("开始进行改名")
        with open(new_info_path, "w") as w:
            with open(old_info_path, "r") as r:
                w.write("#file_path\tsample\twork_dir_path\tseq_num\tbase_num\tmean_length\tmin_length\tmax_length\n")
                r.readline()
                for line in r:
                    line = line.strip()
                    lst = line.split("\t")
                    file_name = os.path.basename(lst[0])
                    sample_name = lst[1]
                    key = file_name + "::" + sample_name
                    if key in self.file_list.keys():
                        new_tool_lst = lst[2].split("/")
                        new_tool_path = self.work_dir + "/" + new_tool_lst[-1]  
                        self.mv(lst[2], new_tool_path, key)
                        w.write(lst[0] + "\t" + self.file_list[key][0] + "\t" +
                                new_tool_path + "\t" + lst[3] + "\t" + lst[4] + "\t" +
                                lst[5] + "\t" + lst[6] + "\t" + lst[7] + "\n")
        self.logger.info("根据file_list改名完成")
        self.create_info(new_info_path)
    
    def mv(self, old_path, new_path, key):
        old_name = self.file_list[key][1]
        new_name = self.file_list[key][0]
        if re.match("^[0-9A-Za-z]+[0-9A-Za-z_]*[0-9A-Za-z]{0,}$", new_name) and not new_name.endswith("_"):  # add by zhujaun 20180428 检测样品名是否含有特殊字符
            pass
        else:
            self.set_error("样本名称中含有除数字、字母和下划线以外的字符，请更改样本名称后再进行工作流分析", code="32706601")
            raise Exception("样本名称中含有除数字、字母和下划线以外的字符，请更改样本名称后再进行工作流分析")
        if new_name == "OUT" or new_name == "IN":
            # raise Exception("样本名称不能为IN与OUT")
            self.set_error("样本名称不能为IN与OUT", code="32706602")
        if not os.path.exists(new_path):
            os.mkdir(new_path)
        output_path = new_path + "/output"
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        output_length_path = output_path + "/length"
        output_fa_path = output_path + "/fa"
        if not os.path.exists(output_length_path):
            os.mkdir(output_length_path)
        if not os.path.exists(output_fa_path):
            os.mkdir(output_fa_path)
        for file in os.listdir(old_path + "/output/length"):
            sample = re.match("(.+)\.length_file", file).group(1)
            if sample == old_name:
                old_file = os.path.join(old_path + "/output/length", file)
                new_file = os.path.join(output_length_path, new_name + ".length_file")
                if os.path.exists(new_file):
                    with open(new_file, "a") as a:
                        with open(old_file, "r") as r:
                            for line in r:
                                a.write(line)
                else:
                    shutil.copy(old_file, new_file)  # 此处不能用os.link()
        for file in os.listdir(old_path + "/output/fa"):
            sample = re.match("(.+)\.fa", file).group(1)
            if sample == old_name:
                old_file = os.path.join(old_path + "/output/fa", file)
                new_file = os.path.join(output_fa_path, new_name + ".fasta")
                if os.path.exists(new_file):
                    with open(new_file, "a") as a:
                        with open(old_file, "r") as r:
                            for line in r:
                                a.write(line)
                else:
                    shutil.copy(old_file, new_file)
                    
    def create_info(self, old_info_path):
        sample_lst = []
        sample_workdir = {}  # sample_workdir[sample] = [path1,path2]  edited by shijin on 20170411
        sample_list_path = {}  # sample_list_path[sample] = sample_workdir[sample][0]
        sample_reads = {}
        sample_bases = {}
        sample_avg = {}
        sample_min = {}
        sample_max = {}
        # if len(eval(self.option("workdir_sample"))) != 1:
        self.logger.info("开始根据样本改名对样本信息进行整理")
        with open(old_info_path, "r") as r:
            with open(self.work_dir + "/info.txt", "w") as w:
                w.write("#file_path\tsample\twork_dir_path\tseq_num\t"
                        "base_num\tmean_length\tmin_length\tmax_length\n")
                r.readline()
                for line in r:
                    line = line.strip()
                    lst = line.split("\t")
                    sample_name = lst[1]
                    if sample_name not in sample_lst:
                        sample_lst.append(sample_name)
                        sample_workdir[sample_name] = [lst[2]]
                        sample_reads[sample_name] = int(lst[3])
                        sample_bases[sample_name] = int(lst[4])
                        sample_avg[sample_name] = int(sample_bases[sample_name]) / int(sample_reads[sample_name])
                        sample_min[sample_name] = int(lst[6])
                        sample_max[sample_name] = int(lst[7])
                    else:
                        sample_workdir[sample_name].append(lst[2])
                        sample_reads[sample_name] += int(lst[3])
                        sample_bases[sample_name] += int(lst[4])
                        sample_avg[sample_name] = int(sample_bases[sample_name]) / int(sample_reads[sample_name])
                        sample_min[sample_name] = min(int(lst[6]), sample_min[sample_name])
                        sample_max[sample_name] = max(int(lst[7]), sample_max[sample_name])
                for sample in sample_lst:
                    sample_list_path[sample] = sample_workdir[sample][0]
                    w.write("-" + "\t" + sample + "\t" + sample_list_path[sample] + "\t"
                            + str(sample_reads[sample]) + "\t" + str(sample_bases[sample]) + "\t"
                            + str(sample_avg[sample]) + "\t" + str(sample_min[sample])
                            + "\t" + str(sample_max[sample]) + "\n")
                if len(sample_lst) == 0:
                    self.set_error("Sample not in file_list, please check sample_check")
        self.logger.info("对样本信息重新整理完成")
        self.cat_fastas(sample_workdir)
        self.logger.info("end")
        self.end()

    def cat_fastas(self, sample_workdir):  # sample_workdir为字典，键值为样本名，值为列表
        self.logger.info("开始对样本名称相同的序列文件进行合并")
        for sample in sample_workdir.keys():
            path_list = sample_workdir[sample]
            path_unrepeat = []
            length_unrepeat = []
            for path in path_list:
                # path = os.path.join(path, sample + ".fasta")
                path = path + "/output/fa/" + sample + ".fasta"
                if path not in path_unrepeat:
                    path_unrepeat.append(path)  # 对列表进行去重，属于同一目录的样本不再进行合并
            if len(path_unrepeat) != 1:
                str = " ".join(path_unrepeat)
                os.system("cat {} > {}/tmp.fa".format(str, self.work_dir))
                os.system("mv {}/tmp.fa {}".format(self.work_dir, path_unrepeat[0]))
                self.logger.info("来自于不同文件的样本的序列文件{}，合并完成".format(sample))
                self.logger.info(path_unrepeat[0])
            for path in path_list:  # 对长度文件进行合并
                path = path + "/output/length/" + sample + ".length_file"
                if path not in length_unrepeat:
                    length_unrepeat.append(path)
            if len(length_unrepeat) != 1:
                str = " ".join(length_unrepeat)
                os.system("cat {} > {}/tmp.length".format(str, self.work_dir))
                os.system("mv {}/tmp.length {}".format(self.work_dir, length_unrepeat[0]))
                self.logger.info("来自于不同文件的样本的长度文件{}，合并完成".format(sample))
        
    def end(self):
        with open(self.work_dir + "/info.txt", "r") as file:
            file.readline()
            for line in file:
                tmp = line.split("\t")
                sample = tmp[1]
                if sample.find(".") != -1 or sample.find(" ") != -1:
                    # raise Exception("样本名称中带.或空格，请更改样本名称为不带.和空格的名称后再进行流程")
                    self.set_error("样本名称中带.或空格，请更改样本名称为不带.和空格的名称后再进行流程", code="32706603")
                if sample == "OUT" or sample == "IN":
                    # raise Exception("样本名称不能为IN与OUT")
                    self.set_error("样本名称不能为IN与OUT", code="32706604")
        self.option("file_sample_list").set_path(self.work_dir + "/info.txt")
        with open(self.option("file_sample_list").prop["path"], "r") as file:
            file.readline()
            try:
                next(file)
            except:
                # raise Exception("样本检测没有找到样本，请重新检查文件的改名信息")
                self.set_error("样本检测没有找到样本，请重新检查文件的改名信息", code="32706605")
        super(SampleRenameTool, self).end()