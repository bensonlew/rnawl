# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from __future__ import division
import os
import re
import json
import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.metaasv.common_function import link_dir, get_sample_check, link_file
from mbio.packages.meta.common_function import get_sg_sample_check
from collections import Counter


class SampleRenameAgent(Agent):
    """
    从fastq或者fastq文件夹里提取样本的信息
    """
    def __init__(self, parent):
        super(SampleRenameAgent, self).__init__(parent)
        options = [
            {"name": "task_id", "type": "string"},##task_id
            {"name": "info_txt", "type": "infile", "format": "sequence.info_txt"}, ##统计的结果表
            {'name': 'denoise_method', 'type': 'string',"default": "DADA2"},  # 降噪方法选择
            {"name": "query_id", "type": "string", "default": ""},## 文件的编号，如果工作流直接运行则必须要传此字段，否则前端无法正常显示样本检测内容
            {"name": "task_type", "type": "string", "default": "qiime2"},## 多样性和qiime2的样本重命名公用模块
            {"name": "file_origin", "type": "string", "default": ""}, ### 文件形式，用于寻找s3的真实文件名称
            {"name": "mapping_file", "type": "string", "default":""}, ## 文件夹形式，用于寻找s3的真实文件名称
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
        self._cpu = 1
        self._memory = "4G" 


class SampleRenameTool(Tool):
    def __init__(self, config):
        super(SampleRenameTool, self).__init__(config)
        self.task_id = self.option("task_id")
        self.info_path = self.option("info_txt").prop["path"]

    def run(self):
        """
        运行
        :return:
        """
        super(SampleRenameTool, self).run()
        self.rename_and_stat()
        self.set_output()
        self.end()

    def rename_and_stat(self):
        """
        功能：根据MongoDB中的记录重命名和统计(质控后的统计表)
        file_sample_dict : 文件名称-origin名称-重命名后的样本名称
        file_sample_dict = {"valid.MJ20180723039.zhaoxiao.338F_806R.format_1539268435.fq": [{"C1": "C2"}]}
        """
        s3_origin = {}
        if self.option("file_origin") != "":
            s3_origin = json.loads(self.option("file_origin"))
        elif self.option("mapping_file") != "":
            mapping_file = self.option("mapping_file")
            self.logger.info(mapping_file)
            with open(mapping_file, 'r') as mm:
                js = mm.read()
                json_dict = json.loads(js)
            if self.option("task_type") in ["meta"]:
                fastq_list = json_dict['in_fastq']
            else:
                fastq_list = json_dict['fastq_file']
            for fastq_dict in fastq_list:
                new_file_name = fastq_dict['alias']
                true_file_name = fastq_dict['file_path']
                if new_file_name not in s3_origin:
                    s3_origin[new_file_name] = true_file_name

        # new_s3_origin = {v.split("/")[-1]:k for k,v in s3_origin.items()}
        # s3_origin = {k:v.split("/")[-1] for k,v in s3_origin.items()}
        new_s3_origin = {v.rstrip(".gz"):k.rstrip(".gz") for k,v in s3_origin.items()} ## 改用文件名称为文件路径名称 qingchen.zhang
        s3_origin = {k.rstrip(".gz"):v.rstrip(".gz") for k,v in s3_origin.items()} ## 改用文件名称为文件路径名称@20201027 qingchen.zhang
        self.logger.info("new_s3_origin: {}".format(new_s3_origin))

        if self.option("task_type") in ["qiime2"]:
            file_sample_dict = get_sample_check(self.task_id, self.option("query_id"))##获取到file_name origin_name和specimen的对应关系
        else:
            file_sample_dict = get_sg_sample_check(self.task_id, self.option("query_id"))## 因为数据库不同，meta专用
        file_sample_dict = {k.rstrip(".gz"): v for k, v in file_sample_dict.items()}
        self.logger.info("info_txt: {}".format(self.option("info_txt").prop['path']))
        self.logger.info("file_sample_dict: {}".format(file_sample_dict))
        last_file = os.path.join(self.work_dir, "info_path.xls")
        temp_all_file = os.path.join(self.work_dir, "info_temp.xls")
        outw = open(temp_all_file, 'w')
        origin_list = []
        origin_sample_dict = {}
        file_list = []
        new_name_origin_name = {}
        file_origin_dict = {}
        file_sample_origin_dict = {} ## 用于记录file+sample的对应关系
        all_origin_list = []
        all_new_sample_name_list = []
        for _file in file_sample_dict:
            file_sample_ = file_sample_dict[_file]
            all_origin_list += [x.keys()[0] for x in file_sample_]
            all_new_sample_name_list += [x.values()[0] for x in file_sample_]
        all_origin_list = list(set(all_origin_list))
        all_new_sample_name_list = list(set(all_new_sample_name_list))
        for new_name in all_new_sample_name_list:
            if re.match("^[0-9A-Za-z]+[0-9A-Za-z_]*[0-9A-Za-z]{0,}$", new_name) and not new_name.endswith("_"):
                pass
            else:
                self.set_error("样本名称中含有除数字、字母和下划线以外的字符，请更改样本名称后再进行工作流分析")
            if new_name == "OUT" or new_name == "IN":
                self.set_error("样本名称不能为IN与OUT")
            if new_name == "NA" or new_name == "GT" or new_name == "OR" or new_name == "na" or new_name == "gt" or new_name == "or":
                self.set_error("样本名称不能为NA(na),GT(gt),OR(or)")
            if new_name.find(".") != -1 or new_name.find(" ") != -1:
                self.set_error("样本名称中带.或空格，请更改样本名称为不带.和空格的名称后再进行流程")

        self.logger.info("all_origin_list: {}".format(all_origin_list))
        self.logger.info("all_new_sample_name_list: {}".format(all_new_sample_name_list))
        with open(self.info_path, 'r') as f, open(last_file,'w') as w:
            w.write("Sample_Name\tRead_Num\tBase_Num\tMean_Length\tMin_Length\tMax_Length\n")
            outw.write("Origin_Name\tSample_Name\tFile_Name\tRead_Num\tBase_Num\tMean_Length\tMin_Length\tMax_Length\n")
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                file_path = line[0]
                sample = line[1]
                read_num = line[3]
                base_num = line[4]
                mean_length = line[5]
                min_length = line[6]
                max_length = line[7]
                file_name_new = os.path.basename(file_path)
                if file_name_new in s3_origin:
                    if file_sample_dict == {}:
                        file_name = s3_origin[file_name_new]
                    else:
                        if "/" in file_sample_dict.keys()[0]:
                            file_name = s3_origin[file_name_new]
                        else:
                            file_name = file_name_new
                else:
                    file_name = file_name_new
                if file_name in file_origin_dict:
                    file_origin_dict[file_name] = sample
                ## 这里需要注意的是文件名称+样本名称为唯一的key
                if file_name not in file_sample_origin_dict:
                    file_sample_origin_dict[file_name] = {sample:read_num + "\t" +base_num + "\t" + mean_length + "\t" + min_length +"\t" + max_length}
                    self.logger.info("origin_sample_dict11111: {}".format(origin_sample_dict))
                    if sample not in origin_list:
                        origin_list.append(sample)
                else:
                    origin_sample_dict = file_sample_origin_dict[file_name]
                    self.logger.info("origin_sample_dict22222: {}".format(origin_sample_dict))
                    if sample not in origin_sample_dict:## 记录下来原始文件+样本对应的统计值和原始样本去冗余后的样本名称
                        origin_sample_dict[sample] = read_num + "\t" +base_num + "\t" + mean_length + "\t" + min_length +"\t" + max_length
                        file_sample_origin_dict[file_name] = origin_sample_dict
                        if sample not in origin_list:
                            origin_list.append(sample)
                if file_name not in file_list and (file_name not in ["#file_path"]):## 记录文件的名称
                    file_list.append(file_name)
                if file_sample_dict == {}:### 没有做样本检测的情况，记录原始名称的对应关系
                    if self.option("denoise_method") in ["Deblur"]:
                        if re.search(r"\_", sample):
                            new_sample = re.sub("_", "", sample)
                        else:
                            new_sample = sample
                            # self.set_error("Deblur软件计算时的样本名称不能含有下划线，否则软件会报错！请重新命名样本名称！")
                    else:
                        new_sample = sample
                    if re.match("^[0-9A-Za-z]+[0-9A-Za-z_]*[0-9A-Za-z]{0,}$", new_sample) and not new_sample.endswith("_"):
                        pass
                    else:
                        self.set_error("样本名称中含有除数字、字母和下划线以外的字符，请更改样本名称后再进行工作流分析")
                    if new_sample == "NA" or new_sample == "GT" or new_sample == "OR" or new_sample == "na" or new_sample == "gt" or new_sample == "or":
                        self.set_error("样本名称不能为NA(na),GT(gt),OR(or)")
                    if sample not in new_name_origin_name:
                        new_name_origin_name[sample] = new_sample
                    if file_name in new_s3_origin:
                        new_file = new_s3_origin[file_name]
                    else:
                        new_file = file_name
                    outw.write("{}\t{}\t{}\t{}\n".format(sample, new_sample, new_file, read_num + "\t" +base_num + "\t" + mean_length + "\t" + min_length +"\t" + max_length))
            self.logger.info("all_origin_list: {}".format(all_origin_list))
            self.logger.info("all_new_sample_name_list: {}".format(all_new_sample_name_list))
            self.logger.info("all_origin_list: {}".format(len(all_origin_list)))
            self.logger.info("all_new_sample_name_list: {}".format(all_new_sample_name_list))
            sample_dict = {}
            if file_sample_dict == {}:### 没有做样本检测的情况，直接通过工作流进行样本检测，将新合并的文件做加和使用
                self.logger.info("直接进行预测得到结果!")
                self.logger.info("file_sample_origin_dict: {}".format(file_sample_origin_dict))
                self.logger.info("origin_list: {}".format(origin_list))
                for ori_name in origin_list:
                    read_num = 0
                    base_num = 0
                    min_length = 0
                    max_length = 0
                    for file_name in file_list:## 对不同文件，相同的样本名称进行合并
                        origin_sample_dict = file_sample_origin_dict[file_name]
                        if ori_name in origin_sample_dict:
                            origin_sample_number_list = origin_sample_dict[ori_name].split("\t")
                            read_num += int(float(read_num)) + int(float(origin_sample_number_list[0]))
                            base_num += int(float(base_num)) + int(float(origin_sample_number_list[1]))
                            mean_length = round(float(float(base_num) / read_num), 6)
                            min_length = min([x for x in [min_length, int(float(origin_sample_number_list[3]))] if x != 0])
                            max_length = max([x for x in [max_length, int(float(origin_sample_number_list[4]))] if x != 0])
                            sample_dict[ori_name] = str(read_num) + "\t" +str(base_num) + "\t" + str(mean_length) + "\t" + str(min_length) +"\t" + str(max_length)
                    if ori_name in sample_dict:
                        stat_number = sample_dict[ori_name]
                        now_name = new_name_origin_name[ori_name]
                        w.write("{}\t{}\n".format(now_name, stat_number))


            elif len(set(all_origin_list)) == len(set(all_new_sample_name_list)):
                ## 做过样本检测，并对样本名称进行重命名，需要对数据进行合并和统计，第一种情况：不同的文件，重命名后的样本名称不会重复
                ## 必须保证原始的样本与改名后的样本存在一一对应的，因为这里如果存在名称相同的情况，可能会造成结果存在问题
                self.logger.info("====================++++++++++++++++++++++")
                for new_sample in all_new_sample_name_list:
                    read_num = 0
                    base_num = 0
                    min_length = 0
                    max_length = 0
                    for file_name2 in file_list:
                        if file_name2 in file_sample_dict.keys():## 进行样本检测的样本与运行工作流的样本是一致的
                            file_sample_list = file_sample_dict[file_name2]
                            file_sample_list.sort()
                            for origin_dict in file_sample_list:
                                origin_name = origin_dict.keys()[0]
                                origin_new_sample = origin_dict.values()[0]
                                origin_sample_dict = file_sample_origin_dict[file_name2]
                                if origin_new_sample in [new_sample]:
                                    rename_name = origin_dict[origin_name]
                                    if self.option("denoise_method") in ["Deblur"]:
                                        if re.search(r"\_", rename_name):
                                            self.set_error("Deblur软件计算时的样本名称不能含有下划线，否则软件会报错！请重新命名样本名称！")
                                    origin_sample_number_list = origin_sample_dict[origin_name].split("\t")
                                    read_num += int(float(origin_sample_number_list[0]))
                                    base_num += int(float(origin_sample_number_list[1]))
                                    mean_length = round(float(float(base_num) / read_num), 6)
                                    min_length = min([x for x in [min_length, int(float(origin_sample_number_list[3]))] if x != 0])
                                    max_length = max([x for x in [max_length, int(float(origin_sample_number_list[4]))] if x != 0])
                                    sample_dict[new_sample] = str(read_num) + "\t" +str(base_num) + "\t" + str(mean_length) + "\t" + str(min_length) +"\t" + str(max_length)
                                # else:
                                #     sample_dict[new_sample] = origin_sample_dict[origin_name]
                # for new_sample in all_new_sample_name_list:
                #     w.write("{}\t{}\n".format(new_sample, sample_dict[new_sample]))
                record_newname_list = []
                for specimen in all_origin_list:## 获取选取标记的样本
                    for file_name3 in file_list:
                        origin_sample_dict = file_sample_origin_dict[file_name3]
                        file_sample_list = file_sample_dict[file_name3]
                        if len(file_sample_list) != 0:
                            new_sample_name_aa = [dict for dict in file_sample_list if dict.keys()[0] == specimen]
                            if len(new_sample_name_aa) != 0:
                                new_sample_name = new_sample_name_aa[0].values()[0]
                                if specimen in origin_sample_dict:
                                    if file_name3 in new_s3_origin:
                                        new_file = new_s3_origin[file_name3]
                                    else:
                                        new_file = file_name3
                                    outw.write("{}\t{}\t{}\t{}\n".format(specimen, new_sample_name, new_file, origin_sample_dict[specimen]))
                                    if specimen not in record_newname_list:## add by qingchen.zhang 20201209原因是为防止不同文件中的样本合并出现重复的情况
                                        record_newname_list.append(specimen)
                                        w.write("{}\t{}\n".format(new_sample_name, sample_dict[new_sample_name]))

            elif len(set(all_origin_list)) != len(set(all_new_sample_name_list)):
                ## 做样本检测，并对样本名称进行重命名，需要对数据进行合并和统计，第二种情况：不同的文件，重命名后的样本名称会重复
                repeated_list = []
                self.logger.info("11111111111111111111111111")
                all_origin_list = []
                all_new_sample_name_list = []
                for _file in file_sample_dict:
                    file_sample_ = file_sample_dict[_file]
                    all_origin_list += [x.keys()[0] for x in file_sample_]
                    all_new_sample_name_list += [x.values()[0] for x in file_sample_]
                all_new_sample_name_dict = Counter(all_new_sample_name_list)
                for name in all_new_sample_name_list:
                    value = all_new_sample_name_dict[name]
                    if int(value) > 1:
                        if name not in repeated_list:
                            repeated_list.append(name)##得到重复样本
                new_origin_sample_dict = {}
                self.logger.info("repeated_list: {}".format(repeated_list))
                if len(repeated_list) == 0:##这里对原始样本有重名的情况进行判断
                    read_num = 0
                    base_num = 0
                    min_length = 0
                    max_length = 0
                    for file_name in file_list:##遍历文件路径
                        origin_sample_dict = file_sample_origin_dict[file_name] ##从工作流得到的检测结果
                        self.logger.info("file_name: {}".format(file_name))
                        all_origin_list = set(all_origin_list)
                        for origin_name in all_origin_list:## 遍历原始样本名称
                            file_sample_list = file_sample_dict[file_name]## 页面检测得到的样本合并结果
                            self.logger.info("origin_name: {}".format(origin_name))
                            if len(file_sample_list) != 0:
                                new_sample_name_li = [dict for dict in file_sample_list if dict.keys()[0] == origin_name]
                                self.logger.info("new_sample_name_li: {}".format(new_sample_name_li))
                                if origin_name in origin_sample_dict:
                                    all_stat = origin_sample_dict[origin_name]
                                # else:
                                #     all_stat = str(read_num) + "\t" +str(base_num) + "\t" + str(0) + "\t" + str(min_length) +"\t" + str(max_length)
                                    if len(new_sample_name_li) != 0:
                                        new_sample_name = new_sample_name_li[0].values()[0]
                                        self.logger.info("new_sample_name: {}".format(new_sample_name))
                                        if new_sample_name in new_origin_sample_dict:## 判断新的样本名称是否存在合并的情况
                                            self.logger.info("开始合并：{}".format(new_sample_name))
                                            all_stat_new = new_origin_sample_dict[new_sample_name]
                                            all_stat_list = all_stat_new.split("\t")
                                            read_num += int(float(all_stat_list[0]))
                                            base_num += int(float(all_stat_list[1]))
                                            mean_length = float(base_num / read_num)
                                            min_length = min([x for x in [min_length, int(float(all_stat_list[3]))] if x != 0])
                                            max_length = max([x for x in [max_length, int(float(all_stat_list[4]))] if x != 0])
                                            new_all_stat = str(read_num) + "\t" +str(base_num) + "\t" + str(mean_length) + "\t" + str(min_length) +"\t" + str(max_length)
                                            new_origin_sample_dict[new_sample_name] = new_all_stat
                                        else:
                                            self.logger.info("新样本名称：{}, 原始样本名称：{}".format(new_sample_name, origin_name))
                                            new_origin_sample_dict[new_sample_name] = all_stat

                for repeated_name in repeated_list:## 对重复的新命名的样本进行合并统计
                    read_num = 0
                    base_num = 0
                    min_length = 0
                    max_length = 0
                    for file_name in file_list:
                        origin_sample_dict = file_sample_origin_dict[file_name]
                        for origin_name in all_origin_list:
                            file_sample_list = file_sample_dict[file_name]
                            self.logger.info("origin_name: {}".format(origin_name))
                            if len(file_sample_list) != 0:
                                new_sample_name_li = [dict for dict in file_sample_list if dict.keys()[0] == origin_name]
                                self.logger.info("new_sample_name_li: {}".format(new_sample_name_li))
                                if origin_name in origin_sample_dict:
                                    all_stat = origin_sample_dict[origin_name]
                                else:
                                    all_stat = str(read_num) + "\t" +str(base_num) + "\t" + str(0) + "\t" + str(min_length) +"\t" + str(max_length)
                                if len(new_sample_name_li) != 0:
                                    new_sample_name = new_sample_name_li[0].values()[0]
                                    self.logger.info("new_sample_name: {}".format(new_sample_name))
                                    self.logger.info("repeated_name: {}".format(repeated_name))
                                    if new_sample_name in [repeated_name]:
                                        all_stat_list = all_stat.split("\t")
                                        read_num += int(float(all_stat_list[0]))
                                        base_num += int(float(all_stat_list[1]))
                                        mean_length = float(base_num / read_num)
                                        min_length = min([x for x in [min_length, int(float(all_stat_list[3]))] if x != 0])
                                        max_length = max([x for x in [max_length, int(float(all_stat_list[4]))] if x != 0])
                                        new_all_stat = str(read_num) + "\t" +str(base_num) + "\t" + str(mean_length) + "\t" + str(min_length) +"\t" + str(max_length)
                                        new_origin_sample_dict[new_sample_name] = new_all_stat
                                    else:
                                        # if new_sample_name in new_origin_sample_dict:
                                        #     all_stat_list = all_stat.split("\t")
                                        #     read_num += int(float(all_stat_list[0]))
                                        #     base_num += int(float(all_stat_list[1]))
                                        #     mean_length = float(base_num / read_num)
                                        #     min_length = min([x for x in [min_length, int(float(all_stat_list[3]))] if x != 0])
                                        #     max_length = max([x for x in [max_length, int(float(all_stat_list[4]))] if x != 0])
                                        #     new_all_stat = str(read_num) + "\t" +str(base_num) + "\t" + str(mean_length) + "\t" + str(min_length) +"\t" + str(max_length)
                                        #     new_origin_sample_dict[new_sample_name] = new_all_stat
                                        # else:
                                        if new_sample_name not in repeated_list:
                                            new_origin_sample_dict[new_sample_name] = all_stat

                for new_sample in set(all_new_sample_name_list):
                    if self.option("denoise_method") in ["Deblur"]:
                        if re.search(r"\_", new_sample):
                            self.set_error("Deblur软件计算时的样本名称不能含有下划线，否则软件会报错！请重新命名样本名称！")

                    w.write("{}\t{}\n".format(new_sample, new_origin_sample_dict[new_sample]))

                for specimen in origin_list:## 获取选取标记的样本
                    for file_name3 in file_list:
                        origin_sample_dict = file_sample_origin_dict[file_name3]
                        file_sample_list = file_sample_dict[file_name3]
                        if len(file_sample_list) != 0:
                            new_sample_name_aa = [dict for dict in file_sample_list if dict.keys()[0] == specimen]
                            if len(new_sample_name_aa) != 0:
                                new_sample_name = new_sample_name_aa[0].values()[0]
                                if specimen in origin_sample_dict:
                                    if file_name3 in new_s3_origin:
                                        new_file = new_s3_origin[file_name3]
                                    else:
                                        new_file = file_name3
                                    outw.write("{}\t{}\t{}\t{}\n".format(specimen, new_sample_name, new_file, origin_sample_dict[specimen]))

            outw.close()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        link_file(os.path.join(self.work_dir, "info_path.xls"), os.path.join(self.output_dir, "info_path.xls"))
        link_file(os.path.join(self.work_dir, "info_temp.xls"), os.path.join(self.output_dir, "info_temp.xls"))
        self.logger.info("设置结果目录完成！")