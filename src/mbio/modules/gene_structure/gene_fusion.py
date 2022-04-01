#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile
import re
import shutil
import math

class GeneFusionModule(Module):
    """
    
    version 1.0
    author: sj
    last_modify: 2017.01.17 by zx
    """
    def __init__(self, work_id):
        super(GeneFusionModule, self).__init__(work_id)
        options = [
            {"name": "sickle_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 质量剪切输出结果文件夹(包括左右段)  /mnt/ilustre/users/sanger-dev/workspace/20161110/RefrnaNoAssemble_sj_ref_rna/QualityControl/output/sickle_dir
            {"name": "seqprep_dir", "type": "infile", "format": "sequence.fastq_dir"},  # PE的去接头输出结果文件  /mnt/ilustre/users/sanger-dev/workspace/20161110/RefrnaNoAssemble_sj_ref_rna/QualityControl/output/seqprep_dir
            {"name": "reads_number", "type": "int", "default": 3},  # 支持融合的reads数目的最小值 3 5 10
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # fq文件夹 /mnt/ilustre/users/sanger-dev/sg-users/shijin/PE_test/paired_end_fq_dir
            {"name": "fq_type", "type": "string", "default": "PE"},  # PE或SE
            {"name": "fusion_result_dir", "type": "outfile", "format": "ref_rna.gene_fusion.fusion_result_dir"} #modify by zx 2016.11.21
            #module 的输出文件夹没有设置内容 zx
        ]
        self.add_option(options)
        self.samples_number = 0
        self.sample_names = []
        self.file_sample = {}
        self.gz_path = {}
        
    def check_options(self):
        """
        检查参数
        """
        if self.option("fq_type") != "PE":
            raise OptionError("非双端测序不能进行基因融合分析")
        if not self.option("sickle_dir").is_set:
            raise OptionError("未检测到sickle文件夹")
        else:
            list_path = os.path.join(self.option("sickle_dir").prop["path"], "list.txt")
            if not os.path.exists(list_path):
                raise OptionError("未检测到list文件")
            list = FileSampleFile()
            list.set_path(list_path)
            list.get_info()
            self.logger.info("samples number: " + str(list.prop["sample_number"]))
            self.logger.info("sample names :" + str(list.prop["sample_names"]))
            self.sample_number = list.prop["sample_number"]
            self.sample_names = list.prop["sample_names"]
            self.file_sample = list.get_list()
            # list.prop["file_sample"] 
            self.logger.info(self.file_sample)
        return True
        
    def get_sample_length(self,sample):
        length_l = 0
        length_r = 0
        sample_l_path = self.file_sample[sample]["l"]
        sample_l_path = os.path.join(self.option("sickle_dir").prop["path"],sample_l_path)
        with open(sample_l_path,"r") as r:
            for line in r:
                m = re.match("^@",line)
                if m:
                    line2 = r.next()
                    length_l = len(line2)
                    break
        sample_r_path = self.file_sample[sample]["r"]
        sample_r_path = os.path.join(self.option("sickle_dir").prop["path"],sample_r_path)
        with open(sample_r_path,"r") as r:
            for line in r:
                m = re.match("^@",line)
                if m:
                    line2 = r.next()
                    length_r = len(line2)
                    break
        lent = str(length_l) + "/" + str(length_r)
        return lent
    
    # def create_list_file(self):
    #     list_path = os.path.join(self.work_dir, "sample_list")
    #     self.logger.info(list_path)
    #     with open(list_path, "w") as w:
    #         for i in range(self.sample_number):
    #             w.write(self.sample_names[i] + "\t" + "Lib-" + chr(97 + i) + "\t" + "Run-" + chr(97+i) + "\t" + self.get_sample_length(self.sample_names[i]) + "\n")
    #     self.logger.info("done!")
    #     return True

    def create_list_file(self):
        all_list_file_path = os.path.join(self.work_dir, "sample_list")
        for i in self.sample_names:
            list_path = os.path.join(self.work_dir, i)
            self.logger.info(list_path)
            with open(list_path, "w") as w:
                w.write(i + "\t" + "Lib-a" + "\t" + "Run-a" + "\t" + self.get_sample_length(
                    i) + "\n")
            with open(all_list_file_path, "a") as f:
                f.write(i + "\t" + "Lib-a" + "\t" + "Run-a" + "\t" + self.get_sample_length(
                    i) + "\n")
        self.logger.info("done!")
        return True
        
    def get_gz_path(self):
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        list = FileSampleFile()
        list.set_path(list_path)
        print(list_path)
        self.gz_path = list.get_list()
        print(list.get_list())
        for sample in self.sample_names:
            new_path = self.gz_path[sample]["l"]
            new_path = new_path + "_seqprep_l.gz"    #  edited by sj on 20161125
            seqprep_path = self.option("seqprep_dir").prop["path"]
            if os.path.exists(seqprep_path + "/" + new_path):
                self.gz_path[sample]["l"] = new_path
            else:
                raise OptionError("不能在seqprep文件夹中找到文件{}".format(new_path))
            new_path = self.gz_path[sample]["r"]
            new_path = new_path + "_seqprep_r.gz"
            seqprep_path = self.option("seqprep_dir").prop["path"]
            if os.path.exists(seqprep_path + "/" + new_path):
                self.gz_path[sample]["r"] = new_path
            else:
                raise OptionError("不能在seqprep文件夹中找到文件{}".format(new_path))
        self.logger.info(self.gz_path)
        return True
    
    def create_dir(self,list_path,mode = "copy"):
        first_dir_path = os.path.join(self.work_dir,"WHOLE_SEQ_DIR")
        if os.path.exists(first_dir_path):
            shutil.rmtree(first_dir_path)
        with open(list_path,"r") as r:
            for line in r:
                sample = line.strip().split("\t")[0]
                dir = line.strip().split("\t")[1]
                name = line.strip().split("\t")[2]
                first_dir_path = os.path.join(self.work_dir,"WHOLE_SEQ_DIR")
                if not os.path.exists(first_dir_path):
                    os.mkdir(first_dir_path)
                second_dir_path = os.path.join(first_dir_path,sample)
                if not os.path.exists(second_dir_path):
                    os.mkdir(second_dir_path)
                third_dir_path = os.path.join(second_dir_path,dir)
                if not os.path.exists(third_dir_path):
                    os.mkdir(third_dir_path)
                if mode == "link":
                    old_path = os.path.join(self.option("seqprep_dir").prop["path"],self.gz_path[sample]["l"])
                    new_path = os.path.join(third_dir_path, name+"_1.fq.gz")
                    os.link(old_path,new_path)
                    old_path = os.path.join(self.option("seqprep_dir").prop["path"],self.gz_path[sample]["r"])
                    new_path = os.path.join(third_dir_path, name+"_2.fq.gz")
                    os.link(old_path,new_path)
                if mode == "copy":  # edited by sj on 20161128 
                    old_path = os.path.join(self.option("seqprep_dir").prop["path"],self.gz_path[sample]["l"])
                    new_path = os.path.join(third_dir_path, name+"_1.fq.gz")
                    shutil.copy(old_path,new_path)
                    old_path = os.path.join(self.option("seqprep_dir").prop["path"],self.gz_path[sample]["r"])
                    new_path = os.path.join(third_dir_path, name+"_2.fq.gz")
                    shutil.copy(old_path,new_path)
        self.logger.info("create dir done")
        return True
    
    # def tool_run(self):
    #     opts = {
    #         "sample_data": self.work_dir + "/" + "WHOLE_SEQ_DIR",
    #         "sample_list": self.work_dir + "/" + "sample_list",
    #         "reads_number": self.option("reads_number")
    #         #"sample_name": 'ERR1464149'  # 修改 zx
    #     }
    #     tool = self.add_tool("ref_rna.gene_fusion.soapfuse")
    #     tool.set_options(opts)
    #     tool.on("end", self.get_message)  # edited by sj
    #     print "111"
    #     tool.run()

    def tool_run(self):
        times = int(math.ceil(self.sample_number/5))
        self.logger.info(times)
        self.tools = []
        n = 0
        for i in self.sample_names:
            tool = self.add_tool("ref_rna.gene_fusion.soapfuse")
            self.step.add_steps('gene_fusion_{}'.format(i))
            sample_data = os.path.join(self.work_dir, "WHOLE_SEQ_DIR")
            sample_list = os.path.join(self.work_dir, i)
            tool.set_options = (
                {
                    "sample_data": sample_data,
                    "sample_list": sample_list,
                    "reads_number": self.option("reads_number"),
                    "sample_name": i,
                }
            )
            a = type(i)
            print(i)
            print(a)
            step = getattr(self.step, 'gene_fusion_{}'.format(i))
            step.start()
            tool.on('end', self.finish_update, 'gene_fusion_{}'.format(i))
            self.tools.append(tool)
            n += 1
        if len(self.tools) == 1:
            self.tools[0].on('end', self.get_message)
        elif 1 < len(self.tools) <= 5:
            self.on_rely(self.tools, self.get_message)
        else:
            if times == 1:
                a = self.tools[0:5]
                b = self.tools[6:-1]
                self.on_rely(a, b)
                if len(self.tools) == 6:
                    self.tools[5].on('end', self.get_message)
                else:
                    self.on_rely(self.tools[6:-1], self.get_message)
            else:
                for i in range(0, times):
                    self.on_rely(self.tools[i*5:i*5+4], self.tools[i*5+5:i*5+9])
                    if i*5+9+5 >= len(self.tools):
                        number = i
                        break
                self.on_rely(self.tools[number*5+5:number*5+9], self.tools[number*5+10:-1])
                if len(self.tools) == number*5+9+2:
                    self.tools[-1].on('end', self.get_message)
                else:
                    self.on_rely(self.tools[number*5+10:-1], self.get_message)
        for tool in self.tools[0:5]:
            tool.run()
            # tool.on("end", self.get_message)  # edited by sj
        # print "111"
        # tool.run()


    def finish_update(self):  #  修改 zx
        self.step.tool.finish()
        self.step.update()

    def get_message(self): #获得可视化数据的过程并将结果复制到output文件夹下 zx
        list = FileSampleFile()
        for i in range(self.sample_number):
            # path_old2 = self.work_dir + "/Soapfuse/output/OUT/" + self.sample_names[i] + "/" + self.sample_names[i] + ".final.Fusion.specific.for.trans"
            # path_new2 = self.output + "/" + self.sample_names[i] + ".final.Fusion.specific.for.trans"
            # shutil.copy(path_old2, path_new2)
            path_old = self.work_dir + "/Soapfuse/output/finalresult/" + self.sample_names[i] + "/" + self.sample_names[i] +".final.Fusion.specific.for.genes"
            # path_new = self.output + "/" + self.sample_names[i] +".final.Fusion.specific.for.genes"
            # shutil.copy(path_old, path_new)
            a = open(path_old, 'r')
            b = open("/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/SOAPfuse-v1.26/get_len_number.txt",'r')
            content = a.readlines()
            ref_number = b.readlines()
            for f in content:
                arr = f.strip().split("\t")
                if (arr[0] == "up_gene") == False:
                    up_name = arr[0]
                    up_chr = arr[1]
                    for c in ref_number:
                        brr = c.strip().split("\t")
                        if brr[0] == up_chr:
                            up_number = int(arr[3]) * int(brr[2]) / int(brr[1]) #去括号
                    down_name = arr[5]
                    down_chr = arr[6]
                    for c in ref_number:
                        brr = c.strip().split("\t")
                        if brr[0] == down_chr:
                            down_number = int(arr[8]) * int(brr[2]) / int(brr[1]) #去括号
                    linecontent = up_name + "\t" + up_chr + "\t" + str(up_number) + "\t" + down_name + "\t" + down_chr + "\t" + str(down_number) + "\n"
                    with open(self.work_dir + "/Soapfuse/output/finalresult" + self.sample_names[i] + "/"+ self.sample_names[i] + ".message.txt" , "a") as w:
                        w.write(linecontent)
        self.logger.info("get message done")
        return True

    def set_output(self): #zx
        self.logger("设置module的结果目录")
        try:
            old = self.work_dir + "/Soapfuse/output/finalresult"
            new = self.output_dir + "/fusion_result"
            shutil.copytree(old , new)
            self.option("fusion_result_dir").set_path(self.output_dir + "/fusion_result")
            self.logger.info("设置基因融合module分析结果目录成功")
        except Exception as e:
            self.logger.info("设置基因融合module分析结果目录失败{}".format(e))
            self.set_error("设置基因融合module分析结果目录失败{}".format(e))

    
    def run(self):
        super(GeneFusionModule,self).run()
        self.create_list_file()
        self.get_gz_path()
        list_path = os.path.join(self.work_dir,"sample_list")
        self.create_dir(list_path)
        self.tool_run()

    def end(self):
        super(GeneFusionModule, self).end()
