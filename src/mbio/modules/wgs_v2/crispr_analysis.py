# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modify:20190321

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import json


class CrisprAnalysisModule(Module):
    """
     crispr接口module，首先使用riserach2.x软件对数据进行处理，随后使用resso软件对不同的样本进行处理。
     当有wt文件的时候，mut文件需要先与wt文件对比后进行统计，如果没有wt文件的话，不匹配直接进行统计。

     本接口只应用于二倍体，如果不是二倍体，需要修改导表

    """
    def __init__(self, work_id):
        super(CrisprAnalysisModule, self).__init__(work_id)
        options = [
            {"name": "project_sn", "type": 'string'},
            {"name": "main_id", "type": 'string'},
            {"name": "task_id", "type": 'string'},
            {"name": "update_info", "type": 'string'},
            {"name": "ref_fa", "type": 'string'},
            {"name": "sg_rna", "type": 'string'},  # 这里需要自己写入文件
            {"name": "bam_list", "type": 'string'},
            {"name": "region", "type": 'string'},  # 格式为chr2:start-end  chr1:112194-112213
            {"name": "gene_name", "type": 'string'},  # 用于最后导表的时候基因名字显示
            {"name": "bam_path", "type": 'string'},  # 用于最后导表的时候基因名字显示
            {"name": "wt", "type": 'string'},  # 作用同bam_list，名字加上bam_path前缀，最后写入文档。
        ]
        self.add_option(options)
        self.vcftools_filter_tools = []
        self.ld_decay_tools = []
        self.group_sample = {}

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("请输入ref.fa参数")
        if not self.option("sg_rna"):
            raise OptionError("请输入sg_rna参数")
        if not self.option("bam_list"):
            raise OptionError("请输入bam_list文件")
        if not self.option("region"):
            raise OptionError("请输入region参数")
        if not self.option("gene_name"):
            raise OptionError("请输入gene_name参数")
        if not self.option("bam_path"):
            raise OptionError("请输入bam_path参数")
        # if not self.option("wt"):
        #     raise OptionError("请输入wt参数")
        return True

    def get_files(self):
        """
        本函数的目的在于整理传入的参数，形成文件，包括sgrna_fa文件和bam_list文件
        :return:
        """
        with open(self.output_dir + "/sgRNA.fa", "w") as w1, open(self.output_dir + "/bam_list", "w") as w2:
            w1.write(">sgrna\n" + self.option("sg_rna"))
            for i in json.loads(self.option("bam_list")):
            # for i in self.option("bam_list").strip().split(","):
                w2.write(i + "\t" + os.path.join(self.option("bam_path"), i + ".sort.bam") + "\n")
            w2.write("WT" + "\t" + os.path.join(self.option("bam_path"), self.option("wt") + ".sort.bam") + "\n")

    def run_risearch(self):
        crispr_risearch = self.add_tool("wgs_v2.crispr_risearch")
        crispr_risearch.set_options({
            "ref_fa": self.option("ref_fa"),
            "sgrna_fa": self.output_dir + "/sgRNA.fa"
        })
        crispr_risearch.on('end', self.set_output, 'crispr_risearch_dir')
        crispr_risearch.on('end', self.run_resso)
        crispr_risearch.run()

    def run_resso(self):
        resso_tools = []
        with open(self.output_dir + "/bam_list") as f:
            lines = f.readlines()
            for line in lines:
                sample_name = line.strip().split("\t")[0]
                sample_path = line.strip().split("\t")[1]
                crispr_resso = self.add_tool("wgs_v2.crispr_resso")
                crispr_resso.set_options({
                    "bam_file": sample_path,
                    "name": sample_name,
                    "ref_fa": self.option("ref_fa"),
                    "sgrna_region": self.work_dir + "/CrisprRisearch/output/sgRNA.region",
                })
                resso_tools.append(crispr_resso)
            for j in range(len(resso_tools)):
                if resso_tools[j].option("name") == "WT":
                    resso_tools[j].on("end", self.set_output, 'crispr_resso_dir_wt')
                else:
                    resso_tools[j].on("end", self.set_output, 'crispr_resso_dir_mut')
            if resso_tools:
                if len(resso_tools) > 1:
                    self.on_rely(resso_tools, self.run_stat)
                elif len(resso_tools) == 1:
                    resso_tools[0].on('end', self.run_stat)
            else:
                raise Exception("resso_tools列表为空！")
            for tool in resso_tools:
                gevent.sleep(1)
                tool.run()

    def run_stat(self):
        stat_tools = []
        mut_dir_list = os.listdir(self.output_dir + "/mut")
        if os.path.exists(self.output_dir + "/wt/"):
            wt_dir = self.output_dir + "/wt/" + os.listdir(self.output_dir + "/wt")[0]
        else:
            os.makedirs(self.output_dir + "/wt/")
            wt_dir = self.output_dir + "/wt/"
        for i in mut_dir_list:
            name = i.strip().split("_on_")[1]
            crispr_stat = self.add_tool("wgs_v2.crispr_calc")
            crispr_stat.set_options({
                "wt_dir": wt_dir,
                "mut_dir": self.output_dir + "/mut/" + i,
                "sgrna_region": self.work_dir + "/CrisprRisearch/output/sgRNA.region",
                "name": name
            })
            """"
            因为各个样本分开跑，所以这里就需要添加一个name字段以区分不同样本的结果"""
            stat_tools.append(crispr_stat)
        for j in range(len(stat_tools)):
            stat_tools[j].on("end", self.set_output, 'crispr_stat_dir')
        if stat_tools:
            if len(stat_tools) > 1:
                self.on_rely(stat_tools, self.run_results)
            elif len(stat_tools) == 1:
                stat_tools[0].on('end', self.run_results)
        else:
            raise Exception("stat_tools列表为空！")
        for tool in stat_tools:
            gevent.sleep(1)
            tool.run()

    def run_results(self):
        R_num = []  # 用于存储R_num的所有值，包含过滤之前的
        R_sim = []
        R_num_new = [] # 用于存储新的R值，来源为CRISPResso文件夹下的。
        chr = self.option("region").strip().split(":")[0]
        span = range(int(self.option("region").strip().split(":")[1].strip().split("-")[0]),
                     int(self.option("region").strip().split(":")[1].strip().split("-")[1])
                     + 1)
        with open(self.work_dir + "/CrisprRisearch/output/sgRNA.region") as f:
            lines = f.readlines()
            for line in lines:
                chr_f = line.strip().split("\t")[0]
                if chr_f == chr:
                    if int(line.strip().split("\t")[1]) in span and int(line.strip().split("\t")[2]) in span:
                        R_num.append(line.strip().split("\t")[3])
                        R_sim.append(line.strip().split("\t")[3] + "_" + line.strip().split("\t")[0] + ":" +
                                     line.strip().split("\t")[1] + "-" + line.strip().split("\t")[2])
        files = os.listdir(self.output_dir + "/mut/" + os.listdir(self.output_dir + "/mut")[0])
        """
        并不是每个样本产生Resso文件都一样
        """
        for file_name in files:
            if re.match("CRISPResso_on_.*", file_name):
                if file_name.strip().split("_")[-1] in R_num:
                    R_num_new.append( file_name.strip().split("_")[-1])
        self.logger.info("这是files文件夹")
        self.logger.info(R_num_new)
        if len(R_num_new) == 1:   # 当大于1或者等于0的时候没有比较的意义了。
            path = self.output_dir +"/crispr_stat_dir/"
            files_list = os.listdir(path)
            with open(self.output_dir + "/off_target.txt", "w") as w, open(self.output_dir + "/on_target.txt", "w") as w1, \
                    open(self.output_dir + "/variant_type.txt", "w") as r:
                w.write("sample_name\tmodified_num\ttotal_num\n")
                w1.write("Pos\tsample_name\twt_sequence\tmut_sequence\tvariant_type\n")
                r.write("\tinsertion\tdeletion\tsubstitution\tcombination\n")
                for file in files_list:
                    with open(path + "/" + file) as m:
                        list_pos = []
                        list_modified = []
                        sequence_wt = []
                        sequence_mut = []
                        sequence_type = []
                        sequence_pos = []
                        sequence_num = []  # 用于统计插入或者确实碱基
                        lines = m.readlines()
                        insertion = 0
                        deletion = 0
                        substitutition = 0
                        combination = 0
                        for line in lines[1:]:
                            item = line.strip().split("\t")
                            try:
                                if not item[6] == "":
                                    if item[0] not in R_num and item[0].strip().split("new")[1] not in R_num:
                                        list_pos.append(item[0])
                                        if item[5] == "modified":
                                            list_modified.append(item[0])
                                    else:
                                        sequence_wt.append(item[2])
                                        sequence_mut.append(item[3])
                                        # sequence_pos.append(item[1])
                                        sequence_pos.append(R_sim[0].strip().split("_")[1])
                                        if item[5] == "unmodified":
                                            sequence_type.append("wt")
                                            sequence_num.append("")
                                        else:
                                            if re.match(".*-.*", item[2]) and not re.match(".*-.*", item[3]):
                                                sequence_type.append("i")
                                                dict1 = {}
                                                for i in item[2]:
                                                    dict1[i] = str(item[2]).count(i)
                                                sequence_num.append(dict1["-"])
                                                insertion += 1
                                            elif re.match(".*-.*", item[3]) and not re.match(".*-.*", item[2]):
                                                sequence_type.append("d")
                                                dict2 = {}
                                                for j in item[3]:
                                                    dict2[j] = str(item[3]).count(j)
                                                sequence_num.append(dict2["-"])
                                                deletion += 1
                                            elif re.match(".*-.*", item[2]) and re.match(".*-.*", item[3]):
                                                sequence_type.append("c")
                                                combination += 1
                                                dict3 = {}
                                                for k in item[2]:
                                                    dict3[k] = str(item[2]).count(k)
                                                dict4 = {}
                                                for n in item[3]:
                                                    dict4[n] = str(item[3]).count(n)
                                                sequence_num.append(dict3["-"] + dict4["-"])
                                            elif not re.match(".*-.*", item[2]) and not re.match(".*-.*", item[3]):
                                                sequence_type.append("s")
                                                sequence_num.append(0)
                                                substitutition += 1
                            except:
                                continue
                        w.write(file + "\t" + str(len(list_modified)) + "\t" + str(len(list_pos)) + "\n")
                        r.write(file + "\t" + str(insertion) + "\t" + str(deletion) + "\t" + str(substitutition) + "\t"
                                + str(combination) + "\n")
                        # print sequence_wt
                        # print sequence_mut
                        # print sequence_type
                        # print sequence_pos
                        print str(file.strip().split(".")[0]) + "_____________________"
                        print "**************************"
                        print sequence_wt
                        print sequence_mut
                        print sequence_type
                        print sequence_pos
                        print sequence_num
                        print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
                        if len(sequence_wt) == 0:
                            sequence_wt.append("No PASS")
                            sequence_mut.append("No PASS")
                            sequence_type.append("No PASS")
                            sequence_pos.append(R_sim[0].strip().split("_")[1])
                            sequence_num.append("No PASS")
                        for (x, y, z, k, s) in zip(sequence_wt, sequence_mut, sequence_type, sequence_pos, sequence_num):
                            w1.write(k + "\t" + str(file.strip().split(".")[0]) + "\t" + x + "\t" + y + "\t" + z + "\t" +
                                     str(s) + "\n")
        self.end()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'crispr_risearch_dir':
            pass
        elif event['data'].startswith("crispr_resso_dir"):
            if event['data'].strip().split("_dir_")[1] == "wt":
                self.softdir(obj.output_dir, "wt")
            elif event['data'].strip().split("_dir_")[1] == "mut":
                self.softdir(obj.output_dir, "mut")
        elif event['data'].startswith("crispr_stat_dir"):
            self.linkdir(obj.output_dir, 'crispr_stat_dir')

    def softdir(self, dirpath, sample_type):
        """
        此处用于建立软链接
        :param dirpath:
        :param dirname:
        :param sample_type:
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = ""
        if sample_type == "wt":
            newdir = os.path.join(self.output_dir, "wt")
        elif sample_type == "mut":
            newdir = os.path.join(self.output_dir, "mut")
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
            os.symlink(oldfiles[i], newfiles[i])

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)    # 目标文件夹下所有的文件
        newdir = os.path.join(self.output_dir, dirname)  # 确定新目录地址
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]  # 所有老文件的路径
        newfiles = [os.path.join(newdir, i) for i in allfiles]   # 所有新文件的路径
        for newfile in newfiles:                                 # 删除新文件夹中已经有的文件，为复制提供地方。
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)

        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        """
        保存结果到mongo
        """
        if self.option("main_id"):
            crispr_analysis = self.api.api("wgs_v2.crispr_analysis")
            crispr_analysis.add_modified_effiency(self.output_dir + "/on_target.txt", self.option("gene_name"),
                                                  json.loads(self.option("bam_list")), self.option("main_id"), self.option("task_id"))
            crispr_analysis.add_crispr_stat(self.output_dir + "/on_target.txt", self.option("main_id"), self.option("task_id"))
            crispr_analysis.add_crispr_off_target(self.output_dir + "/off_target.txt", self.option("main_id"), self.option("task_id"))

    def run(self):
        super(CrisprAnalysisModule, self).run()
        self.get_files()
        self.run_risearch()
        # self.run_results()

    def end(self):
        if os.path.exists(self.output_dir + "/on_target.txt"):
            self.set_db()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(CrisprAnalysisModule, self).end()
