# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""OTUtable格式文件类"""

from biocluster.iofile import File
import subprocess
from biocluster.config import Config
import os
import re
from biocluster.core.exceptions import FileError
from collections import defaultdict


class OtuTableFile(File):
    """
    OTUtable
    """
    def __init__(self):
        """
        """
        super(OtuTableFile, self).__init__()
        self.biom_path = os.path.join(Config().SOFTWARE_DIR, "program/Python/bin/")
        self.otu2shared_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/meta/scripts/otu2shared.pl")

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(OtuTableFile, self).get_info()
        info = self.get_otuinfo()
        self.set_property("form", info[0])
        self.set_property("otu_num", info[1])
        self.set_property("sample_num", info[2])
        self.set_property("metadata", info[3])

    def check(self):
        """
        检测文件是否满足要求
        :return:
        """
        if super(OtuTableFile, self).check():
            self.get_info()
            if self.prop['form']:
                pass
            else:
                raise FileError("文件格式错误", code="42700801")
        return True

    def get_otuinfo(self):
        """
        获取otu表信息
        """
        form, otu_num, sample_num, metadata = True, 0, 0, ''
        with open(self.prop['path'], 'r') as f:
            heads = f.readline().rstrip().split('\t')
            colnum = len(heads)
            # if not re.match(r'#*OTU ID', heads[0]):
            #     form = False
            if colnum < 2:
                form = False
            if form:
                sample_num = colnum - 1
                if heads[colnum - 1] == 'taxonomy':
                    metadata = 'taxonomy'
                    sample_num = colnum - 2
                while 1:
                    line = f.readline().rstrip()
                    if not line:
                        break
                    otu_num += 1
        return (form, otu_num, sample_num, metadata)

    def convert_to_biom(self, biom_filepath):
        """
        转换为biom格式
        """
        cmd = self.biom_path + "biom convert -i " + self.prop['path'] + " -o " + biom_filepath + ' --table-type \"OTU table\" --to-hdf5'
        if 'metadata' not in self.prop:
            self.get_info()
        if self.prop['metadata'] == "taxonomy":
            cmd += " --process-obs-metadata taxonomy"
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            raise FileError("biom convert 运行出错！", code="42700802")
        return True

    def convert_to_shared(self, shared_filepath):
        """
        转换为mothur的shared格式
        """
        # otu2shared.pl -i otutable -l 0.97 -o otu.shared
        if self.prop['metadata'] == "taxonomy":
            raise FileError(u"can not covert otutable with taxon info.", code="42700803")
        cmd = self.otu2shared_path + " -l 0.97 -i " + self.prop['path'] + " -o " + shared_filepath
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            raise FileError("otu2shared.pl 运行出错！",code="42700804")
        return True

    def ahead_taxonmoy(self, ahead_path, otu_path=None):
        if not otu_path:
            otu_path = self.prop["path"]
        with open(otu_path, 'rb') as r, open(ahead_path, 'wb') as w:
            """
            将taxonomy列信息补充道第一列当中
            """
            line = r.next()
            line = line.rstrip('\n')
            line = re.split('\t', line)
            if line[-1] != "taxonomy":
                raise FileError("文件最后一列不为taxonomy，无法将taxonomy列提前", code="42700805")
            line.pop()
            w.write("\t".join(line) + "\n")
            for line in r:
                line = line.rstrip('\n')
                line = re.split('\t', line)
                tax_info = line[-1]
                line.pop()
                line[0] = tax_info + "; " + line[0]
                newline = "\t".join(line)
                w.write(newline + "\n")

    def behind_taxnomy(self, behind_path, otu_path=None):
        """
        将文件的第一列信息拆分开来，只留下otu的信息，将OTU的信息放到每一行的最后，形成taxnomy列
        """
        if not otu_path:
            otu_path = self.prop["path"]
        with open(otu_path, 'rb') as r, open(behind_path, 'wb') as w:
            line = r.next().rstrip()
            line = re.split("\t", line)
            if line[-1] == "taxonomy":
                raise FileError("文件最后一列为taxonomy，无法将taxonomy列放到最后", code="42700806")
            w.write("\t".join(line) + "\ttaxonomy\n")
            for line in r:
                line = line.rstrip('\n')
                line = re.split('\t', line)
                sp = re.split("; ", line[0])
                line.pop(0)
                otu = sp.pop(-1)
                taxonomy = "; ".join(sp)
                newline = otu + "\t" + "\t".join(line) + "\t" + taxonomy
                w.write(newline + "\n")

    def del_taxonomy(self, otu_path, no_tax_otu_path):
        """
        删除taxonomy列  并返回一个字典
        """
        tax_dic = dict()
        with open(otu_path, 'rb') as r, open(no_tax_otu_path, 'wb') as w:
            line = r.next()
            line = line.rstrip('\n')
            line = re.split('\t', line)
            if line[-1] != "taxonomy":
                raise FileError("文件最后一列不为taxonomy，不需要去除taxonomy列", code="42700807")
            line.pop()
            w.write("\t".join(line) + "\n")
            for line in r:
                line = line.rstrip('\n')
                line = re.split('\t', line)
                tax_dic[line[0]] = line[-1]
                line.pop()
                newline = "\t".join(line)
                w.write(newline + "\n")
        return tax_dic

    def add_taxonomy(self, tax_dic, no_tax_otu_path, otu_path):
        """
        根据字典， 为一个otu表添加taxonomy列
        """
        if not isinstance(tax_dic, dict):
            raise Exception("输入的tax_dic不是一个字典")
        with open(no_tax_otu_path, 'rb') as r, open(otu_path, 'wb') as w:
            line = r.next().rstrip('\n')
            w.write(line + "\t" + "taxonomy\n")
            line = re.split('\t', line)
            if line[-1] == "taxonomy":
                raise Exception("输入文件已经含有taxonomy列")
            for line in r:
                line = line.rstrip("\n")
                w.write(line)
                line = re.split("\t", line)
                w.write("\t" + tax_dic[line[0]] + "\n")

    def complete_taxonomy(self, otu_path, complete_path):
        """
        将一个OTU表的taxnomy列补全
        """
        with open(otu_path, 'rb') as r, open(complete_path, 'wb') as w:
            line = r.next().rstrip('\n')
            w.write(line + "\n")
            line = re.split('\t', line)
            if line[-1] != "taxonomy":
                raise Exception("文件taxonomy信息缺失，请调用complete_taxonomy_by_dic方法补全taxonomy")
            for line in r:
                line = line.rstrip("\n")
                line = re.split('\t', line)
                tax = line.pop()
                info = "\t".join(line)
                new_tax = self._comp_tax(tax)
                w.write(info + "\t" + new_tax + "\n")

    def complete_taxonomy_by_dic(self, tax_dic, otu_path, complete_path):
        """
        根据字典, 添加taxnomy列并补全
        """
        if not isinstance(tax_dic, dict):
            raise Exception("输入的tax_dic不是一个字典")
        with open(otu_path, 'rb') as r, open(complete_path, 'wb') as w:
            line = r.next().rstrip('\n')
            w.write(line + "\t" + "taxonomy\n")
            line = re.split('\t', line)
            if line[-1] == "taxonomy":
                raise Exception("输入文件已经含有taxonomy列, 请调用complete_taxonomy方法补全taxonomy")
            for line in r:
                line = line.rstrip("\n")
                w.write(line)
                line = re.split('\t', line)
                tax = tax_dic[line[0]]
                new_tax = self._comp_tax(tax)
                w.write("\t" + new_tax + "\n")

    @staticmethod
    def _comp_tax(tax):
        LEVEL = {
            0: "d__", 1: "k__", 2: "p__", 3: "c__", 4: "o__",
            5: "f__", 6: "g__", 7: "s__"
        }
        tax = re.sub(r'\s', '', tax)
        cla = re.split(';', tax)
        last_tax = cla[-1]
        search_norank_parent = re.search(r'^[dkpcofgs]__Unclassified_([dkpcofgs]__.+$)', last_tax, flags=re.IGNORECASE)
        if search_norank_parent:
            last_tax = search_norank_parent.groups()[0]
        for i in range(len(cla), 8):
            cla.append(LEVEL[i] + "unclassified_" + last_tax)
        return ";".join(cla)

    def sub_otu_sample(self, samples, path,rmsame="F"):
        """
        从一张otu表里挑选几个样本， 在挑选这些样本之后，有些OTU的数目会变成0, 删除这些OTU,可选择同时删除在所有样品中数目完全一样的OTU  add by zouxuan for pca
        """
        self.get_info()
        tmp_dict = dict()
        my_path = self.prop["path"]
        if self.prop["metadata"] == "taxonomy":
            no_tax = os.path.join(os.path.dirname(self.prop["path"]), os.path.basename(self.prop["path"]) + ".no_tax")
            tmp_dict = self.del_taxonomy(self.prop["path"], no_tax)
            my_path = no_tax
        colnum_list = list()
        otu_dict = dict()
        new_head = list()
        if not isinstance(samples, list):
            raise Exception("samples参数必须是列表")
        with open(my_path, 'rb') as r:
            head = r.next().strip("\r\n")
            head = re.split('\t', head)
            for i in range(1, len(head)):
                if head[i] in samples:
                    colnum_list.append(i)
                    new_head.append(head[i])
            for line in r:
                line = line.strip('\r\n')
                line = re.split('\t', line)
                otu_dict[line[0]] = list()
                sum_ = 0
                for i in colnum_list:
                    otu_dict[line[0]].append(line[i])
                all= []
                for num in otu_dict[line[0]]:
                    sum_ += float(num)  # int(num)  modified by guhaidong for metagenome pipline at 2017.10.19
                    all.append(num)
                if sum_ == 0:
                    del otu_dict[line[0]]
                else:
                    if rmsame == "T":          #add by zouxuan for pca
                        if len(set(all)) == 1:
                            del otu_dict[line[0]]		    
        if self.prop["metadata"] == "taxonomy":
            out_path = os.path.join(os.path.dirname(self.prop["path"]), os.path.basename(self.prop["path"]) + ".no_zero")
        else:
            out_path = path
        with open(out_path, 'wb') as w:
            w.write("OTU ID\t")
            w.write("\t".join(new_head) + "\n")
            for otu in otu_dict:
                w.write(otu + "\t")
                w.write("\t".join(otu_dict[otu]) + "\n")
        if self.prop["metadata"] == "taxonomy":
            self.add_taxonomy(tmp_dict, out_path, path)

    def transposition(self, path):
        """
        转置一个otu表
        """
        file_ = list()
        with open(self.prop['path'], 'rb') as r:
            linelist = [l.strip('\r\n') for l in r.readlines()]
        for row in linelist:
            row = re.split("\t", row)
            file_.append(row)
        zip_line = zip(*file_)
        with open(path, 'wb') as w:
            for my_l in zip_line:
                w.write("\t".join(my_l) + "\n")

    def get_sample_info(self):
        """
        获取otu表中样本信息，返回样本列表
        """
        with open(self.prop['path'], 'r') as f:
            sample = f.readline().strip().split('\t')
            del sample[0]
        return sample

    def get_min_sample_num(self):
        """
        获取最小的样本序列数
        """
        sample_name = list()
        sample_num = defaultdict(int)
        self.get_info()
        with open(self.prop['path'], 'rb') as r:
            line = r.next().strip().split("\t")
            line.pop(0)
            if self.prop["metadata"] == "taxonomy":
                line.pop(-1)
            sample_name = line[:]
            for line in r:
                line = line.rstrip().split("\t")
                line.pop(0)
                if self.prop["metadata"] == "taxonomy":
                    line.pop(-1)
                for i in range(len(sample_name)):
                    sample_num[sample_name[i]] += int(line[i])
            min_num = sample_num.values()[0]
            min_sample = sample_num.keys()[0]
            for k in sample_num:
                if sample_num[k] < min_num:
                    min_num = sample_num[k]
                    min_sample = k
        return (min_sample, min_num)

    def extract_info(self):
        """
        解析OTU表，将OTU表解析成一个二维字典
        例如 dict[OTU1][sample1] = 20
        代表sample1有20条序列在OTU1里
        """
        sample_name = list()
        info_dict = defaultdict(dict)
        self.get_info()
        with open(self.prop['path'], 'rb') as r:
            line = r.next().rstrip().split("\t")
            line.pop(0)
            if self.prop["metadata"] == "taxonomy":
                line.pop(-1)
            sample_name = line[:]
            for line in r:
                line = line.rstrip().split("\t")
                otu_name = line.pop(0)
                if self.prop["metadata"] == "taxonomy":
                    line.pop(-1)
                for i in range(len(sample_name)):
                    info_dict[otu_name][sample_name[i]] = int(line[i])
        return info_dict


if __name__ == "__main__":
    #print OtuTableFile._comp_tax('d__Bacteria; k__norank; p__Firmicutes; c__Erysipelotrichia; o__Erysipelotrichales; f__Erysipelotrichaceae; g__norank_f__Erysipelotrichaceae')
    data = OtuTableFile()
    taxon_otu = "/mnt/ilustre/users/sanger-dev/workspace/20170308/MetaBase_meta_test_nt_4444/output/OtuTaxon_summary/otu_taxon.xls"
    taxon_otu_new = "/mnt/ilustre/users/sanger-dev/workspace/20170308/MetaBase_meta_test_nt_4444/output/OtuTaxon_summary/otu_taxon_new.xls"
    data.set_path(taxon_otu)
    data.get_info()
    data.complete_taxonomy(taxon_otu,taxon_otu_new)
