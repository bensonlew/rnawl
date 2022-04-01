# -*- coding: utf-8 -*-
# __author__ = zhujuan
# last_modify：2018.01.30


import re, Bio, urllib2, regex, os
import subprocess
from biocluster.iofile import File
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError

'''
gbk:GenBank
格式说明地址：https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
'''


class GbkFile(File):
    """
    定义Genbank文件
    """

    def __init__(self):
        super(GbkFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(GbkFile, self).get_info()
        seqinfo = self.get_gbk_info()
        self.set_property("file_format", seqinfo[0])
        self.set_property("genome_name", seqinfo[1])
        self.set_property("genome_length", seqinfo[2])
        self.set_property("genome_type", seqinfo[3])
        # self.set_property("genome_source", seqinfo[4])

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        self.get_info()
        if super(GbkFile, self).check():
            if self.prop['file_format'] != '.gbk':
                raise FileError("文件格式错误,必须以\".gbk\"结尾的文件！", code="42200301")
            if self.prop["genome_name"]:
                print "########基因组名称：" + self.prop["genome_name"]
        return True

    def get_gbk_info(self):
        """
        获取GenBank信息
        :return: (file_format,genome_name,genome_length,genome_type,gene_number)
        """
        file_format = os.path.splitext(os.path.basename(self.prop['path']))[1]
        with open(self.prop['path'], 'rb') as r:
            file = r.readlines()
            locus = re.match("^LOCUS\s+(\S+)\s+(\w+) bp\s+\w+\s+(\w+).*$", file[0])
            locus_nu = 0
            gene_nu = 0
            cds_nu = 0
            genome_name = ""
            genome_length = ""
            genome_type = ""
            feature = ""
            if locus:
                genome_name = locus.group(1)
                genome_length = locus.group(2)
                genome_type = locus.group(3)
                locus_nu = 1
            else:
                raise FileError("文件格式错误，Genbank文件应该含有具体LOCUS信息", code="42200302")
            for f in file[1:]:
                # if re.match("^SOURCE\s+(.*)$", f):
                #    genome_source = re.match("^SOURCE\s+(.*)$", f).group(1)
                if re.match("^FEATURES\s+.*$", f):
                    feature = "ok"
                if re.match("^ {5}CDS.*$", f):
                    gene_nu = gene_nu + 1
                if re.match("\s+/translation=\"[-A-Z].*$", f):
                    cds_nu = cds_nu + 1
                locus = re.match("^LOCUS\s+(\w+)\s+(\w+) bp\s+\w+\s+(\w+).*$", f)
                if locus:
                    locus_nu = locus_nu + 1
                    genome_name = genome_name + "," + locus.group(1)
                    genome_length = str(genome_length) + "," + str(locus.group(2))
                    genome_type = genome_type + "," + locus.group(3)
                    # print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                    # print str(gene_nu) + "\t" + str(cds_nu) + "\t" + feature
                    #if gene_nu != cds_nu or feature == "":  gbk文件中可能存在假基因（没有translation），所以不能做这个卡控
                    #   raise FileError("文件格式错误，Genbank文件缺少CDS信息或者CDS信息有误", code="42200303")
                    if feature == "":
                        raise FileError("文件格式错误，缺少FEATURES信息")
                    #if cds_nu == 0:
                    #    raise FileError("gbk文件格式错误，缺少CDS信息")
                    gene_number = ""
                    gene_number = gene_number + "\t" + str(gene_nu)
                    gene_nu = 0
                    cds_nu = 0
                    # if locus_nu == 1:
                    # gene_number = gene_nu
                    # print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                    # print str(gene_nu) + "\t" + str(cds_nu) + "\t" + feature
                    # if gene_number != cds_nu or feature == "":
                    #    raise FileError("文件格式错误，Genbank文件缺少CDS信息或者CDS信息有误")
        return file_format, genome_name, genome_length, genome_type

    def check_cds_info(
            self):  # 当需要从genbank文件中提取基因组序列信息时使用，检查CDS序列是否有对应的protein_id和product，没有则报错(如细菌基因组上传GBK文件作为参考基因组使用)
        with open(self.prop['path'], 'r') as t:
            product_n = 0
            protein_n = 0
            cds_n = 0
            file = t.readlines()
            for line in file:
                if "       /product=" in line:
                    product_n = product_n + 1
                elif "       /protein_id=" in line:
                    protein_n = protein_n + 1
                elif re.match(r'     CDS', line):
                    cds_n = cds_n + 1
            # if product_n == protein_n == cds_n:
            if protein_n == cds_n:
                pass
            else:
                print "cds_n:" + str(cds_n) + "\nproduct_n:" + str(product_n) + "\nprotein_n:" + str(protein_n)
                #raise FileError("文件格式错误，Genbank文件缺少product或者protein_id", code="42200304")  一个cds 和 protein 不能一一对应。允许这种情况情况


if __name__ == "__main__":
    gbk = GbkFile()
    # gbk.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/antismash/two.gbk")
    gbk.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/database/CM001165.gbk-bak")
    # test = gbk.get_info()
    # gbk.check()
    gbk.check_cds_info()
