# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'

import re
import os
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class BsaPathFile(Directory):
    """
    bsa特有路径检查，主要 检查对应文件存在不存在
    __author__ = HONGDONG
    __last_modify__ = 20180611
    对接wgs的目录结构
    """
    def __init__(self):
        super(BsaPathFile, self).__init__()
        self.summary = ""
        self.ref = ""
        self.vcf = ""
        self.is_old = True

    def info_check(self):
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            fastq_dir = os.path.join(self.prop['path'], "01.fastq_qc")
            ref_dir = os.path.join(self.prop['path'], "02.reference")
            map_dir = os.path.join(self.prop['path'], "03.map_stat")
            variant_dir = os.path.join(self.prop['path'], "04.snp_indel")
            annovar_dir = os.path.join(self.prop['path'], "05.annovar")
            if not os.path.exists(fastq_dir + '/qc_stat'):
                self.check_fastq_dir(fastq_dir)
            self.check_ref_dir(ref_dir)
            self.check_map_dir(map_dir)
            self.check_variant_dir(variant_dir)
            self.check_annovar_dir(annovar_dir)
            self.set_property("summary", self.summary)
            self.set_property("ref", self.ref)
            self.set_property("vcf", self.vcf)
        else:
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径--file error!", code="41500101")

    def check_fastq_dir(self, fastq_dir):
        """
        检查01.fastq_qc中文件及文件夹是否存在
        :param fastq_dir:
        :return:
        """
        atgc = fastq_dir + '/cleandata_qc/atgc'
        qual = fastq_dir + '/cleandata_qc/qual'
        stat = fastq_dir + '/cleandata_qc/qc.xls'
        self.is_dir(fastq_dir)
        self.is_dir(atgc)
        self.is_dir(qual)
        for f in os.listdir(atgc):
            self.is_column_num(file_path=atgc + "/" + f, number=6)
        if not os.path.exists(stat):
            raise FileError("缺少文件%s".format(stat), variables=(stat), code="41500102")
        self.is_column_num(file_path=stat, number=5)

    def check_ref_dir(self, ref_dir):
        """
        检查02.reference中文件及文件夹是否存在
        :param ref_dir:
        :return:
        """
        self.is_dir(ref_dir)
        self.ref = ref_dir + '/ref.fa'
        self.is_file(ref_dir + '/ref.fa')
        self.is_file(ref_dir + '/info.log')
        self.is_file(ref_dir + '/ref.chrlist')

    def check_map_dir(self, map_dir):
        """
        检查03.map_stat中文件及文件夹是否存在
        :param map_dir:
        :return:
        """
        self.is_dir(map_dir)
        coverage = map_dir + '/coverage'
        depth = map_dir + '/depth'
        insert = map_dir + '/insert'
        result = map_dir + '/result.stat'
        self.is_dir(coverage)
        self.is_dir(depth)
        self.is_dir(insert)
        self.is_dir(result)
        self.is_file(result + '/Total.mapped.detail.xls')
        for f in os.listdir(coverage):
            self.is_column_num(file_path=coverage + "/" + f, number=3)
        for f in os.listdir(depth):
            self.is_column_num(file_path=depth + "/" + f, number=2)
        for f in os.listdir(insert):
            self.is_column_num(file_path=insert + "/" + f, number=2)

    def check_variant_dir(self, variant_dir):
        """
        检查04.snp_indel中文件及文件夹是否存在
        :param variant_dir:
        :return:
        """
        self.is_dir(variant_dir)
        self.is_dir(variant_dir + '/variant_stat')
        self.is_dir(variant_dir + '/anno_stat')
        self.is_file(variant_dir + '/variant_stat/indel.len')
        self.is_file(variant_dir + '/variant_stat/snp.stat')
        self.is_file(variant_dir + '/variant_stat/indel.stat')
        self.is_file(variant_dir + '/anno_stat/indel.stat')
        self.is_file(variant_dir + '/anno_stat/snp.stat')

    def check_annovar_dir(self, annovar_dir):
        """
        检查05.annovar中文件及文件夹是否存在
        :param annovar_dir:
        :return:
        """
        self.is_dir(annovar_dir)
        self.is_file(annovar_dir + '/combine_variants/pop.final.vcf')
        self.vcf = annovar_dir + '/combine_variants/pop.final.vcf'
        self.is_file(annovar_dir + '/anno_count/pop.summary')
        self.summary = annovar_dir + '/anno_count/pop.summary'
        self.check_pop_summary(annovar_dir + '/anno_count/pop.summary')

    def info_check_old(self):
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            fastq_dir = os.path.join(self.prop['path'], "01.fastq_qc")
            ref_dir = os.path.join(self.prop['path'], "02.reference")
            map_dir = os.path.join(self.prop['path'], "03.map_stat")
            variant_dir = os.path.join(self.prop['path'], "04.variant-stat")
            annovar_dir = os.path.join(self.prop['path'], "05.annovar")
            self.check_fastq_dir_old(fastq_dir)
            self.check_ref_dir_old(ref_dir)
            self.check_map_dir_old(map_dir)
            self.check_variant_dir_old(variant_dir)
            self.check_annovar_dir_old(annovar_dir)
            self.set_property("summary", self.summary)
            self.set_property("ref", self.ref)
            self.set_property("vcf", self.vcf)
        else:
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径--file error!", code="41500103")

    def check_fastq_dir_old(self, fastq_dir):
        """
        检查01.fastq_qc中文件及文件夹是否存在
        :param fastq_dir:
        :return:
        """
        self.is_dir(fastq_dir)
        self.is_dir(fastq_dir + '/atgc')
        self.is_dir(fastq_dir + '/qual')
        self.is_dir(fastq_dir + '/stat')
        atgc = fastq_dir + '/atgc'
        qual = fastq_dir + '/qual'
        stat = fastq_dir + '/stat/QC.xls'
        for f in os.listdir(atgc):
            self.is_column_num(file_path=atgc + "/" + f, number=6)
        if not os.path.exists(stat):
            raise FileError("缺少文件%s", variables=(stat), code="41500104")
        self.is_column_num(file_path=stat, number=5)

    def check_ref_dir_old(self, ref_dir):
        """
        检查02.reference中文件及文件夹是否存在
        :param ref_dir:
        :return:
        """
        self.is_dir(ref_dir)
        self.ref = ref_dir + '/ref.fa.gz'
        self.is_file(ref_dir + '/ref.fa.gz')
        self.is_file(ref_dir + '/ref.log')
        self.is_file(ref_dir + '/ref.chrlist')

    def check_map_dir_old(self, map_dir):
        """
        检查03.map_stat中文件及文件夹是否存在
        :param map_dir:
        :return:
        """
        self.is_dir(map_dir)
        self.is_dir(map_dir + '/coverage')
        self.is_dir(map_dir + '/depth')
        self.is_dir(map_dir + '/insert')
        self.is_dir(map_dir + '/result.stat')
        # self.is_file(map_dir + '/result.stat/3-4.xls')
        # self.is_file(map_dir + '/result.stat/3-5.xls')
        self.is_file(map_dir + '/result.stat/Total.mapped.detail')
        coverage = map_dir + '/coverage'
        for f in os.listdir(coverage):
            self.is_column_num(file_path=coverage + "/" + f, number=3)
        depth = map_dir + '/depth'
        for f in os.listdir(depth):
            self.is_column_num(file_path=depth + "/" + f, number=2)
        insert = map_dir + '/insert'
        for f in os.listdir(insert):
            self.is_column_num(file_path=insert + "/" + f, number=2)

    def check_variant_dir_old(self, variant_dir):
        """
        检查04.variant-stat中文件及文件夹是否存在
        :param variant_dir:
        :return:
        """
        self.is_dir(variant_dir)
        self.is_dir(variant_dir + '/indel')
        self.is_dir(variant_dir + '/snp')
        self.is_file(variant_dir + '/snp/snp.stat')
        self.is_file(variant_dir + '/snp/snp.region')
        self.is_file(variant_dir + '/snp/snp.effects')
        self.is_file(variant_dir + '/indel/indel.stat')
        self.is_file(variant_dir + '/indel/indel.region')
        self.is_file(variant_dir + '/indel/indel.effects')

    def check_annovar_dir_old(self, annovar_dir):
        """
        检查05.annovar中文件及文件夹是否存在
        :param annovar_dir:
        :return:
        """
        self.is_dir(annovar_dir)
        self.is_file(annovar_dir + '/pop.final.vcf.gz')
        self.vcf = annovar_dir + '/pop.final.vcf.gz'
        self.is_file(annovar_dir + '/pop.summary')
        self.summary = annovar_dir + '/pop.summary'
        self.check_pop_summary(annovar_dir + '/pop.summary')

    def is_file(self, file_path):
        """
        检查是否是文件是否存在
        :param file_path:
        :return:
        """
        if not os.path.isfile(file_path) or not os.path.exists(file_path):
            raise FileError("原始文件中不存在%s文件！", variables=(file_path), code="41500105")

    def is_dir(self, dir_path):
        """
        检查文件夹是否存在
        :param dir_path:
        :return:
        """
        if not os.path.isdir(dir_path) or not os.path.exists(dir_path):
            raise FileError("原始文件中不存在%s路径！", variables=(dir_path), code="41500106")

    def is_column_num(self, file_path, number):
        """
        检查文件有多少列
        :param file_path:
        :param number:
        :return:
        """
        with open(file_path, "r") as r:
            for m in r:
                item = m.strip().split("\t")
                if len(item) != number:
                    raise FileError("%s文件列数不为%s，请检查", variables=(file_path, number), code="41500107")
                break

    def check_pop_summary(self, file_path):
        """
        用于检查pop summary是否有问题，重点是检查EggNOG-ID的格式是否正确，会影响到后面的eggnog的注释
        K,T:Transcription;Signal transduction mechanisms
        :param file_path:
        :return:
        """
        with open(file_path, "r") as r:
            data = r.readlines()
            for m in data:
                temp = m.strip().split("\t")
                # print len(temp)
                if re.match(r'^#', temp[0]):
                    if len(temp) not in [21, 25]:
                        raise FileError("pop_summary的列数不对--应该是21列或者25列！", code="41500108")
                    else:
                        continue
                else:
                    # print temp[-2]
                    if temp[19] == "--":
                        continue
                    else:
                        tem1 = temp[19].split(':')
                        if len(tem1) == 2:
                            if len(tem1[0].split(",")) == len(tem1[1].split(';')):
                                pass
                            elif len(tem1[0].split(";")) == len(tem1[1].split(';')):
                                pass
                            else:
                                raise FileError("pop_summary文件EggNOG-ID列：%s命名不正确！".format(temp[19]), variables=(temp[19]), code="41500109")
                        else:
                            raise FileError("pop_summary文件EggNOG-ID列：%s命名不正确！".format(temp[19]), variables=(temp[19]), code="41500110")

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(BsaPathFile, self).check():
            vcf = os.path.join(self.prop['path'], "05.annovar/combine_variants/pop.final.vcf")
            summary = os.path.join(self.prop['path'], "05.annovar/anno_count/pop.summary")
            ref = os.path.join(self.prop['path'], "02.reference/ref.fa")
            if os.path.exists(vcf) and os.path.exists(summary) and os.path.exists(ref):
                self.info_check()
            else:
                self.info_check_old()
            return True

if __name__ == "__main__":
    a = BsaPathFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/workspace/20180528/Wgs_tsg_30041/output")
    # a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/bsa/data_demo/Arabidopsis_thaliana")
    print a.check()
    # print a.summary
    # a.check_pop_summary("/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5aa763ec923bc/Gossypium_hirsutum/05.annovar/pop.summary")
    print "检查通过"