# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180404

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
import gevent.subprocess as subprocess
from biocluster.config import Config


class Pre4sentieonModule(Module):
    """
    用于对bam文件进行call突变位点，该模块中包含了三种方式的call 突变位点，分别是gatks，samtools，freebayes，其中进行call的
    时候是要将所有的bam文件合并在一起，然后将合并后的基因组分成25份进行并行计算，计算完成后再合并到一起
    last modified by fwy @ 20190620 改写为ref_rna_v2专用
    """
    def __init__(self, work_id):
        super(Pre4sentieonModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "in_bam", "type": "infile", "format": "ref_rna_v2.common"},  # 输入用于做前处理的bam文件
            {"name": "call_type", "type": "string", "default": "sentieon"},  # call snp的方式
            {"name": "analysis_format", "type": "string", "default": "bam"},
            {"name": "align_method", "type": "string", "default": "hisat"}
        ]
        self.add_option(options)
        self.file_list = []
        self.ref_fasta=""
        self.ref_index1=""
        self.ref_index2=""
        self.samtools_path = Config().SOFTWARE_DIR + "miniconda2/bin/"
        self.picard_path = Config().SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.java_path=Config().SOFTWARE_DIR+"/program/sun_jdk1.8.0/bin/java"


    def check_options(self):
        if not self.option("ref_fasta").is_set:
            raise OptionError("缺少ref_fasta参数")
        if not self.option("in_bam"):
            raise OptionError("请输入bam_文件")
        if self.option("call_type"):
            if self.option("call_type") not in ["samtools", "gatk", "freebayes", "sentieon"]:
                raise OptionError("call_type方式不合法，必须是samtools, gatk, freebayes, sentieon", code="25600307")
        else:
            raise OptionError("缺少call_type参数", code="25600308")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".fa")[0]
        #if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                #or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            #raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="24500406")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def genome_split(self):
        """
        将所有样本的bam总文件，进行分割，然后进行并行call snp
        :return:
        """
        ref_genome_name = os.path.basename(self.option("ref_fasta").prop["path"])
        if ref_genome_name.endswith(".fasta") or ref_genome_name.endswith(".fa"):
            self.ref_fasta = self.work_dir + "/" + os.path.basename(self.option("ref_fasta").prop["path"])
        else:
            self.ref_fasta = self.work_dir + "/" + os.path.basename(self.option("ref_fasta").prop["path"]) + '.fa'

        if os.path.exists(self.ref_fasta):
            os.remove(self.ref_fasta)
        if os.path.exists(self.ref_fasta + '.fai'):
            os.remove(self.ref_fasta + '.fai')
        if os.path.exists(self.work_dir+"/"+".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1])+ ".dict" ):
            os.remove(self.work_dir+"/"+".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1])+ ".dict" )
        os.link(self.option("ref_fasta").prop["path"], self.ref_fasta)
        self.ref_index1 = self.ref_fasta + ".fai"
        self.ref_index2 = self.work_dir+"/"+".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1])+ ".dict"
        if not os.path.exists(self.option("ref_fasta").prop["path"] + ".fai"):
            cmd = "{} faidx {}".format(self.samtools_path, self.option("ref_fasta").prop["path"])
            self.logger.info("开始进行samtools建索引！")
            subprocess.check_call(cmd, shell=True)
            os.link(self.option("ref_fasta").prop["path"] + '.fai', self.ref_index1)
        else:
            os.link(self.option("ref_fasta").prop["path"] + '.fai', self.ref_index1)
        if not os.path.exists(os.path.dirname(self.option("ref_fasta").prop["path"]) + "/" + ".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1]) + ".dict"):
            dirname = os.path.split(self.option("ref_fasta").prop["path"])[-1]
            dictname = os.path.splitext(dirname)[0] + ".dict"
            self.dict(self.ref_fasta, os.path.join(os.path.split(self.option("ref_fasta").prop["path"])[0],dictname))
            os.link(os.path.join(os.path.split(self.option("ref_fasta").prop["path"])[0], dictname), self.ref_index2)
        else:
            os.link(os.path.dirname(self.option("ref_fasta").prop["path"]) + "/" + ".".join(os.path.basename(self.option("ref_fasta").prop["path"]).split(".")[:-1]) + ".dict", self.ref_index2)
        n = 0
        split_file = "{}/split_file".format(self.work_dir)
        if os.path.exists(split_file) or os.path.isdir(split_file):
            os.system('rm -r %s' % split_file)
            os.mkdir(split_file)
        else:
            os.mkdir(split_file)
        with open(self.ref_index2, "r") as r:
            data = r.readlines()[1:]
            for line in data:
                temp = line.strip().split("\t")
                chr_name = temp[1].strip().split(":")[1]
                length = temp[2].strip().split(':')[1]
                n += 1
                file_name = n % 25
                if self.option("call_type") == "gatk":
                    with open(split_file + "/{}.intervals".format(file_name), 'a+') as w:
                        w.write(chr_name + "\n")
                else:
                    with open(split_file + "/{}.bed".format(file_name), 'a+') as w:
                        w.write(chr_name + "\t" + "0\t" + length + "\n")
        self.get_bed_intervals_file(split_file)

    def dict(self, ref_fasta, dict_name):
        """
        使用picard对参考基因组构建字典
        """

        cmd = "{} -jar {}picard.jar CreateSequenceDictionary R={} O={}" \
            .format(self.java_path,self.picard_path, ref_fasta, dict_name)
        if os.path.exists(dict_name):
            os.remove(dict_name)
        print cmd
        self.logger.info("开始用picard对参考基因组构建字典")
        subprocess.check_call(cmd, shell=True)
        self.logger.info("参考基因组构建dict done!")

    def get_bed_intervals_file(self, split_file):
        """
        获取分割后的bed文件或者intervals文件的个数
        :return:
        """
        files = os.listdir(split_file)
        if len(files) == 0:
            self.set_error("基因组分割后的列表为空，不能进行后续分析下，请检查ref.dict文件！", code="25600307")
        else:
            for m in files:
                self.file_list.append(os.path.join(split_file, m))
            self.logger.info("file_list为{}".format(self.file_list))

    def run_bam_addrg(self):
        self.bam_addrg=self.add_tool("ref_rna_v3.snp.bam_addrg")
        self.bam_addrg.set_options({
            "in_bam": self.option("in_bam").prop["path"],
            "method": "picard",
            "align_method": self.option("align_method")
        })
        if self.option("analysis_format").lower() == "bam":
            self.bam_addrg.on('end', self.run_bam_sort_index)
        elif self.option("analysis_format").lower() == "cram":
            self.bam_addrg.on('end', self.run_bam2cram)
        self.bam_addrg.run()

    def run_bam_sort_index(self):
        self.bam_sort_index=self.add_tool("ref_rna_v3.snp.bam_sort")
        if self.option("analysis_format").lower() == "bam":
            in_bam=self.bam_addrg.option("out_bam")
        elif self.option("analysis_format").lower() == "cram":
            in_bam = self.bam2cram.option("out_bam")
        self.logger.info("未sortbam为{}".format(in_bam))
        self.bam_sort_index.set_options({
            "in_bam": in_bam,
            "method": "samtools",
            "file_format":self.option("analysis_format")
        })
        self.bam_sort_index.on('end',self.run_sentieon_haplptyper)
        self.bam_sort_index.run()

    def run_bam2cram(self):
        self.bam2cram = self.add_tool("ref_rna_v3.snp.bam2cram")
        self.bam2cram.set_options({
            "in_bam": self.bam_addrg.option("out_bam"),
            "fa_file": self.option("ref_fasta").prop["path"],
        })
        self.bam2cram.on("end",self.run_bam_sort_index)
        self.bam2cram.run()


    def run_sentieon_haplptyper(self):
        self.sentieon_haplptyper=self.add_tool("ref_rna_v3.snp.sentieon_haplotyper")
        self.sample_name = os.path.splitext(os.path.basename(self.option("in_bam").prop["path"]))[0]
        self.sentieon_haplptyper.set_options({
            "bam_file": self.bam_sort_index.option("out_bam"),
            "fa_file": self.option("ref_fasta").prop["path"],
            "file_format": self.option("analysis_format"),
            "name":self.sample_name
        })
        self.sentieon_haplptyper.on("end", self.set_output)
        self.sentieon_haplptyper.run()


    def set_output(self):
        allfiles = os.listdir(self.sentieon_haplptyper.output_dir)
        resultfiles = [os.path.join(self.sentieon_haplptyper.output_dir, i) for i in allfiles]
        outfiles = [os.path.join(self.output_dir, i) for i in allfiles]
        for outfile in outfiles:
            if os.path.exists(outfile):
                if os.path.isfile(outfile):
                    os.remove(outfile)
                else:
                    os.system('rm -r %s' % outfile)
                # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            os.link(resultfiles[i], outfiles[i])
        self.end()


    def run(self):
        super(Pre4sentieonModule, self).run()
        self.logger.info("开始准备参考基因组索引文件及相关切割工作")
        self.genome_split()
        self.logger.info("参考基因组索引文件准备及切割工作完成")
        self.run_bam_addrg()
        # if self.option("call_type") == "samtools":
        #     self.samtools_call_run()
        # elif self.option("call_type") == "gatk":
        #     self.haplotype_run()
        # elif self.option("call_type") == "freebayes":
        #     self.freebayes_call_run()
        # elif self.option("call_type") == "sentieon":
        #     self.haplotype_v2_run()

    def set_bam_list(self):
        """
        samtools和freebayes与gatk的bamlist格式不一致，这里进行处理下
        """
        self.logger.info("开始设置bamlist！")
        with open(self.option("bam_list"), "r") as r, open("{}/bam.list".format(self.work_dir), "w") as w:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                w.write(temp[1] + "\n")
        self.logger.info("设置bamlist成功！")

    def end(self):
        super(Pre4sentieonModule, self).end()



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "ref_rna_v3.pre4sentieon",
            "instant": False,
            "options": dict(
                ref_fasta="/mnt/ilustre/users/sanger-dev/workspace/20200114/Single_sentieonnew_majorbio_2319378889/CallSnpIndelSentieon/Pre4sentieon2/Lcu_new_hic.fasta",
                in_bam="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test6/Refrna_majorbio_231937/output/bam/LC_PE_2.bam",
                call_type="sentieon",
                analysis_format="bam"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()