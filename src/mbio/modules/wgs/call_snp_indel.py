# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180404

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os


class CallSnpIndelModule(Module):
    """
    用于对bam文件进行call突变位点，该模块中包含了三种方式的call 突变位点，分别是gatks，samtools，freebayes，其中进行call的
    时候是要将所有的bam文件合并在一起，然后将合并后的基因组分成25份进行并行计算，计算完成后再合并到一起
    last modified by hongdong @ 20190215 add sentieon method
    """
    def __init__(self, work_id):
        super(CallSnpIndelModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "bam_list", "type": "string"},  # 如果是gatk 样本+建库类型+bam文件，否则样本bam列表
            {"name": "ref_dict", "type": "string"},  # 参考基因组的配置文件，里面可以解析出有多少条染色体
            {"name": "call_type", "type": "string", "default": "samtools"}  # call snp的方式
        ]
        self.add_option(options)
        self.file_list = []
        self.samtools_call_tools = []
        self.freebayes_call_tools = []
        self.haplotype_tools = []
        self.gvcf_typing_tools = []
        self.bcftools_vcf = self.add_tool("wgs.bcftool_vcf")
        self.gvcf_typing_v2 = None

    def check_options(self):
        if not self.option("ref_fasta").is_set:
            raise OptionError("缺少ref_fasta参数", code="24500401")
        if not self.option("bam_list"):
            raise OptionError("缺少bam_list参数", code="24500402")
        if not self.option("ref_dict"):
            raise OptionError("缺少ref_dict参数", code="24500403")
        if self.option("call_type"):
            if self.option("call_type") not in ["samtools", "gatk", "freebayes", "sentieon"]:
                raise OptionError("call_type方式不合法，必须是samtools, gatk, freebayes, sentieon", code="24500404")
        else:
            raise OptionError("缺少call_type参数", code="24500405")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".")[0]
        if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="24500406")
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
        n = 0
        split_file = "{}/split_file".format(self.work_dir)
        if os.path.exists(split_file) or os.path.isdir(split_file):
            os.system('rm -r %s' % split_file)
            os.mkdir(split_file)
        else:
            os.mkdir(split_file)
        with open(self.option("ref_dict"), "r") as r:
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

    def get_bed_intervals_file(self, split_file):
        """
        获取分割后的bed文件或者intervals文件的个数
        :return:
        """
        files = os.listdir(split_file)
        if len(files) == 0:
            self.set_error("基因组分割后的列表为空，不能进行后续分析下，请检查ref.dict文件！", code="24500413")
        else:
            for m in files:
                self.file_list.append(os.path.join(split_file, m))
            self.logger.info("file_list为{}".format(self.file_list))

    def samtools_call_run(self):
        """
        按照work_dir中的split_file有多少个文件就运行多少个tool
        :return:
        """
        self.set_bam_list()
        n = 0
        for m in self.file_list:
            samtools_call = self.add_tool("wgs.samtools_call")
            samtools_call.set_options({
                "ref_fasta": self.option("ref_fasta").prop["path"],
                "bed_file": m,
                "bam_list": "{}/bam.list".format(self.work_dir)
            })
            self.samtools_call_tools.append(samtools_call)
            n += 1
        for j in range(len(self.samtools_call_tools)):
            self.samtools_call_tools[j].on('end', self.set_output, 'snp_indel')
        if self.samtools_call_tools:
            if len(self.samtools_call_tools) > 1:
                self.on_rely(self.samtools_call_tools, self.bcftools_vcf_run)
            elif len(self.samtools_call_tools) == 1:
                self.samtools_call_tools[0].on('end', self.bcftools_vcf_run)
            for tool in self.samtools_call_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("samtools_call_tools列表为空！", code="24500414")

    def freebayes_call_run(self):
        self.set_bam_list()
        n = 0
        for m in self.file_list:
            freebayes_call = self.add_tool("wgs.freebayes_call")
            freebayes_call.set_options({
                "ref_fasta": self.option("ref_fasta").prop["path"],
                "bed_file": m,
                "bam_list": "{}/bam.list".format(self.work_dir)
            })
            self.freebayes_call_tools.append(freebayes_call)
            n += 1
        for j in range(len(self.freebayes_call_tools)):
            self.freebayes_call_tools[j].on('end', self.set_output, 'snp_indel')
        if self.freebayes_call_tools:
            if len(self.freebayes_call_tools) > 1:
                self.on_rely(self.freebayes_call_tools, self.bcftools_vcf_run)
            elif len(self.freebayes_call_tools) == 1:
                self.freebayes_call_tools[0].on('end', self.bcftools_vcf_run)
            for tool in self.freebayes_call_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("freebayes_call_tools列表为空！", code="24500415")

    def haplotype_run(self):
        sample_info = self.get_sample_bam()
        n = 0
        for keys in sample_info.keys():
            haplotype = self.add_tool("wgs.haplotype")
            haplotype.set_options({
                "ref_fasta": self.option("ref_fasta").prop["path"],
                "sample_bam": sample_info[keys],
                "sample_id": keys
            })
            self.haplotype_tools.append(haplotype)
            n += 1
        for j in range(len(self.haplotype_tools)):
            self.haplotype_tools[j].on('end', self.set_output, 'haplotype')
        if self.haplotype_tools:
            if len(self.haplotype_tools) > 1:
                self.on_rely(self.haplotype_tools, self.gvcf_typing_run)
            elif len(self.haplotype_tools) == 1:
                self.haplotype_tools[0].on('end', self.gvcf_typing_run)
            for tool in self.haplotype_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("haplotype_tools列表为空！", code="24500416")

    def gvcf_typing_run(self):
        n = 0
        for m in self.file_list:
            gvcf_typing = self.add_tool("wgs.gvcf_typing")
            gvcf_typing.set_options({
                "ref_fasta": self.option("ref_fasta").prop["path"],
                "gvcf_path": self.output_dir + "/haplotype",
                "intervals_file": m
            })
            self.gvcf_typing_tools.append(gvcf_typing)
            n += 1
        for j in range(len(self.gvcf_typing_tools)):
            self.gvcf_typing_tools[j].on('end', self.set_output, 'snp_indel')
        if self.gvcf_typing_tools:
            if len(self.gvcf_typing_tools) > 1:
                self.on_rely(self.gvcf_typing_tools, self.bcftools_vcf_run)
            elif len(self.gvcf_typing_tools) == 1:
                self.gvcf_typing_tools[0].on('end', self.bcftools_vcf_run)
            for tool in self.gvcf_typing_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("gvcf_typing_tools列表为空！", code="24500417")

    def haplotype_v2_run(self):
        sample_info = self.get_sample_bam()
        n = 0
        for keys in sample_info.keys():
            haplotype = self.add_tool("wgs_v2.haplotyper_v2")
            haplotype.set_options({
                "fa_file": self.option("ref_fasta").prop["path"],
                "bam_file": sample_info[keys],
                "name": keys + '.g.vcf'
            })
            self.haplotype_tools.append(haplotype)
            n += 1
        for j in range(len(self.haplotype_tools)):
            self.haplotype_tools[j].on('end', self.set_output, 'haplotype')
        if self.haplotype_tools:
            if len(self.haplotype_tools) > 1:
                self.on_rely(self.haplotype_tools, self.gvcf_typing_v2_run)
            elif len(self.haplotype_tools) == 1:
                self.haplotype_tools[0].on('end', self.gvcf_typing_v2_run)
        else:
            raise Exception("haplotype_tools列表为空！")
        for tool in self.haplotype_tools:
            gevent.sleep(1)
            tool.run()

    def gvcf_typing_v2_run(self):
        """
        DE1_10  /mnt/ilustre/users/caixia.tian/newmdt/Project/MJ20180723036-gaojun/03.ref/08.haplotype/DE1_10.g.vcf
        :return:
        """
        self.make_gvcf_list()
        self.gvcf_typing_v2 = self.add_tool("wgs_v2.gvcf_typing_v2")
        self.gvcf_typing_v2.set_options({
            "vcf_list": self.output_dir + "/haplotype/gvcf.list",
            "fa_file": self.option("ref_fasta").prop["path"]
        })
        self.gvcf_typing_v2.on('end', self.set_output, 'vcf_call')
        self.gvcf_typing_v2.on('end', self.end)
        self.gvcf_typing_v2.run()

    def make_gvcf_list(self):
        """
        DE1_10  /mnt/ilustre/users/caixia.tian/newmdt/Project/MJ20180723036-gaojun/03.ref/08.haplotype/DE1_10.g.vcf
        :return:
        """
        with open(self.output_dir + "/haplotype/gvcf.list", 'w') as w:
            for m in os.listdir(self.output_dir + "/haplotype"):
                n = re.match('(.*)\.g\.vcf$', m)
                if n:
                    w.write('{}\t{}\n'.format(n.group(1), os.path.join(self.output_dir + "/haplotype", m)))

    def get_sample_bam(self):
        """
        生成样本与bam文件对应的字典文件，这里只是gatk call snp需要， bam文件来源于04.bam-sort或者05.bam-mkdup文件夹
        :return:
        """
        sample_info = {}
        with open(self.option("bam_list"), "r") as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                if len(temp) == 3:
                    if temp[1].upper() not in ['WGS', 'RAD', 'GBS', 'WES']:
                        self.logger.warn("{}建库类型不合法!".format(temp[1]))
                    sample_info[temp[0]] = temp[2]
                elif len(temp) == 2:
                    sample_info[temp[0]] = temp[1]
                else:
                    self.set_error("%sbam_list格式不正确！", variables=(self.option("bam_list")), code="24500418")
        self.logger.info("样本与bam文件对应关系：{}".format(sample_info))
        return sample_info

    def get_vcf_list(self):
        gevent.sleep(2)
        path = os.path.join(self.output_dir, 'snp_indel')
        files = os.listdir(path)
        vcf = os.path.join(self.output_dir, "vcf.list")
        with open(vcf, 'w') as w:
            for m in files:
                if self.option("call_type") == "samtools":
                    if re.match(r".*\.bed\.vcf\.gz$", m):
                        w.write("{}\n".format(os.path.join(path, m)))
                elif self.option("call_type") == "gatk":
                    if re.match(r".*\.noid\.vcf$", m):
                        w.write("{}\n".format(os.path.join(path, m)))
                elif self.option("call_type") == "freebayes":
                    if re.match(r".*\.bed\.vcf$", m):
                        w.write("{}\n".format(os.path.join(path, m)))
                else:
                    pass

    def bcftools_vcf_run(self):
        self.get_vcf_list()
        self.bcftools_vcf.set_options({
            "vcf_list": os.path.join(self.output_dir, "vcf.list")
        })
        self.bcftools_vcf.on('end', self.set_output, 'vcf_call')
        self.bcftools_vcf.on("end", self.end)
        self.bcftools_vcf.run()

    def linkdir(self, dirpath, dirname):
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
                # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'snp_indel':
            self.linkdir(obj.output_dir, 'snp_indel')
        elif event['data'] == 'haplotype':
            self.linkdir(obj.output_dir, 'haplotype')
        elif event['data'] == 'gvcf_typing':
            self.linkdir(obj.output_dir, 'gvcf_typing')
        elif event['data'] == "vcf_call":
            self.linkdir(obj.output_dir, 'vcf_call')
        else:
            pass

    def run(self):
        super(CallSnpIndelModule, self).run()
        self.logger.info("开始切割")
        self.genome_split()
        self.logger.info("切割完成")
        if self.option("call_type") == "samtools":
            self.samtools_call_run()
        elif self.option("call_type") == "gatk":
            self.haplotype_run()
        elif self.option("call_type") == "freebayes":
            self.freebayes_call_run()
        elif self.option("call_type") == "sentieon":
            self.haplotype_v2_run()

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
        super(CallSnpIndelModule, self).end()
