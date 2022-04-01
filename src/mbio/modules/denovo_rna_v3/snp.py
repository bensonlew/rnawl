
# -*- coding: utf-8 -*-

import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from Bio import SeqIO
import glob2
import gevent
import pandas as pd
import traceback
import gevent.subprocess as subprocess
import unittest

class SnpModule(Module):
    """
    denovo_rna_v2运行tool下面的snp.py和snpfinal.py
    """

    def __init__(self, work_id):
        super(SnpModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format":"denovo_rna_v2.trinity_fasta"},
            # {"name": "bamlist", "type": "infile", "format": "denovo_rna_v2.bamlist"},
            # {"name": "call_type", "type": "string", "default": "samtools"},  # call snp的方式
            {"name": "fq_type", "type": "string", "default": "PE"},  # fq类型，PE、SE
            {"name": "fq_list", "type": "string"}, # 文件里面根据单双端的不同有3列或者2列
            # {"name": "isoform_unigene", "type": "string"},#实际为一个文件，只是不检查
            {"name": "call_type", "type": "string", "default": "samtools"},
            {"name": "cds_bed", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "allt2g", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "anno", "type": "infile", "format": "denovo_rna_v2.common"},  # 注释结果文件
            {"name": "in_bam", "type": "string", 'default': None},

        ]
        self.sum_tools = []
        self.snp_tools = []
        self.bwa_tools = []
        self.add_option(options)
        self.step.add_steps("snp", "snpfinal")

    def check_options(self):
        """
        检查参数..
        :return:
        """
        # if not self.option('trinity_fa'):
        #     raise OptionError('必须输入trinity_fa作为比对参考文件', code = "22001201")
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def split_ref_run(self):
        split_file = "{}/split_file".format(self.work_dir)
        if os.path.exists(split_file) or os.path.isdir(split_file):
            os.system('rm -r %s' % split_file)
            os.mkdir(split_file)
        else:
            os.mkdir(split_file)
        i = 1
        line = 1
        line_limit = 10000
        w = open(split_file + '/fasta_1', 'wb')
        for seq_record in SeqIO.parse(self.option("ref_fasta").prop['path'], "fasta"):
            if line <= line_limit:
                w.write('{}\t{}\t{}\n'.format(seq_record.id, 0, len(seq_record.seq)))
                line += 1
            else:
                i += 1
                w.close()
                line = 1
                w = open(split_file + '/fasta_{}'.format(i), 'wb')
                w.write('{}\t{}\t{}\n'.format(seq_record.id, 0, len(seq_record.seq)))
                line += 1
        w.close()

    def bwa_run(self, fq_type, fq_list):
        n = 0
        with open(fq_list) as f:
            for line in f:
                self.bwa = self.add_tool("denovo_rna_v2.bwa")
                n += 1
                self.bwa.set_options({
                    "fq_list": line.strip(),
                    "trinity_fa": self.option('ref_fasta'),
                    "fq_type": fq_type,
                    "num": str(n)
                })
                self.bwa_tools.append(self.bwa)
        for j in range(len(self.bwa_tools)):
            self.bwa_tools[j].on('end', self.set_output, 'bwa')
        if self.bwa_tools:
            if len(self.bwa_tools) > 1:
                self.on_rely(self.bwa_tools, self.snp_run)
            elif len(self.bwa_tools) == 1:
                self.bwa_tools[0].on('end', self.snp_run)
        else:
            self.set_error("self.bwa_tools列表为空！", code = "22001202")
        for tool in self.bwa_tools:
            gevent.sleep(1)
            tool.run()

    def make_bamlist(self):
        bam_dict = dict()
        if self.in_bam:
            bamfile_list=sorted(os.listdir(self.in_bam))
            with open(self.work_dir + "/bamlist", "w") as fw:
                for i in bamfile_list:
                    fw.write(os.path.join(self.in_bam,i) + "\n")
        else:
            bamfile_list = glob2.glob(self.work_dir + "/**/ready_sort_bam/*")
            dir = os.path.dirname(bamfile_list[0])
            for i in bamfile_list:
                sample = i.split("/")[-1].split(".bam")[0]
                bam_dict[sample] = i
            df = pd.read_table(self.option('fq_list'), header=None)
            sample_list = df.iloc[:, 0].tolist()
            with open(self.work_dir + "/bamlist", "w") as fw:
                for i in sample_list:
                    fw.write(bam_dict[i] + "\n")

    def get_vcf_list(self):
        gevent.sleep(2)
        path = os.path.join(self.work_dir, 'snp_indel')
        files = os.listdir(path)
        vcf = os.path.join(self.output_dir, "vcf.list")
        with open(vcf, 'w') as w:
            for m in files:
                if "call" in m:
                    w.write("{}\n".format(os.path.join(path, m)))
            else:
                pass

    def bcftools_vcf_run(self):
        self.bcftools_vcf = self.add_tool("denovo_rna_v2.bcftool_vcf")
        self.get_vcf_list()
        self.bcftools_vcf.set_options({
            "vcf_list": os.path.join(self.output_dir, "vcf.list"),
            # "isoform_unigene": self.option('isoform_unigene')
        })
        self.bcftools_vcf.on('end', self.set_output, 'vcf_call')
        self.bcftools_vcf.on("end", self.vcf_filter_run)
        self.bcftools_vcf.run()

    def vcf_filter_run(self):
        self.vcf_filter=self.add_tool('ref_rna_v2.vcf_filter_samtools')
        # with open(self.bcftools_vcf.output_dir + "/pop.variant.vcf") as f, open(
        #         self.bcftools_vcf.output_dir + "/variant.vcf", "w") as w:
        #     for line in f:
        #         if not line.startswith("#CHROM"):
        #             w.write(line)
        #             continue
        #         new_line_1 = line.strip().split("\t")[0:9]
        #         new_line_2 = [x.split("/")[-1].split(".bam")[0] for x in line.strip().split("\t")[9:]]
        #         w.write("\t".join(new_line_1 + new_line_2) + "\n")
        options = {
            "input_file": self.bcftools_vcf.output_dir + "/pop.variant.vcf",
              }
        self.vcf_filter.set_options(options)
        self.vcf_filter.on("end",self.set_output,"vcf_filter")
        if self.inter_inbam:
            self.vcf_filter.on("end", self.snpfinal_run)
        else:
            self.vcf_filter.on("end",self.end)
        self.vcf_filter.run()

    def snp_run(self):
        self.split_ref_run()
        self.make_bamlist()
        files_list = glob2.glob(self.work_dir + "/split_file/*fasta*")
        n = 0
        for m in files_list:
            self.snp = self.add_tool("denovo_rna_v2.snp")
            self.snp.set_options({
                "trinity_fa": self.option('ref_fasta'),
                "bamlist": self.work_dir + "/bamlist",
                "bed_file": m,
                # "isoform_unigene": self.option('isoform_unigene'),
            })
            self.snp_tools.append(self.snp)
            n += 1
        for j in range(len(self.snp_tools)):
            self.snp_tools[j].on('end', self.set_output, 'snp_indel')
        if self.snp_tools:
            if len(self.snp_tools) > 1:
                self.on_rely(self.snp_tools, self.bcftools_vcf_run)
            elif len(self.snp_tools) == 1:
                self.snp_tools[0].on('end', self.bcftools_vcf_run)
        else:
            self.set_error("self.snp_tools列表为空！", code = "22001203")
        for tool in self.snp_tools:
            gevent.sleep(1)
            tool.run()

    def snpfinal_run(self):
        with open(self.work_dir + "/bamlist") as len_bam:
            sample_num = len(len_bam.readlines())
        self.snpfinal = self.add_tool('denovo_rna_v3.snpfinal_new2')
        self.snpfinal.set_options({
            "bamlist": sample_num,
            "call_vcf": self.vcf_filter.output_dir + "/final.vcf",
            'method':'samtools',
            'cds_bed': self.option("cds_bed"),
            'allt2g': self.option('allt2g'),
            'anno': self.option('anno')
            #"qual": self.option('qual'),
           # "dp": self.option('dp'),
        })
        # 这个地方需要在最后跑的那个tool那里加一个on('end', self.end)，这样module才会终止，并且这个
        # self.end和这最后一个tool的set_output是一起跑的，不会先end导致set_output失败
        self.snpfinal.on('end', self.end)
        self.snpfinal.on('end', self.set_output, 'snpfinal')
        self.snpfinal.run()

    def run(self):
        super(SnpModule, self).run()
        if self.option('in_bam'):
            self.in_bam=self.option('in_bam')
            self.inter_inbam = self.option('in_bam')
        else:
            self.in_bam=None
            self.inter_inbam = None
        if self.in_bam:
            self.samtools_index()
            self.snp_run()
        else:
            self.samtools_index()
            self.bwa_run(self.option("fq_type"), self.option("fq_list"))

    def samtools_index(self):
        self.samtools_path = "/bioinfo/align/samtools-1.6/samtools-1.6/samtools"
        if os.path.exists(self.option("ref_fasta").prop['path'] + ".fai"):
            os.remove(self.option("ref_fasta").prop['path'] + ".fai")
        cmd_index = "samtools faidx {}".format(self.option("ref_fasta").prop['path'])
        self.logger.info("开始构建参考序列索引")
        try:
            subprocess.check_call(cmd_index, shell=True)
            self.logger.info(cmd_index + " 运行完成！")
        except subprocess.CalledProcessError:
            print(traceback.format_exc())
            self.logger.info('CMD:{}'.format(cmd_index))
            self.set_error("建立索引失败", code = "22001204")

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
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
                os.link(oldfiles[i], newdir)

    def linkdirs(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(os.path.join(dirpath,"ready_sort_bam"))
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath,"ready_sort_bam",i) for i in allfiles]
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
                os.link(oldfiles[i], newdir)

    def set_output(self, event):
        '''
        设置结果目录
        '''
        obj = event['bind_object']
        self.logger.info('设置目录 {}'.format(event['data']))

        # 这个是创建一个名字叫event['data']，也就是你的tool的名字的文件夹，在这个module里面，这个obj默认的就是这个名字，所有都会产生在那个module下面
        # self.linkdir(obj.output_dir, event['data'])
        if event['data'] == 'snp_indel':
            self.linkdir(obj.output_dir, 'snp_indel')

        if event['data'] == 'bwa':
            self.linkdirs(obj.output_dir, 'bwa')
            # bam_dict=dict()
            # bamfile_list = glob2.glob(self.work_dir + "/**/ready_sort_bam/*")
            # dir = os.path.dirname(bamfile_list[0])
            # for i in bamfile_list:
            #     sample = i.split("/")[-1].split(".bam")[0]
            #     bam_dict[sample] = i
            # df = pd.read_table(self.option('fq_list'), header=None)
            # sample_list = df.iloc[:, 0].tolist()
            # with open(self.work_dir + "/bamlist", "w") as fw:
            #     for i in sample_list:
            #         os.link(bam_dict[i],os.path.join(self.workdir,"bwa",i+".bam"))
            #         fw.write(bam_dict[i] + "\n")
            self.in_bam=os.path.join(self.work_dir,"bwa")

        if event['data'] == "vcf_call":
            self.linkdir(obj.output_dir, 'vcf_call')

        if event['data'] == 'snpfinal':
            if os.path.exists(self.output_dir + "/pop.variant.vcf"):
                    os.remove(self.output_dir + "/pop.variant.vcf")
            os.link(self.bcftools_vcf.output_dir + "/pop.variant.vcf", self.output_dir + "/pop.variant.vcf")
            os.remove(self.output_dir + "/vcf.list")
            if os.path.getsize(obj.work_dir + "/snp_new") > 0:

                snp_detail = os.path.join(obj.output_dir, 'snp_detail')
                out1 = os.path.join(self.output_dir, 'snp_detail')
                if os.path.exists(self.output_dir + "/snp_detail"):
                    os.remove(self.output_dir + "/snp_detail")
                os.link(snp_detail, out1)

                snp_depth_statistics = os.path.join(obj.output_dir, 'snp_depth_statistics')
                out2 = os.path.join(self.output_dir, 'snp_depth_statistics')
                if os.path.exists(self.output_dir + "/snp_depth_statistics"):
                    os.remove(self.output_dir + "/snp_depth_statistics")
                os.link(snp_depth_statistics, out2)

                snp_homo_hete_statistics = os.path.join(obj.output_dir, 'snp_homo_hete_statistics')
                out3 = os.path.join(self.output_dir, 'snp_homo_hete_statistics')
                if os.path.exists(self.output_dir + "/snp_homo_hete_statistics"):
                    os.remove(self.output_dir + "/snp_homo_hete_statistics")
                os.link(snp_homo_hete_statistics, out3)

                snp_transition_tranversion_statistics = os.path.join(obj.output_dir, 'snp_transition_tranversion_statistics')
                out4 = os.path.join(self.output_dir, 'snp_transition_tranversion_statistics')
                if os.path.exists(self.output_dir + "/snp_transition_tranversion_statistics"):
                    os.remove(self.output_dir + "/snp_transition_tranversion_statistics")
                os.link(snp_transition_tranversion_statistics, out4)

                snp_cds_statistics = os.path.join(obj.output_dir,'snp_cds_statistics')
                out5 = os.path.join(self.output_dir, 'snp_cds_statistics')
                if os.path.exists(self.output_dir + "/snp_cds_statistics"):
                    os.remove(self.output_dir + "/snp_cds_statistics")
                os.link(snp_cds_statistics, out5)

                snp_anno_statistics = os.path.join(obj.output_dir, 'snp_anno_statistics')
                out6 = os.path.join(self.output_dir, 'snp_anno_statistics')
                if os.path.exists(self.output_dir + "/snp_anno_statistics"):
                    os.remove(self.output_dir + "/snp_anno_statistics")
                os.link(snp_anno_statistics, out6)

            if os.path.getsize(obj.work_dir + "/indel_new") > 0:
                indel_detail = os.path.join(obj.output_dir, 'indel_detail')
                out7 = os.path.join(self.output_dir, 'indel_detail')
                if os.path.exists(self.output_dir + "/indel_detail"):
                    os.remove(self.output_dir + "/indel_detail")
                os.link(indel_detail, out7)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "SNP分析结果目录"],
            ["snp_detail", " ", "SNP结果详情表"],
            ["indel_detail", " ", "InDel结果详情表"],
            ["snp_depth_statistics", " ", "SNP测序深度统计表"],
            ["snp_homo_hete_statistics", " ", "SNP类型统计表"],
            ["snp_transition_tranversion_statistics", " ", "SNP位点统计表"],
            ["filter_vcf", "", "SNP过滤结果vcf文件"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SnpModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo'
        data = {
            "id": "denovo_snp" + str(random.randint(1, 10000))+"yyyy",
            "type": "module",
            "name": "denovo_rna_v3.snp",
            "instant": False,
            "options": dict(
                ref_fasta="/mnt/ilustre/users/sanger-dev/workspace/20191014/Denovorna_tsg_35796/DenovoAssemble2Filter/output/Trinity.filter.unigene.fasta",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/workspace/20190416/Denovorna_tsg_33857/DenovoAssemble2/output/Trinity.filter.unigene.fasta",
                # fq_list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/denovo_snp/final_snp/test_files/test_bam/fq/fq_list",
                fq_list="/mnt/ilustre/users/sanger-dev/workspace/20191014/Denovorna_tsg_35796/FastpRna/output/fastq/fq_list.txt",
                fq_type="PE",
                # cds_bed="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/denovo_snp/final_snp/cds/all_predicted.bed",
                cds_bed="/mnt/ilustre/users/sanger-dev/workspace/20191014/Denovorna_tsg_35796/AnnotOrfpfam/output/all_predicted.bed",
                allt2g="/mnt/ilustre/users/sanger-dev/workspace/20191014/Denovorna_tsg_35796/AnnotOrfpfam/output/all_tran2gen.txt",
                # allt2g="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/denovo_snp/final_snp/cds/all_tran2gen.txt",
                # anno="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/denovo_snp/final_snp/Snpfinal_new/gene_anno_detail.xls",
                anno="/mnt/ilustre/users/sanger-dev/workspace/20191014/Denovorna_tsg_35796/AnnotClassBeta/output/all_annot.xls",
                # in_bam="/mnt/ilustre/users/sanger-dev/workspace/20190814/Single_denovo_snp8413yyyy/Snp/bwa"
                #qual="20",
                #dp="1"

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

