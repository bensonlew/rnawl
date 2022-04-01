# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.ref_rna.trans_step import *
import re
from mbio.files.sequence.file_sample import FileSampleFile
import unittest


class DenovoAssemble3Module(Module):
    """
    拼接以及新转录本预测
    version v1.0.1
    author: liubinxu
    last_modify: 2019.05.29
    """
    def __init__(self, work_id):
        super(DenovoAssemble3Module, self).__init__(work_id)
        options = [
            ###组装选项
            {"name": "sample_fq_list", "type": "string"}, #样本序列信息表
            {"name": "fq_type", "type": "string", "default": "PE"}, #样本序列类型PE或SE
            {"name": "strand_direct", "type": "string", "default": "none"},  # 链特异性时选择正负链

            {"name": "assemble_fa", "type": "string", "default": ""},  # 组装结果
            {"name": "assemble_g2t", "type": "string", "default": ""},  # 组装结果的转录本基因对应关系文件
            {"name": "filter_sample", "type": "infile", "format": "denovo_rna_v2.common"},  # 删除过滤样本 sample_file
            {"name": "filter_sample_list", "type": "string", "default": ""},  # 删除过滤样本 "sample1,sample2..."
            {'name': 'filter_fa', 'type': 'outfile', "format": "denovo_rna_v2.trinity_fasta"},
            {'name': 'unigene_filter_fa', 'type': 'outfile', "format": "denovo_rna_v2.trinity_fasta"},
            {"name": "assemble_soft", "type": "string", "default": "trinity"},  # 选择拼接软件trinity, spades
            {"name": "assemble_method", "type": "string", "default": "total"}, # 选择拼接软件sample, group, total
            {"name": "group", "type": "infile", "format": "denovo_rna_v2.group_table"},

            # trinity 运行参数
            {"name": "min_contig_length", "type": "int", "default": 200},  # trinity报告出的最短的contig长度。默认为200
            {"name": "kmer_size", "type": "int", "default": 25},  # kmer 长度 <=32
            {"name": "lines", "type": "int", "default": 1000},  # 分布式命令每个文件行数
            {"name": "min_kmer_cov", "type": "int", "default": 5},
            {"name": "jaccard_clip", "type": "bool", "default": False}, #分割高密度基因区间基因
            {"name": "no_normalize_reads", "type": "bool", "default": False}, #不做reads均一化
            {"name": "normalize_max_read_cov", "type": "int", "default": 50}, #reads均一化覆盖倍数
            {"name": "SS_lib_type", "type": "string", "default": 'none'},  # reads的方向，成对的reads: RF or FR; 不成对的reads: F or R，默认情况下，不设置此2参数
            # spades 运行参数
            {'type': 'string', 'name': 'spades_k', 'default': "auto"},
            {'type': 'string', 'name': 'spades_c', 'default': "auto"},

            {'type': 'float', 'name': 'tgicl_p', 'default': 88},
            {'type': 'int', 'name': 'tgicl_c', 'default': 10},


            {"name": "cdhit_filter", "type": "bool", "default": True}, #使用cdhit过滤
            {"name": "cdhit_identity_threshold", "type": "float", "default": 0.95}, #cdhit过滤相似度
            {"name": "TPM_filter", "type": "bool", "default": True}, #使用表达量过滤
            {"name": "TPM_threshold", "type": "float", "default": 0.1}, #使用表达量过滤阈值
            # count 过滤阈值
            {'type': 'int', 'name': 'filter_200', 'default': 0},
            {'type': 'int', 'name': 'filter_500', 'default': 0},
            {'type': 'int', 'name': 'filter_1000', 'default': 0},

            ###统计选项
        ]
        self.add_option(options)
        self.fastq_merge = self.add_tool('denovo_rna_v2.fastq_merge2')
        '''
        self.trinity = self.add_module('denovo_rna_v2.denovo_trinity')
        self.spades = self.add_tool('denovo_rna_v2.spades')
        '''

        self.tgicl = self.add_tool('denovo_rna_v2.tgicl')
        '''
        self.busco = self.add_tool('denovo_rna_v2.busco')
        self.cdhit = self.add_tool('denovo_rna_v2.cdhit_filter')
        self.transrate = self.add_tool('denovo_rna_v2.transrate')
        self.trinity_stat = self.add_tool('denovo_rna_v2.trinity_stat')
        self.exp_filter = self.add_tool('denovo_rna_v2.filterfasta')
        '''
        self.normal_tools = self.add_tool('denovo_rna_v2.normalize_fastq')
        #self.filter_busco = self.add_tool('denovo_rna_v2.busco')
        #self.filter_transrate = self.add_tool('denovo_rna_v2.transrate')
        #self.filter_trinity_stat = self.add_tool('denovo_rna_v2.trinity_stat')
        self.filter_samples = []
        self.right_fq = []
        self.left_fq = []
        self.merge_left = self.work_dir +  '/merge.left.fq'
        self.merge_right = self.work_dir +  '/merge.right.fq'
        self.exp_matrix = ''
        self.assemble_fasta = ''
        self.filtered_fasta = ''
        self.gene2trans = ''

        self.busco_result = ''
        self.transrate_result = ''
        self.filter_busco_result = ''
        self.filter_transrate_result = ''
        self.filter_exp_matrix = ''

        self.transrate_filter_done = False
        self.TPM_filter_done = False

        self.step.add_steps('fastqmerge','trinity', 'busco', 'transrate', 'cdhit',
                            'trinity_stat', 'quntify', 'exp_filter')

    def check_options(self):
        """
        检查参数
        """
        if not self.option('sample_fq_list'):
            raise OptionError("样本列表文件依次为sample_id\tleft.fq\tright.fq", code = "22001701")

        if self.option('fq_type') == "PE":
            with open(self.option('sample_fq_list'), 'r') as f:
                if len(f.readline().strip().split("\t")) < 3:
                    raise OptionError("样本文件格式错误 样本列表文件依次为sample_id\tleft.fq\tright.fq", code = "22001702")
        elif self.option('fq_type') == "SE":
            with open(self.option('sample_fq_list'), 'r') as f:
                if len(f.readline().strip().split("\t")) < 2:
                    raise OptionError("样本文件格式错误 样本列表文件依次为sample_id\tleft.fq", code = "22001703")
        return True

    def run(self):
        '''
        设置运行顺序
        '''
        super(DenovoAssemble3Module, self).run()
        self.run_merge_fastq()

    def run_merge_fastq(self):
        '''
        合并fastq数据
        '''
        options = {
            "sample_fq_list": self.option("sample_fq_list"),
            "filter_sample": self.option("filter_sample"),
            "group": self.option("group"),
            "assemble_method": self.option("assemble_method")
        }
        self.fastq_merge.set_options(options)
        if self.option("assemble_soft") == "trinity":
            self.fastq_merge.on('end', self.run_assemble)
        else:
            self.fastq_merge.on('end', self.run_normalization_fastq)
        self.fastq_merge.run()

    def run_normalization_fastq(self):
        '''
        read insilico normalization
        '''
        options = {
            "l": self.fastq_merge.output_dir + "/merge.left.fq",
            "r": self.fastq_merge.output_dir + "/merge.right.fq",
            "soft": "Trinity",
            "kmer": 25,
        }
        self.normal_tools.set_options(options)
        self.normal_tools.on('end', self.run_assemble)
        self.normal_tools.run()

    def run_assemble(self):
        if self.option("assemble_method") == "sample":
            self.run_assemble_sample()
        elif self.option("assemble_method") == "group":
            self.run_assemble_group()
        elif self.option("assemble_method") == "total":
            self.run_assemble_total()

    def run_assemble_total(self):
        if self.option("assemble_soft") == "trinity":
            ass_tool = self.add_module('denovo_rna_v2.denovo_trinity')
            ass_tool = self.trinity_tool_init(ass_tool, "TRINITY",
                                              self.fastq_merge.output_dir + "/merge.left.fq",
                                              self.fastq_merge.output_dir + "/merge.right.fq")

        elif self.option("assemble_soft") == "spades":
            
            ass_tool = self.add_tool('denovo_rna_v2.spades')
            ass_tool = self.spades_tool_init(ass_tool, "NODE",
                                             self.normal_tools.output_dir + "/normalize_1.fq",
                                             self.normal_tools.output_dir + "/normalize_2.fq")

        ass_tool.on('end', self.set_output, 'total')
        ass_tool.run()

    def run_assemble_group(self):
        self.assemble_tools = dict()
        group_spname =  self.option("group").get_group_spname()
        for group in group_spname.keys():
            if self.option("assemble_soft") == "trinity":
                ass_tool = self.add_module('denovo_rna_v2.denovo_trinity')
                ass_tool = self.trinity_tool_init(ass_tool, group,
                                             self.fastq_merge.output_dir + "/{}.merge.left.fq".format(group),
                                             self.fastq_merge.output_dir + "/{}.merge.right.fq".format(group))
                self.assemble_tools[group] = ass_tool
            elif self.option("assemble_soft") == "spades":
                ass_tool = self.add_tool('denovo_rna_v2.spades')
                ass_tool = self.spades_tool_init(ass_tool, group,
                                             self.fastq_merge.output_dir + "/{}.merge.left.fq".format(group),
                                             self.fastq_merge.output_dir + "/{}.merge.right.fq".format(group))
                self.assemble_tools[group] = ass_tool
        self.on_rely(self.assemble_tools.values(), self.run_tgicl)
        for group, tool in self.assemble_tools.items():
            tool.run()

    def run_assemble_sample(self):
        self.assemble_tools = dict()
        with open(self.option("sample_fq_list")) as fq_list:
            lines = fq_list.readlines()
            for line in lines:
                line = line.strip().split('\t')
                sample = line[0]
                if line[0] in self.filter_samples:
                    pass
                else:
                    if self.option("assemble_soft") == "trinity":
                        ass_tool = self.add_module('denovo_rna_v2.denovo_trinity')
                        ass_tool = self.trinity_tool_init(ass_tool, sample,
                                                          self.fastq_merge.output_dir + "/{}.merge.left.fq".format(sample),
                                                          self.fastq_merge.output_dir + "/{}.merge.right.fq".format(sample))
                        self.assemble_tools[sample] = ass_tool
                    elif self.option("assemble_soft") == "spades":
                        ass_tool = self.add_tool('denovo_rna_v2.spades')
                        ass_tool = self.spades_tool_init(ass_tool, sample,
                                                         self.fastq_merge.output_dir + "/{}.merge.left.fq".format(sample),
                                                         self.fastq_merge.output_dir + "/{}.merge.right.fq".format(sample))
                        self.assemble_tools[sample] = ass_tool
        self.on_rely(self.assemble_tools.values(), self.run_tgicl)
        for sample, tool in self.assemble_tools.items():
            tool.run()


    def run_tgicl(self):
        trans_list = list()
        if self.option("assemble_soft") == "trinity":
            for sample, tool in self.assemble_tools.items():
                trans_list.append(tool.output_dir + '/Trinity.fasta')
        elif self.option("assemble_soft") == "spades":
            for sample, tool in self.assemble_tools.items():
                trans_list.append(tool.output_dir + '/assemble.fa')
        options = {
            "transcript": ",".join(trans_list),
            "p": self.option("tgicl_p"),
            "c": self.option("tgicl_c"),
        }
        self.tgicl.set_options(options)
        self.tgicl.on('end', self.set_output, 'tgicl')
        self.tgicl.run()



    def spades_tool_init(self, ass_tool, name, fq_l, fq_r):
        '''
        运行组装结果
        '''
        trinity_strand = 'none'

        options = {
            "t": 20,
            "k": self.option("spades_k"),
            "c": self.option("spades_c"),
            "n": name
        }
        if self.option("fq_type") == "PE":
            options.update({
                "l": fq_l,
                "r": fq_r,
            })
        elif self.option("fq_type") == "SE":
            options.update({
                "s": fq_l,
            })
        else:
            pass
        ass_tool.set_options(options)
        return ass_tool


    def trinity_tool_init(self, ass_tool, name, fq_l, fq_r):
        '''
        运行组装结果
        '''
        trinity_strand = 'none'
        if self.option('strand_direct') == 'forward':
            if self.option("fq_type") == "PE":
                trinity_strand = "RF"
            else:
                trinity_strand = "R"
        elif self.option('strand_direct') == 'reverse':
            if self.option("fq_type") == "PE":
                trinity_strand = "FR"
            else:
                trinity_strand = "F"
        else:
            pass

        options = {
            "name": name,
            "fq_type": self.option("fq_type"),
            "lines": self.option("lines"),
            "min_contig_length": self.option("min_contig_length"),
            "kmer_size": self.option("kmer_size"),
            "min_kmer_cov": self.option("min_kmer_cov"),
            "jaccard_clip": self.option("jaccard_clip"),
            "no_normalize_reads": self.option("no_normalize_reads"),
            "normalize_max_read_cov": self.option("normalize_max_read_cov"),
            "SS_lib_type": trinity_strand
        }
        if self.option("fq_type") == "PE":
            options.update({
                "left": fq_l,
                "right": fq_r,
            })
        elif self.option("fq_type") == "SE":
            options.update({
                "single": fq_l,
            })
        else:
            pass
        ass_tool.set_options(options)
        return ass_tool


    def set_output(self, event):
        '''
        设置结果目录
        '''
        obj = event['bind_object']
        self.logger.info('设置目录 {}'.format(event['data']))

        # self.linkdir(obj.output_dir, event['data'])

        if event['data'] == 'total':
            if self.option("assemble_soft") == "trinity":
                trinity_fa = os.path.join(obj.output_dir, 'Trinity.fasta')
                out = os.path.join(self.output_dir, 'assemble_raw.fasta')
                gene2trans = os.path.join(obj.output_dir, 'Trinity.fasta.gene_trans_map')
                out1 = os.path.join(self.output_dir, 'assemble_raw.gene_trans_map')
                if os.path.exists(out):
                    os.remove(out)
                os.link(trinity_fa, out)
                if os.path.exists(out1):
                    os.remove(out1)
                os.link(gene2trans, out1)
                self.assemble_fasta = out
            elif self.option("assemble_soft") == "spades":
                trinity_fa = os.path.join(obj.output_dir, 'assemble.fa')
                out = os.path.join(self.output_dir, 'assemble_raw.fasta')
                gene2trans = os.path.join(obj.output_dir, 'assemble.g2t')
                out1 = os.path.join(self.output_dir, 'assemble_raw.gene_trans_map')
                if os.path.exists(out):
                    os.remove(out)
                os.link(trinity_fa, out)
                if os.path.exists(out1):
                    os.remove(out1)
                os.link(gene2trans, out1)
                self.assemble_fasta = out
            self.end()

        if event['data'] == 'tgicl':
            trinity_fa = os.path.join(obj.output_dir, 'merge.fa')
            out = os.path.join(self.output_dir, 'assemble_raw.fasta')
            gene2trans = os.path.join(obj.output_dir, 'merge.g2t')
            out1 = os.path.join(self.output_dir, 'assemble_raw.gene_trans_map')
            if os.path.exists(out):
                os.remove(out)
            os.link(trinity_fa, out)
            if os.path.exists(out1):
                os.remove(out1)
            os.link(gene2trans, out1)
            self.assemble_fasta = out
            self.end()


    def linkdir(self, olddir, newname, mode='link'):
        """
        移动tools目录下的输出文件/文件夹到modules输出文件夹下
        """
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code = "22001701")
        newdir = os.path.join(self.output_dir, newname)
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        os.mkdir(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            else:
                new1 = os.path.join(newdir, os.path.basename(oldfiles[i]))
                os.system("mv {} {}".format(oldfiles[i], new1))

    def end(self):
        #result_dir = self.add_upload_dir(self.output_dir)
        '''
        self.option("filter_fa",  os.path.join(self.output_dir, 'Trinity.filter.fasta'))
        self.option("filter_fa").set_gene2tran(self.trinity_merge.option("gene2trans").prop['path'])
        self.option("filter_fa").get_unigene(os.path.join(self.output_dir, 'Trinity.filter.unigene.fasta'), os.path.join(self.output_dir, 'Trinity.filter'))
        self.option("unigene_filter_fa",  os.path.join(self.output_dir, 'Trinity.filter.unigene.fasta'))
        '''

        super(DenovoAssemble3Module, self).end()

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
            "id": "denovo_assemble3" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "denovo_rna_v2.denovo_assemble3",
            "instant": False,
            "options": dict(
                sample_fq_list = "/mnt/ilustre/users/sanger-dev/workspace/20190617/Denovorna_tsg_34448/HiseqQc/output/sickle_dir/fq_list.txt",
                group = "/mnt/ilustre/users/sanger-dev/workspace/20190617/Denovorna_tsg_34448/remote_input/group/default_group.txt",
                assemble_soft = "trinity",
                assemble_method = "total",
                lines = 10,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
