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


class DenovoAssembleModule(Module):
    """
    拼接以及新转录本预测
    version v1.0.1
    author: liubinxu
    last_modify: 2017.11.20
    """
    def __init__(self, work_id):
        super(DenovoAssembleModule, self).__init__(work_id)
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
            {"name": "assemble_method", "type": "string", "default": "trinity"},  # 选择拼接软件
            {"name": "min_contig_length", "type": "int", "default": 200},  # trinity报告出的最短的contig长度。默认为200
            {"name": "kmer_size", "type": "int", "default": 25},  # kmer 长度 <=32
            {"name": "min_kmer_cov", "type": "int", "default": 2},
            {"name": "jaccard_clip", "type": "bool", "default": False}, #分割高密度基因区间基因
            {"name": "no_normalize_reads", "type": "bool", "default": False}, #不做reads均一化
            {"name": "normalize_max_read_cov", "type": "int", "default": 50}, #reads均一化覆盖倍数
            {"name": "SS_lib_type", "type": "string", "default": 'none'},  # reads的方向，成对的reads: RF or FR; 不成对的reads: F or R，默认情况下，不设置此2参数
            ###过滤选项
            {"name": "transrate_filter", "type": "bool", "default": True}, #使用transrate过滤
            {"name": "species", "type": "string", "default": "All"},# 物种类别目前支持 All(表示Eukaryota) Fungi Animal Plant Protist 用于busco
            {"name": "cdhit_filter", "type": "bool", "default": True}, #使用cdhit过滤
            {"name": "cdhit_identity_threshold", "type": "float", "default": 0.95}, #cdhit过滤相似度
            {"name": "TPM_filter", "type": "bool", "default": True}, #使用表达量过滤
            {"name": "TPM_threshold", "type": "float", "default": 0.1}, #使用表达量过滤阈值

            ###统计选项
        ]
        self.add_option(options)
        self.trinity = self.add_tool('denovo_rna_v2.trinity')
        self.busco = self.add_tool('denovo_rna_v2.busco')
        self.cdhit = self.add_tool('denovo_rna_v2.cdhit_filter')
        self.transrate = self.add_tool('denovo_rna_v2.transrate')
        self.trinity_stat = self.add_tool('denovo_rna_v2.trinity_stat')
        self.exp_filter = self.add_tool('denovo_rna_v2.filterfasta')
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

        self.step.add_steps('trinity', 'busco', 'transrate', 'cdhit',
                            'trinity_stat', 'quntify', 'exp_filter')

    def check_options(self):
        """
        检查参数
        """
        if self.option('assemble_fa') and self.option('assemble_g2t'):
            self.logger.info('使用输入组装结果，跳过组装步骤')
        elif not self.option('sample_fq_list'):
            raise OptionError("样本列表文件依次为sample_id\tleft.fq\tright.fq", code="22000501")
        return True

    def merge_fastq(self):
        '''
        合并fastq数据2
        '''
        if self.option('filter_sample').prop.has_key('path'):
            with open(self.option('filter_sample').prop['path'], 'r') as f:
                lines = f.readlines()
                for line in lines:
                    self.filter_samples.append(line.strip())
        elif self.option('filter_sample_list'):
            self.filter_samples = self.option('filter_sample_list').split(',')
        else:
            pass
        with open(self.option("sample_fq_list")) as fq_list:
            lines = fq_list.readlines()
            for line in lines:
                line = line.strip().split('\t')
                if line[0] in self.filter_samples:
                    pass
                else:
                    self.right_fq.append(line[2])
                    self.left_fq.append(line[1])
            with open(self.merge_left, 'w') as left, open(self.merge_right, 'w') as right:
                for file in self.left_fq:
                    for fq in open(file, 'r'):
                        left.write(fq)
                for file in self.right_fq:
                    for fq in open(file, 'r'):
                        right.write(fq)

    def run_trinity(self):
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
            "fq_type": self.option("fq_type"),
            "fq_l": self.merge_left,
            "fq_r": self.merge_right,
            "min_contig_length": self.option("min_contig_length"),
            "kmer_size": self.option("kmer_size"),
            "min_kmer_cov": self.option("min_kmer_cov"),
            "jaccard_clip": self.option("jaccard_clip"),
            "no_normalize_reads": self.option("no_normalize_reads"),
            "normalize_max_read_cov": self.option("normalize_max_read_cov"),
            "SS_lib_type": trinity_strand
        }

        if self.option('assemble_fa') and self.option('assemble_g2t'):
            self.trinity.option("trinity_fa", self.option('assemble_fa'))
            self.trinity.option("gene2trans", self.option('assemble_g2t'))
            self.logger.info("跳过组装步骤，使用 {} 作为组装结果；使用 {} 作为基因转录本映射文件".format(self.option('assemble_fa'), self.option('assemble_g2t')))
            self.trinity.end()
        else:
            self.trinity.set_options(options)
            self.trinity.run()

    def run_busco(self):
        '''
        busco组装结果评估
        '''
        options = {
            "mode": "tran",
            "species": self.option("species"),
            "fa": self.trinity.option("trinity_fa").prop['path']
        }
        self.busco.set_options(options)
        self.busco.run()
        #self.busco_result = self.busco.option("busco_result")

    def run_transrate(self):
        '''
        transrate组装结果评估
        '''
        options = {
            "left": self.merge_left,
            "right": self.merge_right,
            "assembly": self.trinity.option("trinity_fa").prop['path']
        }
        self.transrate.set_options(options)
        #self.transrate_result = self.transrate.option('result')
        #self.exp_matrix = self.transrate.option('exp_result')
        if self.option('transrate_filter'):
            self.transrate_filter_done = True
            #self.filtered_fasta = self.transrate.option('good_fa')
        else:
            self.logger.info('不使用transrate过滤')

        self.transrate.run()

    def run_trinity_stat(self):
        '''
        组装结果统计
        '''
        #self.option('gene2trans')
        self.logger.info('** gene2trans结果为 {}'.format(self.trinity.option("gene2trans")))
        self.logger.info('** fasta结果为 {}'.format(self.trinity.option("trinity_fa").prop['path']))
        self.logger.info('** busco结果为 {}'.format(self.busco.option("busco_result")))
        self.logger.info('** transrate结果为 {}'.format(self.transrate.option("result")))
        options = {
            "fasta": self.trinity.option("trinity_fa").prop['path'],
            "gene2trans": self.trinity.option("gene2trans").prop['path'],
            "busco_result": self.busco.option("busco_result").prop['path'],
            "transrate_result": self.transrate.option('result').prop['path'],
            "exp_result": self.transrate.option('exp_result').prop['path'],
            "marker": True,
        }
        self.trinity_stat.set_options(options)
        self.trinity_stat.run()



    def run_exp_filter(self):
        '''
        使用TPM过滤低表达序列
        '''
        options = {
            "transcripts": self.get_unfilter_fa(),
            "matrix": self.transrate.option('exp_result').prop['path'],
            "matrix_type": "transrate",
            "min_expr_any": self.option('TPM_threshold')
        }
        if self.option('TPM_filter'):
            self.TPM_filter_done = True
            self.exp_filter.set_options(options)
            self.exp_filter.run()
            #self.filtered_fasta = self.exp_filter.option('exp_filterfasta')
        else:
            self.logger.info('不进行TPM表达量过滤，此步跳过')
            self.exp_filter.end()

    def run_cdhit_filter(self):
        '''
        使用cdhit过滤冗余序列
        '''
        options = {
            "query": self.get_unfilter_fa(),
            "identity": self.option("cdhit_identity_threshold"),
        }
        self.cdhit.set_options(options)
        if self.option('cdhit_filter'):
            self.cdhit.run()
            #self.filtered_fasta = self.cdhit.option('cdhit_filter_fasta').prop['path']
        else:
            self.logger.info('不进行cdhit冗余序列过滤，此步跳过')
            self.cdhit.end()

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def get_filter_fa(self):
        '''
        获取过滤得到的fasta文件
        '''
        filter_fa = self.trinity.option("trinity_fa").prop['path']
        if self.option('cdhit_filter'):
            filter_fa = self.cdhit.option('cdhit_filter_fasta').prop['path']
        elif self.option('TPM_filter'):
            filter_fa = self.exp_filter.option('exp_filterfasta').prop['path']
        elif self.option('transrate_filter'):
            filter_fa = self.transrate.option('good_fa').prop['path']
        else:
            pass

        return filter_fa

    def get_unfilter_fa(self):
        '''
        获取需要过滤的fasta文件
        '''
        un_filter_fa = self.trinity.option("trinity_fa").prop['path']
        if self.TPM_filter_done:
            un_filter_fa = self.exp_filter.option('exp_filterfasta').prop['path']
        elif self.transrate_filter_done:
            un_filter_fa = self.transrate.option('good_fa').prop['path']
        else:
            pass
        return un_filter_fa

    def run(self):
        '''
        设置运行顺序
        '''
        super(DenovoAssembleModule, self).run()
        if self.option('assemble_fa') and self.option('assemble_g2t'):
            self.on_rely([self.busco, self.transrate], self.run_trinity_stat)
        else:
            self.trinity.on('end', self.run_busco)
            self.trinity.on('end', self.run_transrate)
            self.trinity.on('end', self.set_output, 'trinity')
            self.on_rely([self.trinity, self.busco, self.transrate], self.run_trinity_stat)

        if self.option('TPM_filter'):
            self.on_rely([self.trinity_stat], self.run_exp_filter)
            self.exp_filter.on('end', self.set_output, 'exp_filter')
            if self.option('cdhit_filter'):
                self.on_rely([self.exp_filter], self.run_cdhit_filter)
                self.cdhit.on('end', self.set_output, 'cdhit')
            else:
                pass
        else:
            if self.option('cdhit_filter'):
                self.on_rely([self.trinity_stat], self.run_cdhit_filter)
                self.cdhit.on('end', self.set_output, 'cdhit')
            else:
                pass

        self.busco.on('end', self.set_output, 'busco')
        self.transrate.on('end', self.set_output, 'transrate')
        self.trinity_stat.on('end', self.set_output, 'trinity_stat')

        '''
        if self.option("TPM_filter"):
            if self.option("transrate_filter"):
                self.transrate.on('end',self.run_exp_filter)
            else:
                pass
        else:
            pass

        if self.option("cdhit_filter"):
            if self.option("TPM_filter"):
                self.exp_filter.on('end',self.run_cdhit_filter)
            else:
                if self.option("transrate_filter"):
                    self.transrate.on('end',self.run_cdhit_filter)
                else:
                    self.trinity.on('end',self.run_cdhit_filter)
        else:
        '''
        self.check_options()
        self.merge_fastq()

        if self.option('assemble_fa') and self.option('assemble_g2t'):
            self.trinity.option("trinity_fa", self.option('assemble_fa'))
            self.trinity.option("gene2trans", self.option('assemble_g2t'))
            self.logger.info("跳过组装步骤，使用 {} 作为组装结果；使用 {} 作为基因转录本映射文件".format(
                self.option('assemble_fa'), self.option('assemble_g2t')))
            self.run_busco()
            self.run_transrate()
        else:

            self.run_trinity()


    def set_output(self, event):
        '''
        设置结果目录
        '''
        obj = event['bind_object']
        self.logger.info('设置目录 {}'.format(event['data']))

        self.linkdir(obj.output_dir, event['data'])

        if event['data'] == 'trinity':
            trinity_fa = os.path.join(obj.output_dir, 'Trinity.fasta')
            out = os.path.join(self.output_dir, 'Trinity.fasta')
            os.link(trinity_fa, out)
            self.assemble_fasta = out
        if event['data'] == 'trinity_stat':
            if (not self.option('cdhit_filter')) and (not self.option('TPM_filter') and (not self.option('transrate_filter'))):
                self.end()
        if event['data'] == 'transrate':
            if (not self.option('cdhit_filter')) and (not self.option('TPM_filter')) and self.option('transrate_filter'):
                filter_fa = self.get_filter_fa()
                out = os.path.join(self.output_dir, 'Trinity.filter.fasta')
                os.link(filter_fa, out)
                # self.end()
        if event['data'] == 'exp_filter':
            if (not self.option('cdhit_filter') and self.option('TPM_filter')):
                filter_fa = self.get_filter_fa()
                out = os.path.join(self.output_dir, 'Trinity.filter.fasta')
                os.link(filter_fa, out)
                self.end()
        if event['data'] == 'cdhit':
            filter_fa = self.get_filter_fa()
            out = os.path.join(self.output_dir, 'Trinity.filter.fasta')
            os.link(filter_fa, out)
            self.end()
        '''
        if event['data'] == 'trinity_filter_stat':
            filter_fa = self.get_filter_fa()
            #out = os.path.join(self.output_dir, 'Trinity.filter.fasta')
            os.link(filter_fa, out)
            self.end()
        '''

    def linkdir(self, olddir, newname, mode='link'):
        """
        移动tools目录下的输出文件/文件夹到modules输出文件夹下
        """
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code="22000501")
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
        self.option("filter_fa",  os.path.join(self.output_dir, 'Trinity.filter.fasta'))
        self.option("filter_fa").set_gene2tran(self.trinity.option("gene2trans").prop['path'])
        self.option("filter_fa").get_unigene(os.path.join(self.output_dir, 'Trinity.filter.unigene.fasta'), os.path.join(self.output_dir, 'Trinity.filter'))
        self.option("unigene_filter_fa",  os.path.join(self.output_dir, 'Trinity.filter.unigene.fasta'))
        '''
        result_dir.add_relpath_rules([
            [r".", "", "结果输出目录"],
            ["Stringtie", "", "拼接后的各样本文件夹"],
            ["StringtieMerge", "", "拼接组装合并之后结果文件夹"],
            ["StringtieMerge/merged.gtf", "gtf", "样本合并之后的注释文件"],

        ])
        result_dir.add_regexp_rules([
            [r"Stringtie/.*_out\.gtf$", "gtf", "样本拼接之后的注释文件"],
            [r"Stringtie/.*_out\.fa$", "fasta", "样本拼接之后序列文件"],
            [r"Statistics/trans_count_stat_.*\.txt$", "txt", "新转录本步长统计文件"],
            [r"Statistics/old_.*\.txt$", "txt", "统计结果文件"],
            [r"Statistics/new_.*\.txt$", "txt", "统计结果文件"],

        ])
        '''
        super(DenovoAssembleModule, self).end()

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
            "id": "denovo_assemble" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "denovo_rna_v2.denovo_assemble",
            "instant": False,
            "options": dict(
                sample_fq_list = test_dir + "/test_data1/fq.list"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
