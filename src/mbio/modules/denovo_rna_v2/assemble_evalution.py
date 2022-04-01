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


class AssembleEvalutionModule(Module):
    """
    拼接以及新转录本预测
    version v1.0.1
    author: liubinxu
    last_modify: 2017.12.15
    """
    def __init__(self, work_id):
        super(AssembleEvalutionModule, self).__init__(work_id)
        options = [
            ###组装选项
            {"name": "sample_fq_list", "type": "string"}, #样本序列信息表
            {"name": "assemble_fa", "type": "infile", "format": "denovo_rna_v2.trinity_fasta"},  # 组装结果
            {"name": "fq_type", "type": "string", "default": "PE"}, #样本序列类型PE或SE
            {"name": "filter_sample", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "min_contig_length", "type": "int", "default": 200},
            {"name": "assemble_g2t", "type": "string", "default": "none"},  # 组装结果的转录本基因对应关系文件
            {"name": "species", "type": "string", "default": "All"}, # 物种类别目前支持 All(表示Eukaryota) Fungi Animal Plant Protist 用于busco
        ]
        self.add_option(options)
        #self.trinity = self.add_tool('denovo_rna_v2.trinity')
        self.fastq_merge = self.add_tool('denovo_rna_v2.fastq_merge')
        self.busco = self.add_tool('denovo_rna_v2.busco')
        self.transrate = self.add_tool('denovo_rna_v2.transrate')
        self.trinity_stat = self.add_tool('denovo_rna_v2.trinity_stat')
        self.step.add_steps('busco', 'transrate', 'trinity_stat')

        self.right_fq = []
        self.left_fq = []
        self.merge_left = self.work_dir +  '/merge.left.fq'
        self.merge_right = self.work_dir +  '/merge.right.fq'

    def check_options(self):
        """
        检查参数
        """
        if not (self.option('sample_fq_list') and self.option('assemble_fa') and self.option('assemble_g2t')):
            raise OptionError("参数为必须参数", code = "22000101")
        if self.option('assemble_g2t') == "none":
            self.option("assemble_fa").get_gene2tran(self.work_dir + "/gene2tran.txt")
            self.option('assemble_g2t', self.work_dir + '/gene2tran.txt')
        if self.option('fq_type') == "PE":
            with open(self.option('sample_fq_list'), 'r') as f:
                if len(f.readline().strip().split("\t")) < 3:
                    raise OptionError("样本文件格式错误 样本列表文件依次为sample_id\tleft.fq\tright.fq", code = "22000102")
        elif self.option('fq_type') == "SE":
            with open(self.option('sample_fq_list'), 'r') as f:
                if len(f.readline().strip().split("\t")) < 2:
                    raise OptionError("样本文件格式错误 样本列表文件依次为sample_id\tleft.fq", code = "22000103")
        return True

    def run_merge_fastq(self):
        '''
        合并fastq数据
        '''
        options = {
            "sample_fq_list": self.option("sample_fq_list"),
            "filter_sample": self.option("filter_sample"),
        }
        self.fastq_merge.set_options(options)
        self.fastq_merge.run()

    """
    def merge_fastq(self):
        '''
        合并fastq数据
        '''
        with open(self.option("sample_fq_list")) as fq_list:
            lines = fq_list.readlines()
            for line in lines:
                line = line.strip().split('\t')
                self.right_fq.append(line[2])
                self.left_fq.append(line[1])
            with open(self.merge_left, 'w') as left, open(self.merge_right, 'w') as right:
                for file in self.left_fq:
                    for fq in open(file, 'r'):
                        left.write(fq)
                for file in self.right_fq:
                    for fq in open(file, 'r'):
                        right.write(fq)
    """

    def run_busco(self):
        '''
        busco组装结果评估
        '''
        options = {
            "mode": "tran",
            "species": self.option("species"),
            # "g2t": self.option("assemble_g2t"),
            "fa": self.option("assemble_fa")
        }
        self.busco.set_options(options)
        self.busco.run()
        #self.busco_result = self.busco.option("busco_result")

    def run_transrate(self):
        '''
        transrate组装结果评估
        '''
        options = {
            "left": self.fastq_merge.option('merge_left').prop['path'],
            "assembly": self.option("assemble_fa")
        }
        if self.option("fq_type") == "PE":
            options.update({
                "right": self.fastq_merge.option('merge_right').prop['path']
            })
        self.transrate.set_options(options)
        self.transrate.run()

    def run_trinity_stat(self):
        '''
        组装结果统计
        '''
        #self.option('gene2trans')
        self.logger.info('** fasta结果为 {}'.format(self.option("assemble_fa")))
        self.logger.info('** gene2trans结果为 {}'.format(self.option("assemble_g2t")))
        self.logger.info('** busco结果为 {}'.format(self.busco.option("busco_result")))
        self.logger.info('** transrate结果为 {}'.format(self.transrate.option("result")))
        options = {
            "fasta": self.option("assemble_fa"),
            "gene2trans": self.option("assemble_g2t"),
            "busco_result": self.busco.option("busco_result").prop['path'],
            "min_len": self.option("min_contig_length"),
            "transrate_result": self.transrate.option('result').prop['path'],
            "exp_result": self.transrate.option('exp_result').prop['path'],
        }
        self.trinity_stat.set_options(options)
        self.trinity_stat.run()

    def run(self):
        '''
        设置运行顺序
        '''
        super(AssembleEvalutionModule , self).run()

        self.fastq_merge.on('end', self.run_transrate)
        self.on_rely([self.busco, self.transrate], self.run_trinity_stat)

        #设置结果目录
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
        self.run_merge_fastq()
        self.run_busco()


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
            if os.path.exists(self.output_dir + "/Trinity.fasta"):
                os.remove(self.output_dir + "/Trinity.fasta")
                os.link(trinity_fa, out)
            else:
                os.link(trinity_fa, out)
            self.assemble_fasta = out
        if event['data'] == 'trinity_stat':
            self.end()
            #os.link()


    def linkdir(self, olddir, newname, mode='link'):
        """
        移动tools目录下的输出文件/文件夹到modules输出文件夹下
        """
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code = "22000104")
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
        super(AssembleEvalutionModule, self).end()

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
            "id": "assemble_evalution" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "denovo_rna_v2.assemble_evalution",
            "instant": False,
            "options": dict(
                sample_fq_list = test_dir + "/test_data1/fq.list",
                assemble_fa = test_dir + "/test_data1/Trinity.filter.fasta",
                assemble_g2t = test_dir + "/test_data1/Trinity.fasta.gene_trans_map"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
