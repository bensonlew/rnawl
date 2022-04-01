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
from Bio import SeqIO


class DenovoTrinityModule(Module):
    """
    拼接以及新转录本预测
    version v1.0.1
    author: liubinxu
    last_modify: 2019.06.17
    """
    def __init__(self, work_id):
        super(DenovoTrinityModule, self).__init__(work_id)
        options = [
            ###组装选项
            {"name": "sample_fq_list", "type": "string"}, #样本序列信息表
            {"name": "fq_type", "type": "string", "default": "PE"}, #样本序列类型PE或SE
            {"name": "strand_direct", "type": "string", "default": "none"},  # 链特异性时选择正负链

            {"name": "assemble_fa", "type": "string", "default": ""},  # 组装结果
            {"name": "assemble_g2t", "type": "string", "default": ""},  # 组装结果的转录本基因对应关系文件

            {"name": "left", "type": "infile", "format": "sequence.fastq"},
            {"name": "right", "type": "infile", "format": "sequence.fastq"},
            {"name": "single", "type": "infile", "format": "sequence.fastq"},
            {'name': 'filter_fa', 'type': 'outfile', "format": "denovo_rna_v2.trinity_fasta"},
            {'name': 'unigene_filter_fa', 'type': 'outfile', "format": "denovo_rna_v2.trinity_fasta"},
            {"name": "assemble_method", "type": "string", "default": "trinity"},  # 选择拼接软件
            {"name": "min_contig_length", "type": "int", "default": 200},  # trinity报告出的最短的contig长度。默认为200
            {"name": "name", "type": "string", "default": "TRNITY"},

            {"name": "kmer_size", "type": "int", "default": 25},  # kmer 长度 <=32
            {"name": "lines", "type": "int", "default": 1000},  # 分布式命令每个文件行数
            {"name": "min_kmer_cov", "type": "int", "default": 5},
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
            # count 过滤阈值
            {'type': 'int', 'name': 'filter_200', 'default': 0},
            {'type': 'int', 'name': 'filter_500', 'default': 0},
            {'type': 'int', 'name': 'filter_1000', 'default': 0},

            ###统计选项
        ]
        self.add_option(options)
        # self.fastq_merge = self.add_tool('denovo_rna_v2.fastq_merge')
        self.trinity = self.add_tool('denovo_rna_v2.trinity2')
        self.trinity_dis = []
        self.trinity_merge = self.add_tool('denovo_rna_v2.trinity2_merge')

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
        pass
        '''

        if self.option('fq_type') == "PE":
            with open(self.option('sample_fq_list'), 'r') as f:
                if len(f.readline().strip().split("\t")) < 3:
                    raise OptionError("样本文件格式错误 样本列表文件依次为sample_id\tleft.fq\tright.fq", code = "22000602")
        elif self.option('fq_type') == "SE":
            with open(self.option('sample_fq_list'), 'r') as f:
                if len(f.readline().strip().split("\t")) < 2:
                    raise OptionError("样本文件格式错误 样本列表文件依次为sample_id\tleft.fq", code = "22000603")
        return True
        '''


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
            self.logger.info("{} {}".format(self.option('left').prop['path'], self.option('right').prop['path']))
            options.update({
                "fq_l": self.option('left'),
                "fq_r": self.option('right'),
            })
        elif self.option("fq_type") == "SE":
            options.update({
                "fq_s": self.option('single'),
            })
        else:
            pass

        self.trinity.set_options(options)
        self.trinity.run()

    def run_trinity_dis(self):
        '''
        分布投递Trinity命令
        '''
        for f in os.listdir(self.trinity.output_dir):
            if f.endswith("completed"):
                continue
            options={
                "distribute_cmd": self.trinity.output_dir + "/" + f,
                "cpu": 2,
            }
            trinity_dis_tool = self.add_tool('denovo_rna_v2.trinity2_distribute')
            trinity_dis_tool.set_options(options)
            self.trinity_dis.append(trinity_dis_tool)
        self.on_rely(self.trinity_dis, self.run_trinity_merge)
        for tool in self.trinity_dis:
            tool.run()

    def run_trinity_merge(self):
        '''
        分布投递结果合并
        '''
        options = {
            "partion_dir": self.trinity.option('partion_path')
        }
        self.trinity_merge.set_options(options)
        self.trinity_merge.on('end', self.set_output, 'trinity')
        self.trinity_merge.run()

    def run(self):
        '''
        设置运行顺序
        '''
        super(DenovoTrinityModule, self).run()
        self.trinity.on('end', self.run_trinity_dis)
        self.run_trinity()

    def rename(self, trinity_fa, gene2trans):
        ass_fa = trinity_fa
        with open(self.work_dir + '/Trinity.fasta', 'w') as fa_fo:
            for seq in SeqIO.parse(ass_fa, "fasta"):
                new_name = self.option("name") + seq.name.split("TRINITY")[-1]
                fa_fo.write(">{}\n{}\n".format(new_name, seq.seq))
        with open(gene2trans, 'r') as g2t_fi, open(self.work_dir + '/Trinity.fasta.gene_trans_map', 'w') as g2t_fo:
            for line in g2t_fi:
                g2t_fo.write(line.replace("TRINITY", self.option("name")))
        return self.work_dir + '/Trinity.fasta', self.work_dir + '/Trinity.fasta.gene_trans_map'


    def set_output(self, event):
        '''
        设置结果目录
        '''
        obj = event['bind_object']
        self.logger.info('设置目录 {}'.format(event['data']))

        # self.linkdir(obj.output_dir, event['data'])

        if event['data'] == 'trinity':

            trinity_fa = os.path.join(obj.output_dir, 'Trinity.fasta')
            out = os.path.join(self.output_dir, 'Trinity.fasta')
            gene2trans = os.path.join(obj.output_dir, 'Trinity.fasta.gene_trans_map')
            out1 = os.path.join(self.output_dir, 'Trinity.fasta.gene_trans_map')
            trinity_fa, gene2trans = self.rename(trinity_fa, gene2trans)
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
            self.set_error('需要移动到output目录的文件夹不存在。', code = "22001601")
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

        super(DenovoTrinityModule, self).end()

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
            "id": "denovo_trinity" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "denovo_rna_v2.denovo_trinity",
            "instant": False,
            "options": dict(
                left="/mnt/ilustre/users/sanger-dev/workspace/20190617/Single_FastqMerge1089/FastqMerge/output/CL.merge.left.fq",
                right="/mnt/ilustre/users/sanger-dev/workspace/20190617/Single_FastqMerge1089/FastqMerge/output/CL.merge.right.fq",
                lines="20"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
