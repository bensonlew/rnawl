# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# @20200226

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.toolapps.common import link_dir


class CdhitUnigeneModule(Module):
    """
    宏转录组小工具：cdhit去冗余
    数据可以是单个核酸文件或者单个蛋白文件或者是所有序列的文件夹，需要指定是对蛋白还是核酸去冗余
    上传文件可以压缩也可以不压缩
    """
    def __init__(self, work_id):
        super(CdhitUnigeneModule, self).__init__(work_id)
        options = [
            {"name": "gene_fa_dir", "type": "infile", "format": "toolapps.fasta_dir"},  # 基因预测核酸文件
            {"name": "gene_faa_dir", "type": "infile", "format": "toolapps.fasta_dir"},  # 基因预测蛋白文件
            {"name": "cdhit_identity", "type": "float", "default": 0.9},  # 给出cdhit的参数identity
            {"name": "cdhit_coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
            {"name": "ana_type", "type": "string", "default": "nucl"},  # 输入分析类型，对蛋白还是对核酸聚类
            {"name": "clstr", "type": "bool", "default": True},  # 是否生成cluster文件，False不生成，True生成
            # {'name': 'is_unzip', 'type': 'string', 'default': 'false'},  # 是否是压缩格式， 默认为非压缩格式
        ]
        self.add_option(options)
        self.cdhit = self.add_module("toolapps.cdhit_unigene_sample")
        self.merge = self.add_tool("toolapps.cdhit_merge")
        self.cat_reads = self.add_tool("sequence.cat_table")
        self.translate = self.add_tool("toolapps.translate")
        self.step.add_steps("merge", "cdhit", "cdhit_merge")
        self.sequence_list = []

    def check_options(self):
        """
        进行参数二次检查
        """
        if not 0.75 <= self.option("cdhit_identity") <= 1:
            raise OptionError("cdhit identity必须在0.75，1之间", code="24400601")
        if not 0 <= self.option("cdhit_coverage") <= 1:
            raise OptionError("cdhit coverage必须在0,1之间", code="24400602")
        if self.option("ana_type") in ['nucl']:
            if (not self.option("gene_fa_dir").is_set):
                raise OptionError("必须提供分析所需的核酸文件夹", code="24400603")
        elif self.option("ana_type") in ['prot']:
            if (not self.option("gene_fa_dir").is_set) and (not self.option("gene_faa_dir").is_set):
                raise OptionError("必须提供分析所需的蛋白序列文件夹或者核酸文件夹", code="24400604")
        else:
            raise OptionError("不存在的分析类型", code="24400605")


    def merge_sequence(self):
        """
        合并sequence_dir文件夹中的所有序列，序列格式为fasta
        """
        self.logger.info("start cat fasta")
        opts = {
            "prefix": "all.gene"
            }
        if self.option("gene_fa_dir").is_set:
            list_dirs = os.listdir(self.option("gene_fa_dir").prop['fasta_dir'])
            if "list.txt" in list_dirs:
                list_path = os.path.join(self.option("gene_fa_dir").prop['fasta_dir'], 'list.txt')
                os.remove(list_path)
                opts["fa_dir"] = self.option("gene_fa_dir").prop['fasta_dir']
        elif self.option("gene_faa_dir").is_set:
            list_dirs = os.listdir(self.option("gene_faa_dir").prop['fasta_dir'])
            if "list.txt" in list_dirs:
                list_path = os.path.join(self.option("gene_faa_dir").prop['fasta_dir'], 'list.txt')
                os.remove(list_path)
                opts["fa_dir"] = self.option("gene_faa_dir").prop['fasta_dir']
        self.cat_reads.set_options(opts)
        self.cat_reads.run()

    def run_translate(self):
        """
        对上传的序列进行翻译
        :return:
        """
        if self.option("gene_fa_dir").is_set:
            fa_path = os.path.join(self.cat_reads.output_dir, "all.gene.xls")
        self.logger.info("start translate!")
        if os.path.exists(os.path.join(self.work_dir, 'all.gene.fa')):
            os.remove(os.path.join(self.work_dir, 'all.gene.fa'))
            os.link(fa_path, os.path.join(self.work_dir, 'all.gene.fa'))
        else:
            os.link(fa_path, os.path.join(self.work_dir, 'all.gene.fa'))
        opts = {
            "fasta" : os.path.join(self.work_dir, 'all.gene.fa'),
            'is_dir': 'false'
            }
        self.translate.set_options(opts)
        self.translate.run()

    def run_cd_hit(self):
        """
        cdhit去冗余，可以用两种不同的方法进行分析（nucl和prot）
        """
        self.logger.info("start cdhit")
        if self.option("ana_type") in ["nucl"]: ##设置切分文件大小为500M
            if self.option("gene_fa_dir").is_set:
                self.number1 = os.path.getsize(os.path.join(self.cat_reads.output_dir, 'all.gene.xls')) / 300000000 + 1
        else:
            if self.option("gene_faa_dir").is_set:
                self.number1 = os.path.getsize(os.path.join(self.cat_reads.output_dir, 'all.gene.xls')) / 300000000 + 1
            else:
                self.number1 = os.path.getsize(os.path.join(self.translate.output_dir, 'gene.uniGeneset.faa')) / 300000000 + 1
        opts = {
            "number": int(self.number1),
            "out_dir": self.work_dir,
            "identity": self.option("cdhit_identity"),
            "coverage": self.option("cdhit_coverage"),
            "ana_type": self.option("ana_type"),
            }
        if self.option("ana_type") in ["nucl"]:
            if self.option("gene_fa_dir").is_set:
                if os.path.exists(os.path.join(self.work_dir, 'all.gene.fa')):
                    os.remove(os.path.join(self.work_dir, 'all.gene.fa'))
                    os.link(os.path.join(self.cat_reads.output_dir, 'all.gene.xls'), os.path.join(self.work_dir, 'all.gene.fa'))
                else:
                    os.link(os.path.join(self.cat_reads.output_dir, 'all.gene.xls'), os.path.join(self.work_dir, 'all.gene.fa'))
                opts["gene_tmp_fa"] = os.path.join(self.work_dir, 'all.gene.fa')
        else:
            if self.option("gene_faa_dir").is_set:
                if os.path.exists(os.path.join(self.work_dir, 'all.gene.faa')):
                    os.remove(os.path.join(self.work_dir, 'all.gene.faa'))
                    os.link(os.path.join(self.cat_reads.output_dir, 'all.gene.xls'), os.path.join(self.work_dir, 'all.gene.faa'))
                else:
                    os.link(os.path.join(self.cat_reads.output_dir, 'all.gene.xls'), os.path.join(self.work_dir, 'all.gene.faa'))
                opts["gene_tmp_fa"] = os.path.join(self.work_dir, 'all.gene.faa')
            else:
                input_fa = os.path.join(self.translate.output_dir, 'gene.uniGeneset.faa')
                new_input_faa = os.path.join(self.work_dir, 'all.gene.faa')
                if os.path.exists(new_input_faa):
                    os.remove(new_input_faa)
                    os.link(input_fa, new_input_faa)
                else:
                    os.link(input_fa, new_input_faa)
                opts["gene_tmp_fa"] = new_input_faa
        self.cdhit.set_options(opts)
        self.cdhit.run()

    def run_merge(self):
        """
        合并去冗余后的序列文件
        :return:
        """
        self.logger.info("star merge fa and translate")
        opts = {
            "compare_dir": self.work_dir + '/gene.uniGeneset.fa.cd-hit-para-tmp',
            "num1":self.number1,
            'ana_type': self.option("ana_type"),
            'table': 11
            }
        if self.option("clstr"):
            opts["clstr"] =  1
        else:
            opts["clstr"] =  0
        if self.option("ana_type") in ['prot']:
            if self.option("gene_faa_dir").is_set:
                opts["is_trans"] = "false"
            else:
                if self.option("gene_fa_dir").is_set:
                    opts["is_trans"] = "true"
                    opts["gene_fa"] = os.path.join(self.work_dir, 'all.gene.fa')
        else:
            if self.option("gene_faa_dir").is_set:
                opts["is_trans"] = "false"
                opts['gene_faa'] = os.path.join(self.work_dir, 'all.gene.faa')
            else:
                if self.option("gene_fa_dir").is_set:
                    opts["gene_fa"] = os.path.join(self.work_dir, 'all.gene.fa')
                    opts["is_trans"] = "true"
        self.merge.set_options(opts)
        self.merge.run()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        output_path = os.path.join(self.output_dir, "cdhit")
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        else:
            os.mkdir(output_path)
        link_dir(self.merge.output_dir, output_path)
        self.logger.info("设置结果文件目录完成")
        self.end()

    def run(self):
        """
        运行开始
        :return:
        """
        super(CdhitUnigeneModule, self).run()
        merge_file = os.path.join(self.work_dir, 'gene.uniGeneset.fa.cd-hit-para-tmp')
        if os.path.exists(merge_file):
            pass
        else:
            os.mkdir(merge_file)
        if self.option("ana_type") in ['nucl']: ##用nucl核酸进行构建非冗余基因集
            if self.option("gene_fa_dir").is_set:##已经解压过了
                self.cat_reads.on("end", self.run_cd_hit)
                self.cdhit.on("end", self.run_merge)
                self.merge.on("end", self.set_output)
                self.merge_sequence()
        else: ##用prot蛋白进行构建非冗余基因集
            if self.option("gene_faa_dir").is_set:
                self.cat_reads.on("end", self.run_cd_hit)
                self.cdhit.on("end", self.run_merge)
            else:
                if self.option("gene_fa_dir").is_set:
                    self.cat_reads.on("end", self.run_translate)
                    self.translate.on("end", self.run_cd_hit)
                    self.cdhit.on("end", self.run_merge)
                else:
                    raise OptionError("没有提供正确的文件夹")
            self.merge.on("end", self.set_output)
            self.merge_sequence()

    def end(self):
        """
        运行结束
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        self.logger.info("Now start upload!")
        result_dir.add_relpath_rules([
            [".", "", "非冗余基因集输出目录"],
            ["./geneCatalog_stat.xls", "xls", "非冗余基因集统计结果"],
            ["./gene.uniGeneset.fa", "fa", "非冗余基因集核酸序列"],
            ["./gene.uniGeneset.faa", "faa", "非冗余基因集蛋白序列"],
        ])
        super(CdhitUnigeneModule, self).end()

