# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:2018.01.09

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import unittest


class UniGeneModule(Module):
    def __init__(self, work_id):
        super(UniGeneModule, self).__init__(work_id)
        options = [
            {"name": "gene_tmp_fa", "type": "infile", "format": "sequence.fasta"},  # 样品单拼序列
            {"name": "gene_tmp_fa_mix", "type": "infile", "format": "sequence.fasta"},  # 混拼序列
            {"name": "number1", "type": "int", "default": 0},  # 样品单拼序列切分为几份，默认0表示按文件大小自动计算，指定某个整数时则按指定数量分割
            {"name": "number2", "type": "int", "default": 0},  # 混拼序列切分为几份，默认0表示按文件大小自动计算，指定某个整数时则按指定数量分割
            {"name": "uni_fasta", "type": "outfile", "format": "sequence.fasta"},  # 非冗余基因集核酸序列
            {"name": "uni_fastaa", "type": "outfile", "format": "sequence.fasta"},  # 非冗余基因集蛋白序列
            {"name": "cdhit_identity", "type": "float", "default": 0.9},  # 给出cdhit的参数identity
            {"name": "cdhit_coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
            {"name": "insertsize", "type": "infile", "format": "sample.insertsize_table"},  # 插入片段文件
            {"name": "QC_dir", "type": "infile", "format": "sequence.fastq_dir"},  # qc后reads文件夹
            {"name": "reads_abundance", "type": "outfile", "format": "sequence.profile_table"},  # reads_abundance
            {"name": "rpkm_abundance", "type": "outfile", "format": "sequence.profile_table"},  # rpkm_abundance
            {"name": "gene_length_list", "type": "outfile", "format": "sequence.profile_table"},
            # gene长度统计表  add by GHD @ 20180109
            {"name": "seed", "type": "int", "default": 35},
            # align the initial n bps as a seed means whole lengths of read
            {"name": "mode", "type": "int", "default": 4},
            # match mode for each read or the seed part of read, which shouldn't contain more than 2 mismatches: 0 for exact mathc only; 1 for 1 mismatch; 2 for 2 mismatch; 4 for find the best hits
            {"name": "processors", "type": "int", "default": 6},
            {"name": "mismatch", "type": "int", "default": 20},  # maximum number of mismatches allowed on a read
            {"name": "repeat", "type": "int", "default": 1},  # how to report repeat hits, 0=none, 1=random one, 2=all
            {"name": "soap_identity", "type": "float", "default": 0.95},  # soap aligner identity
            {"name": "map_type", "type": "int", "default": 1},  # 默认1从拼接开始，type2从非冗余基因开始
            {"name": "ana_type", "type": "string", "default": "nucl"},  # 输入分析类型，是对核酸聚类还是随蛋白聚类add by qingchen.zhang@20190423
            {"name": "gene_tmp_faa", "type": "infile", "format": "sequence.fasta"},  # 基因预测蛋白序列add by qingchen.zhang@20190428
            {"name": "gene_tmp_faa_mix", "type": "infile", "format": "sequence.fasta"},  # 基因预测的蛋白序列add by qingchen.zhang@20190428
        ]
        self.add_option(options)
        self.cdhit = self.add_module("cluster.cdhit_unigene_sample")
        self.cdhit2 = self.add_module("cluster.cdhit_unigene_mix")
        self.soap_aligner = self.add_module("align.split_map")
        self.merge = self.add_tool("cluster.cdhit_merge")
        self.length_tool = self.add_tool("sequence.length_distribute")
        self.stat = self.add_tool("statistical.unigene_stat")
        self.step.add_steps("cdhit", "cdhit2", "soap", "merge", "length", "stat")

    def check_options(self):
        if self.option("map_type") == 1:
            if not self.option("gene_tmp_fa").is_set:
                raise OptionError("必须提供样品单拼序列", code="21600501")
            if not 0.75 <= self.option("cdhit_identity") <= 1:
                raise OptionError("cdhit identity必须在0.75，1之间", code="21600502")
            if not 0 <= self.option("cdhit_coverage") <= 1:
                raise OptionError("cdhit coverage必须在0,1之间", code="21600503")
            if self.option("number1") < 0:
                raise OptionError("number1必须大于等于0", code="21600504")
            if self.option("number2") < 0:
                raise OptionError("number2必须大于等于0", code="21600505")
            if not self.option("repeat") in [0, 1, 2]:
                raise OptionError("repeat必须为0,1,或2", code="21600506")
            if not self.option("mode") in [0, 1, 2, 4]:
                raise OptionError("repeat必须为0,1,2,或4", code="21600507")
            if not 0 < self.option("seed") <= 256:
                raise OptionError('seed参数必须设置在1-256之间:%s', variables=(self.option('seed')), code="21600508")
            if not 0 < self.option("soap_identity") < 1:
                raise OptionError("soap identity必须在0，1之间", code="21600509")
            if not self.option("QC_dir").is_set:
                raise OptionError("必须提供质控后的fq文件夹", code="21600510")
            if not self.option("insertsize").is_set:
                raise OptionError("必须提供插入片段文件", code="21600511")


    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def cd_hit(self):
        self.cdhit.set_options({
            "gene_tmp_fa": self.option("gene_tmp_faa"),
            "number": self.number1,
            "out_dir": self.work_dir,
            "identity": self.option("cdhit_identity"),
            "coverage": self.option("cdhit_coverage"),
            "ana_type": "prot",
        })
        self.cdhit.on("start", self.set_step, {"start": self.step.cdhit})
        self.cdhit.on("end", self.set_step, {"end": self.step.cdhit})
        self.cdhit.on("end", self.set_step, {"start": self.step.soap})
        self.cdhit.run()

    def cd_hit2(self):
        self.cdhit2.set_options({
            "gene_tmp_fa": self.option("gene_tmp_faa_mix"),
            "number": self.number2,
            "out_dir": self.work_dir,
            "identity": self.option("cdhit_identity"),
            "coverage": self.option("cdhit_coverage"),
            "ana_type": "prot",
        })
        self.cdhit2.on("start", self.set_step, {"start": self.step.cdhit2})
        self.cdhit2.on("end", self.set_step, {"end": self.step.cdhit2})
        self.cdhit2.on("end", self.set_step, {"start": self.step.soap})
        self.cdhit2.run()

    def unigene_stat(self):
        faa_seq = ""
        if self.option("map_type") == 1:
            fa_seq = self.merge.option("fa")
            faa_seq = self.merge.option("faa")
        elif self.option("map_type") == 2:
            fa_seq = self.option("uni_fasta")
        opts = {
            "gene_lenth": self.soap_aligner.output_dir + "/gene_profile/gene.uniGeneset.fa.length.txt",
            "fafile": fa_seq
        }
        if faa_seq:
            opts["faafile"] = faa_seq
        self.stat.set_options(opts)
        self.stat.on("start", self.set_step, {"start": self.step.stat})
        self.stat.on("end", self.set_step, {"end": self.step.stat})
        self.stat.on("end", self.length_state)
        self.stat.run()

    def soap(self):
        if self.option("map_type") == 1:
            fa_seq = self.merge.option("fa")
        elif self.option("map_type") == 2:
            fa_seq = self.option("uni_fasta")
        self.soap_aligner.set_options({
            "fafile": fa_seq,
            "insertsize": self.option("insertsize"),
            "QC_dir": self.option("QC_dir"),
            "seed": self.option("seed"),
            "mode": self.option("mode"),
            "processors": self.option("processors"),
            "mismatch": self.option("mismatch"),
            "repeat": self.option("repeat"),
            "identity": self.option("soap_identity")
        })
        self.soap_aligner.on("start", self.set_step, {"start": self.step.soap})
        self.soap_aligner.on("end", self.set_step, {"end": self.step.soap})
        self.soap_aligner.on("end", self.unigene_stat)
        self.soap_aligner.run()

    def set_output(self):
        # self.linkdir(self.merge.output_dir, os.path.join(self.output_dir, 'uniGeneset'))
        # self.linkdir(self.length_tool.output_dir, os.path.join(self.output_dir, 'length_distribute'))
        self.option('uni_fasta', self.stat.option("fa"))
        self.option('uni_fastaa', self.stat.option("faa"))
        self.linkdir(self.length_tool.output_dir, "length_distribute")
        self.linkdir(self.stat.output_dir, "uniGeneset")
        self.linkdir(self.soap_aligner.output_dir + "/gene_profile", "gene_profile")
        self.option('reads_abundance', self.soap_aligner.option("reads_abundance"))
        self.option('rpkm_abundance', self.soap_aligner.option("rpkm_abundance"))
        self.option('gene_length_list').set_path(
            self.output_dir + '/gene_profile/gene.uniGeneset.fa.length.txt')  # add option('gene_length_list') by GHD @ 20180109
        self.end()

    def merge_run(self):
        if self.option("gene_tmp_fa_mix").is_set:
            opts = {
                "compare_dir": self.work_dir + '/gene.uniGeneset.fa.cd-hit-para-tmp',
                "clstr": 0,
                "num1":self.number1,
                "num2":self.number2,
                "gene_fa_mix": self.option('gene_tmp_fa_mix'),
                "gene_fa": self.option('gene_tmp_fa'),
                "is_trans": False
            }
        else:
            opts = {
                "compare_dir": self.work_dir + '/gene.uniGeneset.fa.cd-hit-para-tmp',
                "clstr": 0,
                "num1":self.number1,
                "num2":self.number2,
                "gene_fa": self.option('gene_tmp_fa'),
                "is_trans": False
            }
        self.merge.set_options(opts)
        self.merge.on("start", self.set_step, {'start': self.step.merge})
        self.merge.on("end", self.set_step, {'end': self.step.merge})
        self.merge.on('end', self.soap)
        self.merge.run()

    def length_state(self):
        opts = {
            "fasta_dir": os.path.split(self.stat.option("fa").prop['path'])[0],
            "len_range": "200,300,400,500,600,800",
        }
        self.length_tool.set_options(opts)
        self.length_tool.on("start", self.set_step, {'start': self.step.length})
        self.length_tool.on("end", self.set_step, {'end': self.step.length})
        self.length_tool.on('end', self.set_output)
        self.length_tool.run()

    def run(self):
        super(UniGeneModule, self).run()
        if self.option("map_type") == 1:
            if os.path.exists(self.work_dir + '/gene.uniGeneset.fa.cd-hit-para-tmp'):
                pass
            else:
                os.mkdir(self.work_dir + '/gene.uniGeneset.fa.cd-hit-para-tmp') #add by zouxuan 防止之后存在跑mix和sample可能会由于文件夹的问题报错
            self.number1 = self.option("number1")
            if self.option("number1") == 0:
                self.number1 = os.path.getsize(self.option("gene_tmp_faa").prop['path']) / 500000000 + 1
            self.number2 = self.option("number2")
            if self.option("gene_tmp_fa_mix").is_set:
                if self.option("number2") == 0:
                    self.number2 = os.path.getsize(self.option("gene_tmp_faa_mix").prop['path']) / 500000000 + 1
                rely = [self.cdhit, self.cdhit2]
                self.on_rely(rely, self.merge_run)
                self.cd_hit2()
            else:
                self.cdhit.on("end", self.merge_run)
            self.cd_hit()
        else:
            self.soap()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["uniGeneset", "", "非冗余基因集输出目录"],
            ["uniGeneset/geneCatalog_stat.xls", "xls", "非冗余基因集统计结果"],
            ["uniGeneset/gene.uniGeneset.fa", "fa", "非冗余基因集核酸序列"],
            ["uniGeneset/gene.uniGeneset.faa", "faa", "非冗余基因集蛋白序列"],
            ["length_distribute", "", "非冗余基因集长度分布统计目录"],
            ["gene_profile", "", "非冗余基因丰度目录"]
        ])
        result_dir.add_regexp_rules([
            ['length_distribute/gene_step_.*\.txt$', 'txt', '长度分布统计结果'],
            ['gene_profile/*\.xls$', 'xls', '基因丰度表']
        ])
        super(UniGeneModule, self).end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "CreatTable" + str(random.randint(1, 10000)),
            "id": "Uni_gene",
            "type": "module",
            "name": "cluster.uni_gene",
            "instant": True,
            "options": dict(
                uni_fasta="/mnt/ilustre/users/sanger-dev/workspace/20181120/MetaGenomic_metag_20181120/UniGene/CdhitMerge/output/gene.uniGeneset.fa",
                insertsize="/mnt/ilustre/users/sanger-dev/workspace/20181120/MetaGenomic_metag_20181120/UniGene/insertsize",
                QC_dir="/mnt/ilustre/users/sanger-dev/workspace/20181120/MetaGenomic_metag_20181120/QcAndStat/output/after_qc_dir",
                map_type=2
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
