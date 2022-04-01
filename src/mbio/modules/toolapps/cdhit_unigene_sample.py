# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:2017.08.22

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class CdhitUnigeneSampleModule(Module):
    def __init__(self, work_id):
        super(CdhitUnigeneSampleModule, self).__init__(work_id)
        options = [
            {"name": "gene_tmp_fa", "type": "infile", "format": "sequence.fasta"},  # 输出改名并合并的序列
            {"name": "number", "type": "int", "default": 0},  # 切分为几份，默认0表示按文件大小自动计算，指定某个整数时则按指定数量分割
            {"name": "size", "type": "int", "default": 500000000},  # 切分文件大小，当不选择number时，根据size进行切分，默认300M一份  add by zouxuan 2017.11.28
            {"name": "out_dir", "type": "string"},  # 输出路径
            # {"name": "uni_fasta", "type": "outfile", "format": "sequence.fasta"},  # 非冗余基因集核酸序列
            # {"name": "uni_fastaa", "type": "outfile", "format": "sequence.fasta"},  # 非冗余基因集蛋白序列
            {"name": "identity", "type": "float", "default": 0.95},  # 给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
            {"name": "ana_type", "type": "string", "default": "nucl"},  # 输入分析类型，是对核酸聚类还是随蛋白聚类
        ]
        self.add_option(options)
        self.split = self.add_tool("sequence.cdhit_split_fasta")
        self.single = self.add_tool("cluster.cdhit_compare_single")
        # self.merge = self.add_tool("cluster.cdhit_merge")
        self.para = []
        # self.length_tool = self.add_tool("sequence.length_distribute")
        self.step.add_steps('split', 'single')

    def check_options(self):
        if not 0.75 <= self.option("identity") <= 1:
            raise OptionError("identity必须在0.75，1之间", code="21600401")
        if not 0 <= self.option("coverage") <= 1:
            raise OptionError("coverage必须在0,1之间", code="21600402")
        if self.option("number") < 0:
            raise OptionError("number必须大于等于0", code="21600403")
        if self.option("number") == 0:
            self.number = os.path.getsize(self.option("gene_tmp_fa").prop['path']) / self.option("size") + 1
        else:
            self.number = self.option("number")

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def div_fasta(self):
        self.split.set_options({"gene_tmp_fa": self.option('gene_tmp_fa'),
                                "ou_dir": self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp',
                                "number": self.number
                                })
        self.logger.info(self.split)
        self.split.on("start", self.set_step, {'start': self.step.split})
        self.split.on("end", self.set_step, {'end': self.step.split})
        self.split.on("end", self.single_compare)
        self.split.run()

    def single_compare(self):
        self.single.set_options(
            {"query": self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp' + "/gene.geneset.tmp.fa.div-0",
             "qunum": 0,
             "compare": self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp',
             "identity": self.option("identity"),
             "coverage": self.option("coverage"),
             "ana_type": self.option("ana_type"),
             "memory_limit": 40000,
             "num_threads": 4
             })
        self.single.on("start", self.set_step, {'start': self.step.single})
        self.single.on("end", self.set_step, {'end': self.step.single})
        if self.number > 1:
            self.single.on("end", self.add_para)
        else:
            self.single.on("end", self.end)
        self.single.run()

    #    def run_para(self,event):
    #        self.para[event["data"]].run

    def add_para(self):
        n = 1
        for i in range(1, self.number):
            para = self.add_module("toolapps.cdhit_para")
            self.step.add_steps('para_{}'.format(n))
            opts = {
                "first": i,
                "last": self.number,
                "in_dir": self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp',
                "coverage": self.option("coverage"),
                "identity": self.option("identity"),
                "ana_type": self.option("ana_type")
            }
            para.set_options(opts)
            step = getattr(self.step, 'para_{}'.format(n))
            step.start()
            para.on("end", self.finish_update, "para_{}".format(n))
            n += 1
            self.para.append(para)
            if i >= 2:
                self.para[i - 2].on("end", self.para[i - 1].run)
        if len(self.para) == 1:
            self.para[0].on('end', self.end)
        else:
            self.para[self.number - 2].on("end", self.end)
        self.para[0].run()

    # def merge_run(self):
    #     opts = {
    #         "compare_dir": self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp',
    #     }
    #     self.merge.set_options(opts)
    #     self.merge.on("start", self.set_step, {'start': self.step.merge})
    #     self.merge.on("end", self.set_step, {'end': self.step.merge})
    #     self.merge.on('end', self.length_state)
    #     self.merge.run()
    #
    # def length_state(self):
    #     opts = {
    #         "fasta_dir": os.path.split(self.merge.option("fa").prop['path'])[0],
    #         "len_range": "200,300,400,500,600,800",
    #     }
    #     self.length_tool.set_options(opts)
    #     self.length_tool.on("start", self.set_step, {'start': self.step.length})
    #     self.length_tool.on("end", self.set_step, {'end': self.step.length})
    #     self.length_tool.on('end', self.set_output)
    #     self.length_tool.run()
    #
    # def set_output(self):
    #     self.linkdir(self.merge.output_dir, os.path.join(self.output_dir, 'uniGeneset'))
    #     self.linkdir(self.length_tool.output_dir, os.path.join(self.output_dir, 'length_distribute'))
    #     self.option('uni_fasta', self.merge.option("fa"))
    #     self.option('uni_fastaa', self.merge.option("faa"))
    #     self.end()
    #
    # def linkdir(self, dirpath, dirname):
    #     """
    #     link一个文件夹下的所有文件到本module的output目录
    #     """
    #     allfiles = os.listdir(dirpath)
    #     newdir = os.path.join(self.output_dir, dirname)
    #     if not os.path.exists(newdir):
    #         os.mkdir(newdir)
    #     oldfiles = [os.path.join(dirpath, i) for i in allfiles]
    #     newfiles = [os.path.join(newdir, i) for i in allfiles]
    #     for newfile in newfiles:
    #         if os.path.exists(newfile):
    #             os.remove(newfile)
    #     for i in range(len(allfiles)):
    #         os.link(oldfiles[i], newfiles[i])

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"],
        #     ["uniGeneset", "", "非冗余基因集输出目录"],
        #     ["uniGeneset/geneCatalog_stat.xls", "xls", "非冗余基因集统计结果"],
        #     ["uniGeneset/gene.uniGeneset.fa", "fa", "非冗余基因集核酸序列"],
        #     ["uniGeneset/gene.uniGeneset.faa", "faa", "非冗余基因集蛋白序列"],
        #     ["length_distribute", "", "非冗余基因集长度分布统计目录"],
        # ])
        # result_dir.add_regexp_rules([
        #     [r'length_distribute/gene_step_.*\.txt$', 'txt', '长度分布统计结果']
        # ])
        super(CdhitUnigeneSampleModule, self).end()

    def run(self):
        super(CdhitUnigeneSampleModule, self).run()
        # self.on_rely(self.para,self.merge_run)
        if self.number > 1:
            self.div_fasta()
        else:
            if os.path.exists(self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp'):
                pass
            else:
                os.mkdir(self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp')
            if os.path.exists(self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp/gene.geneset.tmp.fa.div-0'):
                os.remove(self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp/gene.geneset.tmp.fa.div-0')
            shutil.copyfile(self.option("gene_tmp_fa").prop['path'],
                            self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp/gene.geneset.tmp.fa.div-0')
            self.single_compare()
