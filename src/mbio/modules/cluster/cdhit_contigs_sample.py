# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:2017.08.22

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class CdhitContigsSampleModule(Module):
    def __init__(self, work_id):
        super(CdhitContigsSampleModule, self).__init__(work_id)
        options = [
            {"name": "gene_tmp_fa", "type": "infile", "format": "sequence.fasta"},  # 输出改名并合并的序列
            {"name": "number", "type": "int", "default": 0},  # 切分为几份，默认0表示按文件大小自动计算，指定某个整数时则按指定数量分割
            {"name": "size", "type": "int", "default": 300000000},  # 切分文件大小，当不选择number时，根据size进行切分，默认300M一份
            {"name": "out_dir", "type": "string"},  # 输出路径
            {"name": "identity", "type": "float", "default": 0.95},  # 给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
        ]
        self.add_option(options)
        self.split = self.add_tool("sequence.cdhit_split_fasta")
        self.single = self.add_tool("cluster.cdhit_compare_single")
        self.para = []
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
                                "number": self.number,
                                "order":0,
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
             })
        self.single.on("start", self.set_step, {'start': self.step.single})
        self.single.on("end", self.set_step, {'end': self.step.single})
        if self.number > 1:
            self.single.on("end", self.add_para)
        else:
            self.single.on("end", self.end)
        self.single.run()

    def add_para(self):
        n = 1
        for i in range(1, self.number):
            para = self.add_module("cluster.cdhit_para")
            self.step.add_steps('para_{}'.format(n))
            opts = {
                "first": i,
                "last": self.number,
                "in_dir": self.option('out_dir') + '/gene.uniGeneset.fa.cd-hit-para-tmp',
                "coverage": self.option("coverage"),
                "identity": self.option("identity"),
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

    def end(self):
        super(CdhitContigsSampleModule, self).end()

    def run(self):
        super(CdhitContigsSampleModule, self).run()
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