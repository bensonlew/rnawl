#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,time
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module



class GcContentModule(Module):
    """
    次级代谢产物的预测
    last_modify: 2018.04.14
    """

    def __init__(self, work_id):
        super(GcContentModule, self).__init__(work_id)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"}
        ]
        self.step.add_steps('gc_content')
        self.add_option(options)
        self.modules = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option('fasta_dir').is_set:
            raise OptionError("请传入序列文件夹！", code="")

    def get_list(self):
        if os.path.exists(self.option("fasta_dir").prop['path'] + '/' + "list.txt"):
            sample_dict = {}
            with open(self.option("fasta_dir").prop['path'] + '/' + "list.txt") as t:
                for i in t.readlines()[1:]:
                    if len(i.split("\t")) != 2:
                        raise OptionError("list文件应为两列数据！", code="")
                    if i.split("\t") not in sample_dict.keys():
                        sample_dict[i.split("\t")[0]] = [i.split("\t")[1].strip()]
                    else:
                        sample_dict[i.split("\t")[0]].append(i.split("\t")[1].strip())
        else:
            raise OptionError("list文件不存在！", code="")
        return sample_dict

    def run_gc_content(self):
        """
        fasta文件夹GC统计
        """
        if self.option('fasta_dir').is_set:
            self.sample_dict = self.get_list()
            n = 0
            for sample in self.sample_dict.keys():
                self.gc_content = self.add_tool('tool_lab.gc_content')
                self.step.add_steps('gc_content{}'.format(n))
                if len(self.sample_dict[sample]) >1:
                    cat_file = self.work_dir + '/' + sample + '_cat.fasta'
                    os.system("cat {} >> {}".format(" ".join(self.sample_dict[sample]),cat_file))
                    time.sleep(10)
                    opts = {
                        "fasta": cat_file,
                        "sample_name": sample
                    }
                    self.gc_content.set_options(opts)
                    step = getattr(self.step, 'gc_content{}'.format(n))
                    step.start()
                    self.step.update()
                    self.gc_content.on('end',  'gc_content{}'.format(n))
                    self.modules.append(self.gc_content)
                    n += 1
                else:
                    opts = {
                        "fasta": self.option("fasta_dir").prop['path'] + '/' + self.sample_dict[sample][0],
                        "sample_name": sample
                    }
                    self.gc_content.set_options(opts)
                    step = getattr(self.step, 'gc_content{}'.format(n))
                    step.start()
                    self.step.update()
                    self.gc_content.on('end', 'gc_content{}'.format(n))
                    self.modules.append(self.gc_content)
                    n += 1
            self.logger.info(self.modules)
            if len(self.modules) >1:
               self.on_rely(self.modules, self.set_output)
            elif len(self.modules) == 1:
                self.modules[0].on('end',self.set_output)
            self.step.update()
            for module in self.modules:
                module.run()

    def run(self):
        self.run_gc_content()
        super(GcContentModule, self).run()

    def set_output(self):
        self.logger.info("设置结果目录")
        all_result_file = self.output_dir + '/all.result.xls'
        if os.path.exists(self.work_dir + "/tmp"):
            shutil.rmtree(self.work_dir + "/tmp")
        os.mkdir(self.work_dir + "/tmp")
        for module in self.modules:
            files =os.listdir(module.output_dir)
            for file in files:
                print file
                print module.output_dir + '/' + file
                if os.path.exists(self.work_dir + '/tmp/' + file):
                    os.remove(self.work_dir + '/tmp/' + file)
                os.link(module.output_dir + '/' + file,self.work_dir + '/tmp/' + file)
        with open(all_result_file,"w") as t:
            t.write("sample name" + "\t" + "seq count" + "\t" + "GC content (%)" + "\t" + "N ratio(%)" + "\t" + "AT/GC ratio" + "\n")
            for sample in self.sample_dict.keys():
                if os.path.exists(self.work_dir + '/tmp/' + sample + '.all.result.xls'):
                    with open(self.work_dir + '/tmp/' + sample + '.all.result.xls') as g:
                            t.write(g.readlines()[1])
                    if os.path.exists(self.output_dir + '/' + sample + '.detail.result.xls'):
                        os.remove(self.output_dir + '/' + sample + '.detail.result.xls')
                    os.link(self.work_dir + '/tmp/' + sample + '.detail.result.xls', self.output_dir + '/' + sample + '.detail.result.xls')
        self.end()

    def end(self):
        super(GcContentModule, self).end()