#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import unittest
import shutil
from collections import OrderedDict


class SrnaModule(Module):
    """
    对所有样本进行定量
    """
    def __init__(self, work_id):
        super(SrnaModule, self).__init__(work_id)
        options = [
            dict(name="fna", type="string"),
            dict(name="input_file", type="string"),
            dict(name="type", type="string", default="feature"),
            dict(name="group_list", type="string"),
            dict(name="trimPairFq", type="string"),
            dict(name="evalue", type="string", default="0.00001"),
            dict(name="special", type="string", default="true"),
            dict(name="rfam", type="string", default="infernal"),
        ]
        self.add_option(options)
        self.rockhopper = self.add_tool('prok_rna.rockhopper')
        self.annot = self.add_tool('prok_rna.sRNA_annotation')
        self.fold = self.add_tool('prok_rna.sRNA_fold')
        self.target = self.add_tool('prok_rna.sRNA_target')
        self.out_time = 0

    def check_options(self):
        if self.option("type").lower() not in ["feature", "gff", "gtf"]:
            raise OptionError("the type of infile is incorrect", code = "25001201")
        if not os.path.exists(self.option('fna')):
            raise OptionError("基因组文件不存在", code = "25001202")
        if not os.path.exists(self.option('input_file')):
            raise OptionError("feature文件不存在", code = "25001203")
        if not os.path.exists(self.option('group_list')):
            raise OptionError("group_list文件不存在", code = "25001204")
        if not os.path.exists(self.option('trimPairFq')):
            raise OptionError("trimPairFq_list文件不存在", code = "25001205")
        if not 0 <= float(self.option('evalue')) <= 1:
            raise OptionError("evalue的值应该在0至1之间", code = "25001206")
        if self.option("special").lower() not in ["true", "false"]:
            raise OptionError("只有特异性和非特异性两种选项", code = "25001207")

    def set_step(self, event):
        if 'start' in event['data']:
            event['data']['start'].start()
        if 'end' in event['data']:
            event['data']['end'].finish()
        self.step.update()

    def run_rockhopper(self):
        options = {
            'fna': self.option('fna'),
            'input_file': self.option('input_file'),
            'type': self.option('type'),
            'group_list': self.option('group_list'),
            'trimPairFq': self.option('trimPairFq'),
            'special': self.option('special'),
        }
        self.rockhopper.set_options(options)
        # self.rockhopper.on('start', self.set_step, {'start': self.step.rockhopper})
        # self.rockhopper.on('end', self.set_step, {'end': self.step.rockhopper})
        self.rockhopper.on('end', self.set_output, 'rockhopper')
        if self.option('rfam') == 'infernal':
            self.rockhopper.on('end', self.run_rfam)
        else:
            self.rockhopper.on('end', self.run_annot)
        self.rockhopper.on('end', self.run_fold)
        self.rockhopper.on('end', self.run_target)
        self.rockhopper.run()

    def run_rfam(self):
        self.rfam = self.add_tool('prok_rna.infernal_rfam')
        options = {
            'query': self.rockhopper.option('predict_fa').prop['path'],
            'evalue': self.option('evalue'),
        }
        self.rfam.set_options(options)
        self.rfam.on('end', self.run_annot)
        self.rfam.run()

    def run_annot(self):
        options = {
            'predict_fa': self.rockhopper.option('predict_fa').prop['path'],
            'evalue': self.option('evalue'),
        }
        if self.option('rfam') == 'infernal':
            options.update({'rfam': self.rfam.option('rfam_result').prop['path']})
        self.annot.set_options(options)
        # self.annot.on('start', self.set_step, {'start': self.step.annot})
        # self.annot.on('end', self.set_step, {'end': self.step.annot})
        self.annot.on('end', self.set_output, 'srna_annot')
        self.annot.run()

    def run_fold(self):
        options = {
            'predict_fa': self.rockhopper.option('predict_fa').prop['path'],
        }
        self.fold.set_options(options)
        # self.fold.on('start', self.set_step, {'start': self.step.fold})
        # self.fold.on('end', self.set_step, {'end': self.step.fold})
        self.fold.on('end', self.set_output, 'srna_fold')
        self.fold.run()

    def run_target(self):
        options = {
            'predict_fa': self.rockhopper.option('predict_fa').prop['path'],
            'genome_bed': self.rockhopper.option('genome_bed').prop['path'],
            'genome_fa': self.option('fna'),
        }
        self.target.set_options(options)
        # self.target.on('start', self.set_step, {'start': self.step.target})
        # self.target.on('end', self.set_step, {'end': self.step.target})
        self.target.on('end', self.set_output, 'srna_target')
        self.target.run()

    def run(self):
        # tools = [self.rockhopper, self.annot, self.fold, self.target]
        # self.on_rely(tools, self.end())
        super(SrnaModule, self).run()
        self.run_rockhopper()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'rockhopper':
            self.linkdir(obj.output_dir, 'rockhopper')
            self.out_time += 1
        elif event['data'] == 'srna_annot':
            self.linkdir(obj.output_dir, 'srna_annot')
            self.out_time += 1
        elif event['data'] == 'srna_fold':
            self.linkdir(obj.output_dir, 'srna_fold')
            self.out_time += 1
        elif event['data'] == 'srna_target':
            self.linkdir(obj.output_dir, 'srna_target')
            self.out_time += 1
        if self.out_time == 4:
            self.end()

    def linkdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        if not os.path.isdir(olddir):
            self.set_error("需要移动到output目录的文件夹不存在", code = "25001208")
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
                os.system("cp -r {} {}".format(oldfiles[i], new1))

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "sRNA分析结果目录"],
        ])
        super(SrnaModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Srna_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "prok_rna.srna",
            "instant": False,
            "options": dict(
                fna="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/pipline/ref/GCF_000009345.1_ASM934v1_genomic.fna",
                input_file="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/pipline/ref/GCF_000009345.1_ASM934v1_feature_table.txt",
                type="feature",
                group_list='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/pipline/ref/index_test/group_list',
                trimPairFq='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/pipline/trimPairFq.list',
                evalue='0.001',
                special = 'true'
            )
        }
        # data['options']['method'] = 'rsem'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        # #
        # data['id'] += '1'
        # data['options']['method'] = 'salmon'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        #
        data['id'] += '_fyt'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
