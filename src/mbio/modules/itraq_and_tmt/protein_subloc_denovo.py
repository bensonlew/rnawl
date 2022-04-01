# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from __future__ import division
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os
import shutil
import unittest
import sys


class ProteinSublocDenovoModule(Module):
    """
    module for denovorna annotation
    """
    def __init__(self, work_id):
        super(ProteinSublocDenovoModule, self).__init__(work_id)
        options = [
            {"name": "fasta", "type": "string", "default": "/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/subloc/exp.fasta"},
            {"name": "subloc", "type": "bool", "default": True},
            {"name": "taxonomy", "type": "string", "default": 'Animals,Fungi,Bacteria'},   # kegg数据库物种分类, Animals/Plants/Fungi/Protists/Archaea/Bacteria
            {"name": "annot_type", "type": "string", "default": "protein" },
        ]
        self.add_option(options)
        self.diamond_nr = self.add_module("itraq_and_tmt.diamond")
        self.go_annot = self.add_tool('itraq_and_tmt.annotation.go_annotation')
        self.sublocs = list()
        species = self.option('taxonomy')
        self.species = species.split(',')
        for specie in self.species:
            subloc = self.add_tool('itraq_and_tmt.annotation.multiloc')
            self.sublocs.append(subloc)


    def check_options(self):
        pass

    def set_step(self, event):
        if 'start' in event['data']:
            event['data']['start'].start()
        if 'end' in event['data']:
            event['data']['end'].finish()
        self.step.update()

    def run_diamond(self):
        self.logger.info("开始运行diamond比对")
        print(self.option("fasta"))
        self.blast_modules = []
        blast_opts = {
            'query': self.option("fasta"),
            'query_type': 'prot',
            'database': "nr",
            'blast': 'blastp',
            'evalue': None,
            'outfmt': 5,
        }
        self.diamond_nr.set_options(blast_opts)
        print(blast_opts)
        self.diamond_nr.on('end', self.run_go_anno)
        self.diamond_nr.run()


    def run_go_anno(self):
        """
        """
        options = {
            'blastout': self.diamond_nr.option('outxml'),
            'protein_fasta': self.option('fasta')
        }
        self.go_annot.set_options(options)
        self.go_annot.on('end', self.run_subloc)
        self.go_annot.run()

    def run_subloc(self):
        go_path = os.path.join(self.go_annot.output_dir, 'query_gos.list')
        options = {
            'fa': self.option('fasta'),
            'go': go_path,
        }
        for n, subloc in enumerate(self.sublocs):
            options.update(
                {
                    'species': self.species[n]
                }
            )
            subloc.set_options(options)
            subloc.run()
        self.on_rely(self.sublocs, self.set_output)

    def run(self):
        super(ProteinSublocDenovoModule, self).run()
        self.run_diamond()

    def set_output(self, event):
        self.end()

    def linkdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        if not os.path.isdir(olddir):
            raise Exception('需要移动到output目录的文件夹不存在。')
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
        super(ProteinSublocDenovoModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import sys

        # fasta = os.path.abspath(sys.argv[1])
        # taxonomy = sys.argv[2]
        fasta = '/mnt/lustre/users/sanger/sg-users/fengyitong/subloc/chendongli_20190319/exp.fasta'
        taxonomy = 'Animals'
        # sub_dir = '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/subloc/'
        data = {
            "id": "protein_subloc_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "itraq_and_tmt.protein_subloc_denovo",
            "instant": False,
            "options": dict(
                fasta=fasta,
                taxonomy=taxonomy,
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
        data['id'] += 'Bl31LJ__PEAKS_9-3xR'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
