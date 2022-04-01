# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20200914

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import glob
import pandas as pd
import unittest
import gevent.subprocess as subprocess
import json
from collections import OrderedDict
from biocluster.config import Config
import shutil
from mbio.packages.medical_transcriptome.copy_file import CopyFile

class DiffReactomeEnrichModule(Module):
    """
    该Module用于基因融合分析，默认使用方法
    """
    def __init__(self, work_id):
        super(DiffReactomeEnrichModule, self).__init__(work_id)
        options = [
            {"name": "task_id", "type": "string", "default": None},
            {"name": "regulate", "type": "string", "default": None},
            {'name': 'geneset_reactome', 'type': 'string'},
            {'name': 'reactome_annot', 'type': 'string'},
            {'name': 'geneset_list', 'type': 'string'},
            {'name': 'all_list', 'type': 'string'},
            {"name": "level", "type": "string", "default": None},
            {'name': 'reactome_version', 'type': 'string', 'default': None},
        ]
        self.add_option(options)
        self.geneset_prepare_tools=[]
        self.genest_infos = OrderedDict()
        # self.common = self.add_tool("medical_transcriptome.diff_geneset.annot_prepare")
        self.total_dict = OrderedDict()


    def check_options(self):
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_reactome_rich(self):
        self.enrich_tool = self.add_tool('medical_transcriptome.diff_geneset.diff_reactome_enrich')
        self.enrich_tool.set_options({
            'reactome_list': self.option('reactome_annot'),
            'diff_list': self.option('geneset_list'),
        })
        self.enrich_tool.on('end', self.run_reactome_class)
        self.enrich_tool.run()

    def run_reactome_class(self):
        self.reactome_class = self.add_tool('medical_transcriptome.diff_geneset.diff_reactome_class')
        self.reactome_class.set_options({
            'geneset_ids': self.option('geneset_reactome'),
            'reactome_annot': self.option('reactome_annot'),
            "reactome_version": self.option('reactome_version'),
        })
        self.reactome_class.on('end', self.set_output)
        self.reactome_class.run()

    def set_output(self):
        self.replace_reactome_link()
        svg_path = self.reactome_class.output_dir + '/svg'
        svg_gz_path = self.reactome_class.output_dir + '/svg.tar.gz'
        if os.path.exists(os.path.join(self.output_dir, "svg")):
            shutil.rmtree(os.path.join(self.output_dir, "svg"))
        try:
            CopyFile().linkdir(svg_path, os.path.join(self.output_dir, "svg"))
            # os.system('ln -s {} {}'.format(svg_path, os.path.join(self.output_dir, "svg")))
            # shutil.copytree(svg_path, os.path.join(self.output_dir, "svg"))
        except:
            self.set_error("svg文件夹移动失败")
        if os.path.exists(os.path.join(self.output_dir, "svg.tar.gz")):
            os.remove(os.path.join(self.output_dir, "svg.tar.gz"))
            # shutil.rmtree(os.path.join(self.output_dir, "svg.tar.gz"))
        try:
            # CopyFile().linkdir(svg_path, os.path.join(self.output_dir, "svg"))
            os.link(svg_gz_path, os.path.join(self.output_dir, "svg.tar.gz"))
            # os.system('ln -s {} {}'.format(svg_gz_path, os.path.join(self.output_dir, "svg.tar.gz")))
        except:
            self.set_error("svg文件夹移动失败")
        self.end()

    def replace_reactome_link(self):
        '''
        替换reactome链接
        '''
        enrich_result = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))[0]
        enrich_result_out = self.output_dir + '/' + os.path.basename(enrich_result)
        annot_result = self.reactome_class.output_dir + '/reactome_path.xls'
        annot_df = pd.read_table(annot_result, sep='\t',header=0)
        map2link = dict(zip(annot_df['Pathway ID'], annot_df['link']))
        map2des = dict(zip(annot_df['Pathway ID'], annot_df['Description']))
        map2category = dict(zip(annot_df['Pathway ID'], annot_df['category']))

        enrich_df = pd.read_table(enrich_result, sep='\t',header=0)
        enrich_df['Description'] = [map2des.get(x, "") for x in enrich_df['Pathway ID']]
        enrich_df['link'] = [map2link.get(x, "") for x in enrich_df['Pathway ID']]
        enrich_df['category'] = [map2category.get(x, "") for x in enrich_df['Pathway ID']]

        enrich_df.to_csv(enrich_result_out, sep = '\t', index=False)


    def run(self):
        super(DiffReactomeEnrichModule, self).run()
        self.logger.info("开始运行差异基因集kegg富集分析")
        self.run_reactome_rich()

    def end(self):
        super(DiffReactomeEnrichModule, self).end()



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "diff_geneset_enrich" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "medical_transcriptome.diff_geneset.diff_reactome_enrich",
            "instant": False,
            "options": dict(
                geneset_reactome="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/DiffGenesetPrepare4/output/kegg_class/multi_geneset_list",
                reactome_annot="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/AnnotPrepare/output/all_reactome.list",
                # kegg_table_2 = "/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/AnnotPrepare/output/gene_kegg_level_table.xls",
                geneset_list="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/DiffGenesetPrepare4/output/H1_vs_H3_all_12_gene.list",
                all_list="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/AnnotPrepare/output/all_gene.list",
                task_id ="medical_transcriptome",
                level = "G",
                # add_info="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/AnnotPrepare/output/add_info.txt",
                reactome_version ="202007"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()