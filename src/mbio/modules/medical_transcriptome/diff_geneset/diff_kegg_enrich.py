# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20200915

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import glob
import re
import os
import unittest
import gevent.subprocess as subprocess
import shutil
import json
from collections import OrderedDict
from mbio.packages.medical_transcriptome.copy_file import CopyFile
from biocluster.config import Config


class DiffKeggEnrichModule(Module):
    """
    该Module用于基因融合分析，默认使用方法
    """
    def __init__(self, work_id):
        super(DiffKeggEnrichModule, self).__init__(work_id)
        options = [
            {"name": "task_id", "type": "string", "default": None},
            {"name": "regulate", "type": "string", "default": None},
            {'name': 'geneset_kegg', 'type': 'string'},
            {'name': 'kegg_table', 'type': 'string'},
            {'name': 'geneset_list', 'type': 'string'},
            {'name': 'all_list', 'type': 'string'},
            {'name': 'regulate', 'type': 'string', 'default': None},
            {"name": "level", "type": "string", "default": None},
            {'name': 'add_info', 'type': 'string', 'default': None},
            {'name': 'kegg_table2', 'type': 'string'},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
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

    def run_kegg_rich(self):
        self.enrich_tool = self.add_tool('medical_transcriptome.diff_geneset.diff_kegg_rich')
        self.enrich_tool.set_options({
            'kegg_table': self.option('kegg_table'),
            'diff_list': self.option('geneset_list'),
            'kegg_table2': self.option('kegg_table2'),
            'add_info': self.option('add_info'),
            "kegg_version": self.option('kegg_version'),
        })
        self.enrich_tool.on('end', self.run_kegg_class)
        self.enrich_tool.run()

    def run_kegg_class(self):
        self.kegg_class = self.add_tool('medical_transcriptome.diff_geneset.diff_kegg_class')
        self.kegg_class.set_options({
            'geneset_kegg': self.option('geneset_kegg'),
            'kegg_table': self.option('kegg_table'),
            'kegg_table2': self.option('kegg_table2'),
            "regulate":self.option("regulate"),
            'level': self.option('level'),
            'background_links': self.option('add_info'),
            'task_id': self.option('task_id'),
            "kegg_version": self.option('kegg_version'),
        })
        self.kegg_class.on('end', self.set_output)
        self.kegg_class.run()

    def set_output(self):
        self.replace_kegg_link()
        output_dir1 = self.enrich_tool.output_dir
        output_dir2 = self.kegg_class.output_dir
        if os.path.exists(os.path.join(self.output_dir,"enrich")):
            shutil.rmtree(os.path.join(self.output_dir,"enrich"))
        if os.path.exists(os.path.join(self.output_dir,"class")):
            shutil.rmtree(os.path.join(self.output_dir,"class"))
        try:
            # os.system('ln -s {} {}'.format(output_dir1, os.path.join(self.output_dir,"enrich")))
            # os.system('ln -s {} {}'.format(output_dir2, os.path.join(self.output_dir, "class")))

            # 20210125修改
            CopyFile().linkdir(output_dir1, os.path.join(self.output_dir, "enrich"))
            CopyFile().linkdir(output_dir2, os.path.join(self.output_dir, "class"))
            # shutil.copytree(output_dir1,os.path.join(self.output_dir,"enrich"))
            # shutil.copytree(output_dir2, os.path.join(self.output_dir, "class"))
        except:
            raise Exception("移动kegg富集结果文件夹失败")

        self.end()


    def run(self):
        super(DiffKeggEnrichModule, self).run()
        self.logger.info("开始运行差异基因集kegg富集分析")
        self.run_kegg_rich()

    def replace_kegg_link(self):
        '''
        替换kegg链接
        '''
        enrich_result = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))[0]
        annot_result = self.kegg_class.output_dir + '/kegg_stat.xls'
        with open(annot_result, 'rb') as f:
            map2link = {line.split("\t")[0]:line.strip().split("\t")[-1] for line in f.readlines()[1:] if line.split("\t")[-1].startswith("http")}
        if os.path.exists(enrich_result + "relink"):
            os.remove(enrich_result + "relink")
        with open(enrich_result, 'rb') as f, open(enrich_result + "relink", 'w') as fo:
            header = f.readline()
            fo.write(header)
            for line in f:
                cols = line.split("\t")
                if cols[3] in map2link:
                    cols[9] = map2link[cols[3]]
                fo.write("\t".join(cols))
        os.remove(enrich_result)
        os.link(enrich_result + "relink", enrich_result)


    def end(self):
        super(DiffKeggEnrichModule, self).end()



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
            "name": "medical_transcriptome.diff_geneset.diff_kegg_enrich",
            "instant": False,
            "options": dict(
                geneset_kegg="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/DiffGenesetPrepare4/output/kegg_class/multi_geneset_list",
                kegg_table="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/AnnotPrepare/output/gene_kegg_table.xls",
                kegg_table_2 = "/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/AnnotPrepare/output/gene_kegg_level_table.xls",
                geneset_list="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/DiffGenesetPrepare4/output/H1_vs_H3_all_12_gene.list",
                all_list="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/AnnotPrepare/output/all_gene.list",
                task_id ="medical_transcriptome",
                level = "G",
                add_info="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare6287/DiffGenesetPrepare/AnnotPrepare/output/add_info.txt",
                kegg_version ="202007"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()