# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class VcfTreeWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VcfTreeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "group_file","type": "infile","format":"denovo_rna_v2.common"},
            {"name": "max_missing", "type": "float", "default": 0.3},  # 缺失率
            {"name": "maf", "type": "float", "default": 0.05},  # 次要等位基因频率min
            {"name": "bs_trees", "type": "int", "default": 1000},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.filter = self.add_tool("tool_lab.vcf_filter")
        # self.plink = self.add_tool("tool_lab.vcf_plink")
        self.tree = self.add_module('dna_evolution.tree_generic')

    def check_option(self):
        """
        参数检查
        """
        if not self.option("vcf_file"):
            raise OptionError("必须输入vcf文件")
        if not self.option("group_file"):
            raise OptionError("必须输入分组文件")
        if not self.option("max_missing"):
            raise OptionError("必须设置最大缺失率")
        if not self.option("maf"):
            raise OptionError("必须设置maf")
    
    def run_filter(self):
        self.filter.set_options({
            'vcf_file':self.option("vcf_file"),
            'max_missing':self.option("max_missing"),
            "maf":self.option('maf')
        })
        self.filter.on('end',self.run_tree)
        self.filter.run()

    # def run_plink(self):
    #     self.filter_vcf = os.path.join(self.filter.output_dir,"pop.recode.vcf")
    #     self.plink.set_options({
    #         'vcf_file': self.filter_vcf,
    #     })
    #     self.plink.on('end',self.run_tree)
    #     self.plink.run()

    def run_tree(self):
        self.filter_vcf = os.path.join(self.filter.output_dir,"pop.recode.vcf")
        self.tree.set_options({
            'recode_vcf_path':self.filter_vcf,
            'tree_type': "ml",
            'bs_trees': self.option("bs_trees")
        })
        self.tree.on('end',self.set_output)
        self.tree.run()
    
    def set_output(self):
        self.linkdir(self.tree.output_dir, 'tree')
        self.set_db()
    
    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        self.logger.info("开始导表")
        tree_path = self.output_dir + "/tree/tree/pop.phylip.raxml.bestTree"
        api_vcf_tree = self.api.api("tool_lab.vcf_tree")
        api_vcf_tree.add_tree_detail(self.option("main_id"),tree_path,self.option("group_file").prop["path"])
        self.end()

    def run(self):
        self.run_filter()
        super(VcfTreeWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(VcfTreeWorkflow, self).end()

if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    import datetime
    add_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f") +\
        str(random.randint(1000,10000))
    data = {
        'name': 'test_vcf_tree',
        'id': 'vcf_tree_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        'options': {
        "vcf_file": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/test_WGS/5.tree/pop11.vcf.recode.vcf",
        "group_file":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/test_WGS/5.tree/group.list",
        "bs_trees": 2,
        "main_id" : add_time
        }
    }
    wsheet = Sheet(data=data)
    wf = VcfTreeWorkflow(wsheet)
    wf.run()