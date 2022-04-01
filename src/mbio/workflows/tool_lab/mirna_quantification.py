# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
import re
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
from biocluster.config import Config


class MirnaQuantificationWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MirnaQuantificationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "category", "type": "string", "default": ""},  # 物种分类，Animal or Plant
            {"name": "input_fa", "type": "infile", "format": "small_rna.fasta"},  # 输入miRNA序列文件
            {"name": "database", "type": "string", "default": "mirbase"},  # 参考miRNA数据库
            {"name": "org_choose", "type": "string", "default": ""},  # 物种选择，[all, custom, auto]
            {"name": "species", "type": "string", "default": ""},  # 具体物种,多选时分号分隔
            {"name": "species_1", "type": "string", "default": ""},  # 当database为miRBase时
            {"name": "species_2", "type": "string", "default": ""},  # 当database为PmiREN时
            {"name": "organism_list", "type": "infile", 'format': "small_rna.common"},  # known mirna鉴定的物种列表
            {"name": "mismatch", "type": "int", "default": 0},  # 允许错误匹配
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "dir", "type": "string", "default": "dir"},  # 输出目录后缀
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.known_mirna")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(MirnaQuantificationWorkflow, self).run()

    def check_options(self):
        if self.option('org_choose').lower() not in ['all', 'custom', 'auto']:
            raise OptionError('请输入物种选择参数')
        if self.option('database').lower() == 'mirbase':
            self.option('species', self.option("species_1"))
        else:
            self.option('species', self.option("species_2"))
        if self.option("org_choose").lower() == "auto":
            if not self.option("organism_list").is_set:
                raise OptionError("必须指定物种列表")
        elif self.option("org_choose").lower() == "all":
            self.option('species', 'all')
        else:
            if self.option("species") == "":
                raise OptionError("物种列表不能为空")
        if self.option('database').lower() not in ['mirbase', 'pmiren']:
            raise OptionError('数据库只能输入miRBase/PmiREN')
        if self.option('category').lower() not in ['plant', 'animal']:
            raise OptionError('物种类别只支持plant/animal')
        elif self.option('category').lower() == 'animal' and self.option('database').lower() == 'pmiren':
            raise OptionError('当物种类别为animal时，只能选择miRBase数据库')
        if not self.option("input_fa").is_set:
            raise OptionError("必须输入miRNA序列文件")
        with open(self.option("input_fa").prop["path"], "r") as f:
            line = f.readline()
            if not re.compile(r'^>\S+_x\d+').match(line):
                raise OptionError("输入miRNA序列文件格式不对")
        return True

    def run_tool(self):
        if self.option("org_choose").lower() == "all":
            if self.option("category").lower() == 'animal':
                species = "all_animal"
            else:
                species = "all_plant"
            self.option("species", species)
        elif self.option("org_choose").lower() == "custom":
            species = self.option("species").split(",")
            species_list = list()
            for i in species:
                i = i.strip()
                if len(i) >= 3:
                    if self.option("database").lower() == "pmiren":
                        pmiren = \
                            Config().get_mongo_client(mtype='small_rna')[Config().get_mongo_dbname('small_rna')][
                                'pmiren']
                        try:
                            organism = pmiren.find_one({"Name": i})["organism"]
                            i = organism
                        except:
                            pass
                    else:
                        mirbase = \
                            Config().get_mongo_client(mtype='small_rna')[Config().get_mongo_dbname('small_rna')][
                                'mirbase']
                        try:
                            organism = mirbase.find_one({"Name": i})["organism"]
                            i = organism
                        except:
                            pass
                if i not in species_list:
                    species_list.append(i)
            self.option("species", ",".join(species_list))
        options = {
            "input_fa": self.option("input_fa"),
            "mismatch": self.option("mismatch"),
            "database": self.option("database"),
            "org_choose": self.option("org_choose"),
            "category": self.option("category"),
        }
        if self.option("org_choose").lower() == "auto":
            options.update({'organism_list': self.option("organism_list")})
        else:
            options.update({'species': self.option("species")})
        if self.option("config").is_set:
            options.update({'config': self.option("config")})
        self.tool.set_options(options)
        self.tool.on('end', self.set_output)
        self.tool.run()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "miRNA定量结果目录", 0],
            ["known_mirna_count.xls", "", "miRNA定量表(count)", 0],
            ["known_mirna_detail.xls", "", "miRNA定量详情表", 0],
            ["known_mirna_norm.xls", "", "miRNA定量表(normalized)", 0],
        ])
        super(MirnaQuantificationWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.mirna_quantification import MirnaQuantificationWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            "id": "mirna_quantification_" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.mirna_quantification",
            "options": dict(
                # org_choose='auto',
                # organism_list="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/mirdp2_test/ath_test/orgnism_list",
                # category='plant',
                # database='pmiren',
                # input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/quantifier_test/rfam_trimed.fa",
                org_choose='custom',
                species="dme,hsa,mmu",
                category='animal',
                database='mirbase',
                input_fa="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/quantifier_test/rfam_trimed.fa",
            )
        }
        wsheet = Sheet(data=data)
        wf =MirnaQuantificationWorkflow(wsheet)
        wf.sheet.id = 'mirna_quantification'
        wf.sheet.project_sn = 'mirna_quantification'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
