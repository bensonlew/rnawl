# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError


class GwasWorkflow(Workflow):
    """
    GWAS
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GwasWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": 'trait_file', 'type': 'infile', 'format': "ref_rna_v2.common"},
            {"name": "max_missing", "type": "string", 'default': '30%'},  # --max-missing
            {"name": "maf", "type": "string", 'default': '0.05'},  # --maf 0.05
            {"name": "method", "type": "string", 'default': 'mlm'},
            {"name": 'alpha', "type": "string", 'default': '0.05'},
            {"name": 'advanced_params', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.filter = self.add_tool("tool_lab.gwas.vcftools_filter")
        self.plink = self.add_tool("tool_lab.gwas.plink")
        self.gwas = self.add_tool("tool_lab.gwas.gwas")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_filter()
        super(GwasWorkflow, self).run()

    def check_options(self):
        if not self.option("vcf_path").is_set:
            raise OptionError("必须设置输入VCF文件")
        if not self.option('trait_file').is_set:
            raise OptionError("请设置输入选择性状文件")
        return True

    def run_filter(self):
        missing = self.option('max_missing').split('%')[0]
        opts = {
            'vcf_path': self.option('vcf_path').prop['path'],
            'max_missing': float(missing)/100,
            'maf': float(self.option('maf')),
        }
        self.filter.set_options(opts)
        self.filter.on('end', self.run_plink)
        self.filter.run()

    def run_plink(self):
        opts = {
            'vcf_file': self.filter.option('pop_recode').prop['path'],
            'chrom_map': self.filter.option('chrom_map').prop['path']
        }
        self.plink.set_options(opts)
        self.plink.on('end', self.run_gwas)
        self.plink.run()

    def run_gwas(self):
        if self.option('method').lower() == 'glm':
            method = 'GLM'
        elif self.option('method').lower() == 'farmcpu':
            method = 'FarmCPU'
        else:
            method = 'MLM'
        opts = {
            'plink': self.plink.output_dir,
            'trait': self.option('trait_file').prop['path'],
            'method': method,
            'alpha': float(self.option('alpha')),
            'chrom_map': self.filter.option('chrom_map').prop['path']
        }
        self.gwas.set_options(opts)
        self.gwas.on('end', self.set_output)
        self.gwas.run()

    def set_output(self):
        files = glob.glob(os.path.join(self.gwas.output_dir, '*.pdf'))
        files += glob.glob(os.path.join(self.gwas.output_dir, '*.xls'))
        for file in files:
            if file.endswith('.pdf'):
                file_new = os.path.join(self.output_dir, '.'.join(os.path.basename(file).split('.')[:-2]) + '.pdf')
            else:
                file_new = os.path.join(self.output_dir, os.path.basename(file))
            if os.path.exists(file_new):
                os.remove(file_new)
            os.link(file, file_new)
        self.set_db()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        gwas_api = self.api.api('tool_lab.gwas')
        pngs = glob.glob(os.path.join(self.output_dir, '*.pdf'))
        pic_data = dict()
        for each in pngs:
            png = os.path.join(self._sheet.output, os.path.basename(each))
            pic_type = os.path.basename(each).split('.')
            trait = pic_type[:-3]
            if len(trait) > 1:
                trait = '_'.join(trait)
            else:
                trait = trait[0]
            if trait not in pic_data.keys():
                pic_data[trait] = dict()
            if 'Rectangular-Manhattan' in pic_type:
                pic_data[trait]['manhattan'] = png
            elif 'QQplot' in pic_type:
                pic_data[trait]['qq'] = png
        result = glob.glob(os.path.join(self.output_dir, '*.xls'))
        gwas_api.add_gwas_main(result=result, pic_data=pic_data, trait=pic_data.keys(), main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "GWAS分析结果文件",0],
            ['*xls', 'xls', '差异位点表', 0],
            ['*png', 'png', 'GWAS分析结果图', 0],
            ['*pdf', 'pdf', 'GWAS分析结果图', 0],
        ])
        super(GwasWorkflow, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the workflow. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitoollabtest.py '
        cmd += 'post toollabpipeline '
        cmd += '-c {} '.format("client03")
        cmd += "-b http://bcl.tsg.com "
        cmd += "-n \"params;basis\" -d \"{"
        args = dict(
            vcf_path='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/gwas/test_short.vcf',
            trait_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/gwas/test_trait1.txt'
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.gwas",
            main_table_name="sg_gwas",
            task_id="gwas",
            project_sn="gwas",
            submit_location="gwas"
        )
        for arg in args:
            cmd += "\\\""
            cmd += arg
            cmd += "\\\":\\\""
            cmd += args[arg]
            cmd += "\\\","
        cmd = cmd.rstrip(",")
        cmd += "};{"
        for arg in config:
            cmd += "\\\""
            cmd += arg
            cmd += "\\\":\\\""
            cmd += config[arg]
            cmd += "\\\","
        cmd = cmd.rstrip(",")
        cmd += "}\""

        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()