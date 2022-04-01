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


class BsaWorkflow(Workflow):
    """
    carry out selective sweep analysis  using VCFtools
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BsaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "max_missing", "type": "string", 'default': '30%'},  # --max-missing
            {"name": "maf", "type": "string", 'default': '0.05'},  # --maf 0.05
            {"name": "win_size", "type": "string", "default": '2'},  # --window-pi,窗口大小
            {"name": "win_step", "type": "string", "default": '0.1'},  # --window-pi-step,窗口步长
            {"name": "wp", "type": "string"},  # wild parent
            {"name": "mp", "type": "string"},  # mutant parent
            {"name": "wb", "type": "string"},  # wild bulk
            {"name": "mb", "type": "string"},  # mutant bulk
            {"name": "method", "type": "string", 'default': 'variants_index'},
            {"name": 'advanced_params', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.filter = self.add_tool("tool_lab.bsa.vcftools_filter")
        self.bsa_calc = self.add_tool("tool_lab.bsa.bsa_calc")
        self.model_con = self.add_tool("tool_lab.bsa.model_construction")
        self.manhattan = self.add_tool("tool_lab.bsa.manhattan")
        self.region_calc = self.add_tool("tool_lab.bsa.region_calc")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_filter()
        super(BsaWorkflow, self).run()

    def check_options(self):
        if not self.option("vcf_path").is_set:
            raise OptionError("必须设置输入VCF文件")
        if not self.option('wb'):
            raise OptionError("请设置输入野生型混池信息")
        if not self.option('mb'):
            raise OptionError("请设置输入突变型混池信息")
        # check sample names
        input_list = list()
        input_list.extend([self.option('wb'), self.option('mb')])
        if self.option('wp'):
            input_list.append(self.option('wp'))
        if self.option('mp'):
            input_list.append(self.option('mp'))
        with open(self.option('vcf_path').prop['path'], 'r') as vcf:
            for line in vcf:
                if line.startswith('#CHROM'):
                    samples = line.strip().split('\t')[8:]
                    break
        for i in input_list:
            if i not in samples:
                raise OptionError('选择的野生型或者突变型与输入的vcf文件不符，请检查。')
        return True

    def run_filter(self):
        missing = self.option('max_missing').split('%')[0]
        opts = {
            'vcf_path': self.option('vcf_path').prop['path'],
            'max_missing': float(missing)/100,
            'maf': float(self.option('maf')),
        }
        self.filter.set_options(opts)
        self.filter.on('end', self.calc_bsa)
        self.filter.run()

    def calc_bsa(self):
        opts = {
            'vcf_file': self.filter.option('filtered_vcf').prop['path'],
            'wb': self.option('wb'),
            'mb': self.option('mb'),
        }
        if self.option('wp') and self.option('mp'):
            opts.update({'wp': self.option('wp'), 'mp': self.option('mp')})
        self.bsa_calc.set_options(opts)
        self.bsa_calc.on('end', self.vcf2table_check)
        self.bsa_calc.run()

    def vcf2table_check(self):
        if self.bsa_calc.option('pop_index').is_set:
            self.construct_model()
        else:
            self.set_db()

    def construct_model(self):
        opts = {
            'pop_index': self.bsa_calc.option('pop_index').prop['path'],
            'method': self.option('method'),
            'win_size': int(float(self.option('win_size'))*10**6),
            'win_step': int(float(self.option('win_step'))*10**6),
        }
        self.model_con.set_options(opts)
        self.model_con.on('end', self.draw_manhattan)
        self.model_con.on('end', self.calc_region)
        self.on_rely([self.manhattan, self.region_calc], self.set_output)
        self.model_con.run()

    def draw_manhattan(self):
        opts = {
            'model_file': self.model_con.option('model_file').prop['path'],
            'method': self.option('method'),
        }
        self.manhattan.set_options(opts)
        self.manhattan.run()

    def calc_region(self):
        opts = {
            'model_file': self.model_con.option('model_file').prop['path'],
            'method': self.option('method'),
        }
        self.region_calc.set_options(opts)
        self.region_calc.run()

    def set_output(self):
        files = glob.glob(os.path.join(self.manhattan.output_dir, 'pop.chr.*'))
        files += glob.glob(os.path.join(self.manhattan.output_dir, 'pop.G.*'))
        files += glob.glob(os.path.join(self.region_calc.output_dir, '*'))
        for each in files:
            fname = os.path.basename(each)
            out = os.path.join(self.output_dir, fname)
            if os.path.exists(out):
                os.remove(out)
            os.link(each, out)
        self.set_db()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        bsa_api = self.api.api('tool_lab.bsa')
        result = glob.glob(os.path.join(self.region_calc.output_dir, 'pop.result.*'))
        if len(result) == 0:
            bsa_api.add_bsa_main(main_id=self.option('main_id'))
        else:
            png = glob.glob(os.path.join(self.output_dir, '*.png'))[0]
            png = os.path.join(self._sheet.output, os.path.basename(png))
            bsa_api.add_bsa_main(result=result[0], png=png, main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "基于VCF进行BSA分析",0],
            ['./pop.result*', '', '关联区域变异位点统计表', 0],
            ['./*png', 'png', 'Manhattan图', 0],
            ['./*pdf', 'pdf', 'Manhattan图', 0],
        ])
        super(BsaWorkflow, self).end()


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
            vcf_path='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/bsa/pop.final.vcf',
            # wp='',
            # mp='HQS1',
            mb='XS11_1',
            wb='F44_mix',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.bsa",
            main_table_name="sg_bsa",
            task_id="bsa",
            project_sn="bsa",
            submit_location="bsa"
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