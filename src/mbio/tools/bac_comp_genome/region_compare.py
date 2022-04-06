# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class RegionCompareAgent(Agent):
    def __init__(self, parent):
        super(RegionCompareAgent, self).__init__(parent)
        options = [
            {'name': 'seq_dir', 'type': 'string'},
            {'name': 'samples', 'type': 'string'},  # 样本序列文件路径，同逗号分隔，参考序列在第一个
            {'name': 'region', 'type': 'string'},  # 选择的参考基因组的区域，
        ]
        self.add_option(options)
        self.step.add_steps('region_comp')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.region_comp.start()
        self.step.update()

    def stepfinish(self):
        self.step.region_comp.finish()
        self.step.update()

    def check_options(self):
        if not os.path.isdir(self.option('seq_dir')):
            raise OptionError('请正确给出序列文件夹，{}')
        if set(self.option('samples').split(',')) - set(os.listdir(self.option('seq_dir'))):
            raise OptionError('samples参数中存在样本没有序列，请确定')
        if ',' not in self.option('region'):
            raise OptionError('参数 region 必须以逗号分隔区域')
        elif len(self.option('region').split(',')) != 3:
            raise OptionError('参数region的格式为"chr,start,end"')

    def set_resource(self):
        self._cpu = len(self.option('samples').split(',')) + 1
        self._memory = str(self._cpu) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '', 'region compare结果目录'],
            ['regions.xls', '', '基因组间的同源性区域'],
            ['mash_out.xls', '', '参考基因组和其它基因组的距离关系']
        ])
        super(RegionCompareAgent, self).end()


class RegionCompareTool(Tool):
    def __init__(self, config):
        super(RegionCompareTool, self).__init__(config)
        self.python = '/miniconda2/bin/python'
        self.mum_path = self.config.SOFTWARE_DIR +\
            '/bioinfo/compare_genome/software/MUMmer3.23'
        self.mash_path = self.config.SOFTWARE_DIR +\
            '/bioinfo/compare_genome/software/mash/mash'
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/region_compare.py'

    def run(self):
        super(RegionCompareTool, self).run()
        self.region_comp_run()
        self.set_output()
        self.end()

    def region_comp_run(self):
        samples = self.comb_genome()
        cmd = self.python + ' ' + self.package_path +\
            ' -s {} -r {} -m {} -n {}'.format(samples, self.option('region'),
                                              self.mum_path, self.mash_path)
        command = self.add_command('region_comp', cmd).run()
        self.wait(command)
        map(os.remove, samples.split(','))

        if command.return_code == 0:
            self.logger.info('package region_compare.py 运行成功')
        else:
            self.set_error('package region_compare.py 运行出错')

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for f in files:
                os.remove(os.path.join(root, f))
        os.link('mash_out.xls', self.output_dir + '/mash_out.xls')
        os.link('regions.xls', self.output_dir + '/regions.xls')

    def end(self):
        super(RegionCompareTool, self).end()

    def comb_genome(self):
        sp_list = self.option('samples').split(',')
        samples = []
        i = 0
        for f in sp_list:
            fasta_dir = os.path.join(self.option('seq_dir'), f)
            cmd = '/bin/cat {}/* > {}.fna'.format(fasta_dir, f)
            samples.append(f + '.fna')
            self.add_command('region_compre-comb_genome_' + str(i), cmd, shell=True).run()
            i += 1
        self.wait()
        return ','.join(samples)
