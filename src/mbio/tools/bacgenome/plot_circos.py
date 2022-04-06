# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import subprocess
import cairosvg

class PlotCircosAgent(Agent):
    """

    author: zouxuan
    last_modify: 20180402
    """

    def __init__(self, parent):
        super(PlotCircosAgent, self).__init__(parent)
        options = [
            {"name": "k", "type": "infile", "format": "sequence.profile_table"},  # karyotype.txt
            {"name": "c", "type": "infile", "format": "sequence.profile_table"},  # sense_strand_cog.txt
            {"name": "t", "type": "infile", "format": "sequence.profile_table"},  # temp.txt
            {"name": "ac", "type": "infile", "format": "sequence.profile_table"},  # antisense_strand_cog.txt
            {"name": "pgc", "type": "infile", "format": "sequence.profile_table"},  # positive_gc_count.txt
            {"name": "ngc", "type": "infile", "format": "sequence.profile_table"},  # negative_gc_count.txt
            {"name": "pgs", "type": "infile", "format": "sequence.profile_table"},  # positive_gc_skew.txt
            {"name": "ngs", "type": "infile", "format": "sequence.profile_table"},  # negative_gc_skew.txt
            {"name": "type", "type": "int", "default": 1},  # 类型调节画图参数，质粒为2，其余为1
            {"name": "f", "type": "string", "default": ""},  # 文件路径以逗号分隔
            {"name": "labs", "type": "string", "default": ""}, #分号分割
            {"name": "seq_type", "type": "string", "default": "Circular"}, # 开闭环
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("k").is_set:
            raise OptionError("必须设置karyotype文件", code="31402601")
        if not self.option("c").is_set:
            raise OptionError("必须设置sense_strand_cog文件", code="31402602")
        if not self.option("t").is_set:
            raise OptionError("必须设置temp文件", code="31402603")
        if not self.option("ac").is_set:
            raise OptionError("必须设置antisense_strand_cog文件", code="31402604")
        if not self.option("pgc").is_set:
            raise OptionError("positive_gc_count文件", code="31402605")
        if not self.option("ngc").is_set:
            raise OptionError("必须设置negative_gc_count文件", code="31402606")
        if not self.option("pgs").is_set:
            raise OptionError("必须设置positive_gc_skew文件", code="31402607")
        if not self.option("ngs").is_set:
            raise OptionError("必须设置negative_gc_skew文件", code="31402608")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(PlotCircosAgent, self).end()


class PlotCircosTool(Tool):
    def __init__(self, config):
        super(PlotCircosTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.config_path = self.config.PACKAGE_DIR + '/bacgenome/circos_confi.py'
        self.circos_path = self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/circos-0.69-6/bin/circos'
        self.python_path = '/miniconda2/bin/python'
        self.legend_py = self.config.PACKAGE_DIR + '/bacgenome/circos_legend.py'

    def run(self):
        """
        运行
        :return:
        """
        super(PlotCircosTool, self).run()
        self.run_circos()
        self.end()

    def run_circos(self):
        if self.option('f'):
            file = "-f " + self.option('f')
        else:
            file = ""
        cmd = '{} {} -k {} -c {} -t {} -ac {} -pgc {} -ngc {} -pgs {} -ngs {} {} -o {} -type {} -ci {}'. \
            format(self.python_path, self.config_path, self.option('k').prop['path'], self.option('c').prop['path'],
                    self.option('t').prop['path'], self.option('ac').prop['path'], self.option('pgc').prop['path'],
                    self.option('ngc').prop['path'], self.option('pgs').prop['path'],
                    self.option('ngs').prop['path'],
                    file, self.work_dir + '/circos.conf', str(self.option('type')), str(self.option('seq_type')))
        command = self.add_command('config_file', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("config_file生成成功")
        else:
            self.set_error("config_file生成失败", code="31402601")
            self.set_error("config_file生成失败", code="31402602")
        cmd1 = '{} {} -conf {} -png -svg'.format(self.perl_path, self.circos_path, self.work_dir + '/circos.conf')   #-param image/radius=500p
        command1 = self.add_command('plot_circos', cmd1).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("circos图片生成成功")
        else:
            self.set_error("circos图片生成失败", code="31402603")
            self.set_error("circos图片生成失败", code="31402604")

        ## 生成图例  zouguanqing 20190410
        if self.option('labs')!='':
            trans = {
                'gi':'island',
                'prephage':'prephage',
                '5S_rRNA':'rna',
                '16S_rRNA':'rna',
                '23S_rRNA':'rna',
                'tRNA':'rna',
                'ncrna':'rna',
                'cog' : 'COG',
                'is': 'isf',
                'integron': 'integron',
                'para_def':'def_circle',
            }
            labs_old = self.option('labs').split(';')
            #labs = 'COG'   #COG在第一圈
            labs_new_list = []
            for each in labs_old:
                tmp = []
                if each:   # by zzg 20201209 特定区段为空的情况
                    for each1 in each.split(','):
                        if trans[each1] not in tmp:
                            tmp.append(trans[each1])
                    labs_new_list.append(','.join(tmp))
            labs = ';'.join(labs_new_list)

            cmd2 = '{} {} -i {} -o legend.svg'.format(self.python_path, self.legend_py,labs)
            command2 = self.add_command('plot_legend', cmd2).run()
            self.wait(command2)
            if command2.return_code == 0:
                self.logger.info("circos图例生成成功")
            else:
                self.set_error("circos图例生成失败")
            ## 圈图的svg和图例svg合并,并由svg生成png
            self.merge_svg(self.work_dir + '/circos.svg', self.work_dir + '/legend.svg')
            self.convert_svg_to(self.work_dir + '/circos.svg',self.work_dir + '/circos.png')

        if os.path.exists(self.output_dir + '/circos.png'):
            os.remove(self.output_dir + '/circos.png')
        os.link(self.work_dir + '/circos.png', self.output_dir + '/circos.png')
        if os.path.exists(self.output_dir + '/circos.svg'):
            os.remove(self.output_dir + '/circos.svg')
        os.link(self.work_dir + '/circos.svg', self.output_dir + '/circos.svg')

    #zouguanqing 20190410
    def merge_svg(self, circos_svg,legend):  #legend 的svg 有1行或2行
        legend_pk = ''
        with open(legend) as f:
            for line in f:
                if re.match('<\?xml version=', line):
                    continue
                if re.match('<svg',line):
                    legend_pk = re.sub('<svg\s*[^>]*>','',line)

        with open(circos_svg) as f2:
            tmp = f2.read()
            tmp = tmp.replace('<svg width="1000px"','<svg width="1500px"')
            tmp = tmp.replace('<rect x="0" y="0" width="1000px"','<rect x="0" y="0" width="1500px"')
            tmp = tmp.replace('</svg>', legend_pk)

        with open('circos_new.svg','w') as fw:
            fw.write(tmp)
        os.rename('circos.svg','circos_no_legend.svg')
        os.rename('circos_new.svg','circos.svg')

    #zouguanqing 20190410
    def convert_svg_to(self,old,new):
        #self.image_magick = '/program/ImageMagick/bin/convert'
        m = {'png':new}
        for i in ['png']:  #for循环，便于添加其他格式的转化
            #cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + old + ' ' + m[i]
            cairosvg.svg2png(url='circos.svg',write_to=m[i])
            # try:
            #     subprocess.check_output(self.config.SOFTWARE_DIR + '/'+cmd,shell=True)
            # except:
            #     self.logger.info('%s convert failed,maybe has result'%cmd)


