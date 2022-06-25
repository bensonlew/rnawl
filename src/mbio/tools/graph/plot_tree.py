# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
# import copy
import collections


class PlotTreeAgent(Agent):
    """
    R
    version 1.0
    author shenghe
    last_modified:2016.10.18
    """

    def __init__(self, parent):
        super(PlotTreeAgent, self).__init__(parent)
        options = [
            {"name": "abundance_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "newicktree", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "leaves_group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "sample_group", "type": "infile", "format": "meta.otu.group_table"},
        ]
        self.add_option(options)
        self.step.add_steps('plot_tree')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.plot_tree.start()
        self.step.update()

    def step_end(self):
        self.step.plot_tree.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('abundance_table').is_set:
            raise OptionError('必须提供输入文件:物种丰富文件')
        else:
            otulist = [line.split('\t')[0] for line in open(self.option('abundance_table').prop['path'])][1:]  # 获取所有OTU/物种名
            sample = open(self.option('abundance_table').prop['path']).readline().rstrip().split('\t')[1:]
        if not self.option('newicktree').is_set:
            raise OptionError('必须提供输入文件:进化树')
        for i in otulist:
            if i not in self.option('newicktree').prop["sample"]:
                raise OptionError("丰度表中的物种在进化树中不存在:{}".format(i))
        if self.option("leaves_group").is_set:
            for i in self.option("leaves_group").prop['sample']:
                if i not in self.option("newicktree").prop['sample']:
                    raise OptionError("叶分组文件有叶信息不再进化树中")
        if self.option('sample_group').is_set:
            for i in sample:
                if i not in self.option("sample_group").prop['sample']:
                    raise OptionError("otu/物种丰度表中样本不在样本分组文件中")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 6
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "距离矩阵计算结果输出目录"],
            ["fan.png", "png", "环形树图结果文件"],
            ["bar.png", "png", "带有bar图的树结果文件"],
            ["fan.pdf", "pdf", "环形树图结果文件"],
            ["bar.pdf", "pdf", "带有bar图的树结果文件"]
        ])
        super(PlotTreeAgent, self).end()


class PlotTreeTool(Tool):

    def __init__(self, config):
        super(PlotTreeTool, self).__init__(config)
        self.plot_tree_bar_path = self.config.SOFTWARE_DIR + '/bioinfo/meta/scripts/plot-treeygbar.pl'
        self.plot_tree_fan_path = self.config.SOFTWARE_DIR + '/bioinfo/meta/scripts/plot-tree.pl'
        self.perl_path = 'miniconda2/bin/perl '
        self.colors_num = 21  # 目前支持的颜色数量
        self.max_species_group_name = 0

    def run(self):
        """
        运行
        """
        super(PlotTreeTool, self).run()
        self.run_plot_tree()

    def leaf_group_filter(self):
        """筛选leaf的颜色分类数量，当数量多余23个时去除多余的"""
        if not self.option('leaves_group').is_set:
            self.leaf_group_num = 0
            return None
        leaf_group_filter_file = self.work_dir + '/leaf_group_filter.temp'
        with open(self.option('leaves_group').path) as f, open(leaf_group_filter_file, 'w') as w:
            lines = f.readlines()
            headers = lines[0]
            lines = lines[1:]
            # lines_copy = copy.copy(lines)
            groups = [i.split('\t')[1] for i in lines]
            self.max_species_group_name = max([len(i) for i in groups])
            count_dict = collections.defaultdict(int)
            for i in groups:
                count_dict[i] = count_dict[i] + 1
            if len(count_dict) <= self.colors_num:
                self.leaf_group_num = len(count_dict)
                return self.option('leaves_group').path
            self.leaf_group_num = self.colors_num
            sorted_count = sorted(count_dict.items(), key=lambda x: x[1], reverse=True)
            sorted_count = sorted_count[: self.colors_num - 1]  # 此处本来不应该减 1 ，但是脚本存在问题，减一， 相当于少了一个颜色
            sorted_count = [i[0] for i in sorted_count]
            w.write(headers)
            for i in lines:
                if i.split('\t')[1] in sorted_count:
                    w.write(i)
        return leaf_group_filter_file


    def run_plot_tree(self):
        leaf_group_file = self.leaf_group_filter()
        cmd_fan = self.perl_path + self.plot_tree_fan_path + ' -i ' + self.option('newicktree').path + ' -o ' + self.output_dir + '/fan.pdf' + ' -tretype fan '
        if self.option('leaves_group').is_set:
            cmd_fan += ' -d ' + leaf_group_file
        cmd_bar = self.perl_path + self.plot_tree_bar_path + ' -i ' + self.option('newicktree').path + ' -o ' + self.output_dir + '/bar.pdf'
        cmd_bar += ' -t ' + self.option('abundance_table').path
        if self.option("sample_group").is_set:
            cmd_bar += ' -g ' + self.option('sample_group').path
            group_num = len(set([i.split('\t')[1] for i in open(self.option('sample_group').path).readlines()[1:]]))
        else:
            group_num = len(open(self.option('abundance_table').path).readline().split('\t'))
            cmd_bar += ' -g ALL'
        all_legend_num = group_num + self.leaf_group_num
        if all_legend_num > len(self.option('newicktree').prop['sample']):
            heigth = int(all_legend_num * 0.26)
        else:
            heigth = int(len(self.option('newicktree').prop['sample']) * 0.26)
        if heigth < 3:
            heigth = 3
        if heigth > 80:
            heigth = 80
        if self.option('leaves_group').is_set:
            cmd_bar += ' -d ' + leaf_group_file
        else:
            pass
            # cmd_bar += ' -d ALL'
        fan_label_size = 1
        if len(self.option('newicktree').prop['sample']) > 500:
            fan_label_size = 0.5
        elif len(self.option('newicktree').prop['sample']) > 100:
            fan_label_size = 1 - len(self.option('newicktree').prop['sample']) / 100.0 / 10.0
        else:
            fan_label_size = 1
        max_name_length = 0
        for i in self.option('newicktree').prop['sample']:
            name_length = len(i)
            if name_length > max_name_length:
                max_name_length = name_length
        print max_name_length
        width = int(max_name_length * 0.4)
        if width < 10:
            width = 10
        width = int(heigth / 10) + width
        if self.max_species_group_name > 9:
            lw = '6-1-2'
        else:
            lw = '6-1-1'
        cmd_bar += ' -h {} -w {} -lw {} -cex 1.2'.format(heigth, width, lw)
        # fan_size = int(width * 0.7)
        # if fan_size < 8:
        #     fan_size = 8
        # temp_size = heigth / 2
        fan_size = int(heigth * 0.6) + int(width * 0.65)
        cmd_fan += ' -h {} -w {} -cex {}'.format(fan_size, fan_size, fan_label_size)
        fan_command = self.add_command('fan', cmd_fan)
        bar_command = self.add_command('bar', cmd_bar)
        fan_command.run()
        bar_command.run()
        self.wait()
        if bar_command.return_code == 0 and fan_command.return_code == 0:
            if os.path.isfile(self.output_dir + '/fan.pdf') and os.path.isfile(self.output_dir + '/bar.pdf'):
                if self.convert_pdf_to_png([self.output_dir + '/fan.pdf', self.output_dir + '/bar.pdf'],
                                           [self.output_dir + '/fan.png', self.output_dir + '/bar.png']):
                    self.end()
                else:
                    self.set_error("结果格式转换发生错误")
            else:
                self.set_error('结果文件没有正确生成')
        else:
            self.logger.info(bar_command.return_code)
            self.logger.info(fan_command.return_code)
            self.set_error('程序计算错误')

    def convert_pdf_to_png(self, olds, news):
        self.image_magick = '/program/ImageMagick/bin/convert'
        convert_commands = []
        for index, i in enumerate(olds):
            cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + i + ' ' + news[index]
            command = self.add_command('convert_{}'.format(index), cmd)
            command.run()
            convert_commands.append(command)
        self.wait()
        for i in convert_commands:
            if i.return_code == 0:
                pass
            else:
                return False
        return True
