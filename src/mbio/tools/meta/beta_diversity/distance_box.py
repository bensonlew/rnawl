# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from biocluster.core.exceptions import OptionError


class DistanceBoxAgent(Agent):
    """
    qiime
    version 1.0
    author shenghe
    last_modified:2016.3.24
    """

    def __init__(self, parent):
        super(DistanceBoxAgent, self).__init__(parent)
        options = [
            {"name": "dis_matrix", "type": "infile",
                "format": "meta.beta_diversity.distance_matrix"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "grouplab", "type": "string", "default": ""},
            {"name": "permutations", "type": "int", "default": 999},
        ]
        self.add_option(options)
        self.step.add_steps('distancebox')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.distancebox.start()
        self.step.update()

    def step_end(self):
        self.step.distancebox.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        samplelist = []
        if 10000 >= self.option('permutations') >= 10:
            pass
        else:
            raise OptionError('随机置换次数:%s不在正常范围内[10, 10000]' % self.option('permutations'))
        if not self.option('dis_matrix').is_set:
            raise OptionError('必须提供距离矩阵文件')
        else:
            self.option('dis_matrix').get_info()
            samplelist = self.option('dis_matrix').prop['samp_list']
        if not self.option('group').is_set:
            raise OptionError('必须提供分组信息文件')
        else:
            self.option('group').get_info()
            if self.option('grouplab'):
                if self.option('grouplab') not in self.option('group').prop['group_scheme']:
                    raise OptionError('选定的分组方案名:%s在分组文件中不存在' % self.option('grouplab'))
            else:
                pass
            if len(samplelist) < len(self.option('group').prop['sample']):
                raise OptionError('分组文件中样本数量:%s多于距离矩阵中的样本数量:%s' % (len(self.option('group').prop['sample']),
                                  len(samplelist)))
            for sample in self.option('group').prop['sample']:
                if sample not in samplelist:
                    raise OptionError('分组文件的样本(%s)在距离矩阵的样本中不存在' % sample)
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "距离统计和统计检验分析结果目录"],
            ["./Stats.xls", "xls", "分组统计检验结果"],
            ["./Distances.xls", "xls", "组内组间距离值统计结果"]
        ])
        # print self.get_upload_files()
        super(DistanceBoxAgent, self).end()


class DistanceBoxTool(Tool):
    def __init__(self, config):
        super(DistanceBoxTool, self).__init__(config)
        self._version = '1.9.1'  # qiime版本
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.cmd_path = 'miniconda2/bin/make_distance_boxplots.py'
        self.grouplab = self.option('grouplab') if self.option('grouplab') else self.option('group').prop['group_scheme'][0]
        self.dis_matrix = self.get_matrix()

    def get_matrix(self):
        if len(self.option('dis_matrix').prop['samp_list']) == len(self.option('group').prop['sample']):
            return self.option('dis_matrix').path
        else:
            self.option('dis_matrix').create_new(self.option('group').prop['sample'],
                                                 os.path.join(self.work_dir, 'dis_matrix_filter.temp'))
            return os.path.join(self.work_dir, 'dis_matrix_filter.temp')

    def run(self):
        """
        运行
        """
        super(DistanceBoxTool, self).run()
        self.run_box()

    def linkfile(self, oldfile, newname):
        """
        link文件到output文件夹
        :param oldfile: 资源文件路径
        :param newname: 新的文件名
        :return:
        """
        newpath = os.path.join(self.output_dir, newname)
        if os.path.exists(newpath):
            os.remove(newpath)
        os.link(oldfile, newpath)

    def run_box(self):
        """
        运行qiime:make_distance_boxplots.py
        """
        cmd = self.cmd_path
        cmd += ' -m %s -d %s -o %s -f %s --save_raw_data -n %s' % (self.option('group').path,
                                                                   self.dis_matrix,
                                                                   self.work_dir,
                                                                   self.grouplab,
                                                                   self.option('permutations'))
        self.logger.info('运行qiime:make_distance_boxplots.py程序')
        box_command = self.add_command('box', cmd)
        box_command.run()
        self.wait(box_command)
        if box_command.return_code == 0:
            self.logger.info('运行qiime/make_distance_boxplots.py完成')
            self.format_boxdata(self.work_dir + '/' + self.grouplab + '_Distances.txt',
                                self.work_dir + '/' + self.grouplab + '_Distances_format.txt')
            self.linkfile(self.work_dir + '/' + self.grouplab + '_Distances_format.txt',
                          'Distances.xls')
            self.format_stats(self.work_dir + '/' + self.grouplab + '_Stats.txt',
                              self.work_dir + '/' + self.grouplab + '_Stats_format.txt')
            self.linkfile(self.work_dir + '/' + self.grouplab + '_Stats_format.txt', 'Stats.xls')
            self.end()
        else:
            self.set_error('运行qiime/make_distance_boxplots.py出错')

    def format_boxdata(self, boxdata, outfile):
        """
        """
        datalist = []
        with open(boxdata, 'rb') as f, open(outfile, 'wb') as w:
            for line in f:
                templist = line.rstrip().split('\t')
                values = []
                for i in templist[1:]:
                    try:
                        values.append(float(i))
                    except ValueError:
                        self.set_error('qiime计算箱线数据结果值存在异常，数据：%s应该为数字' % i)
                datalist.append([templist[0]] + self.calculate_boxdata(values))
            w.write('#Group\tmax\tmin\tq3\tq1\tmedian\tfliers\tmean\tstd\n')
            for box in datalist:
                # print box, type(box)
                w.write('\t'.join(box) + '\n')

    def format_stats(self, statfile, outfile):
        """
        去除Stats.txt文件中的名称的空格为下划线“_”
        """
        try:
            with open(statfile, 'rb') as f, open(outfile, 'wb') as w:
                for line in f:
                    if line[0] == '#':
                        continue
                    elif re.match('Group 1\tGroup 2\t', line):
                        w.write(line)
                    else:
                        split_line = line.split('\t')
                        split_line[0] = split_line[0].replace(' ', '_')
                        split_line[1] = split_line[1].replace(' ', '_')
                        w.write('\t'.join(split_line))
        except IOError:
            raise Exception('格式化Stats文件无法打开，或者写入文件无法成功')

    def calculate_boxdata(self, datalist):
        """
        """
        import pandas as pd
        import numpy as np
        data = pd.DataFrame({'name': datalist})
        result = data.boxplot(return_type='dict')
        medians = result['medians'][0].get_data()[1][0]
        max_v = result['whiskers'][1].get_data()[1][1]
        min_v = result['whiskers'][0].get_data()[1][1]
        q3_v = result['whiskers'][1].get_data()[1][0]
        q1_v = result['whiskers'][0].get_data()[1][0]
        fliers = ','.join([str(i) for i in list(result['fliers'][0].get_data()[1])])
        means_v = np.mean(datalist, axis=0)
        std_v = np.std(datalist, axis=0)
        return [str(max_v), str(min_v), str(q3_v), str(q1_v), str(medians), fliers, str(means_v), str(std_v)]
