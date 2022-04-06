# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20190307

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import os
import re
import random
from collections import defaultdict


class CnvCompareAgent(Agent):
    """
    用于交互分析中cnv变异检测--cnv比较分析
    新版本
    """
    def __init__(self, parent):
        super(CnvCompareAgent, self).__init__(parent)
        options = [
            {"name": "infile_path", "type": "string"},
            {"name": "region", "type": "string"},  # chr1:1-500,chr2:2-4
            {"name": "marktype", "type": "string"},  # same or diff or all same,diff
            {"name": "samples", "type": "string"},  # 样本对a|b,b|c
            {"name": "region_type", "type": "string", 'default': 'real_region'},  # 用于后面区域过滤是过滤范
            # 围还是具体的点
            {"name": "analysis_model", "type": "string", 'default': "multiple"}
        ]
        self.add_option(options)
        self.step.add_steps('cnvdiff')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cnvdiff.start()
        self.step.update()

    def step_end(self):
        self.step.cnvdiff.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("samples"):
            raise OptionError("缺少samples参数")
        if not self.option("infile_path"):
            raise OptionError("缺少infile_path参数")
        if not self.option("region"):
            raise OptionError("缺少region参数")
        if not self.option("marktype"):
            raise OptionError("缺少marktype参数")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '30G'
        
    def end(self):
        super(CnvCompareAgent, self).end()


class CnvCompareTool(Tool):
    def __init__(self, config):
        super(CnvCompareTool, self).__init__(config)
        self.cnv_diff_path = self.config.PACKAGE_DIR + "/wgs_v2/cnv_compare.py"
        self.python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/python '
        
    def cnv_diff(self):
        """
        python cnv_compare.py -i ./ -r all -m 'same,diff' -s 'AH03|AH19,CZ02|AH03' -o ./
        :return:
        """
        cmd1 = "{}{} -i {} -r {} -m '{}' -s '{}' -o {} -rt {}"\
            .format(self.python_path, self.cnv_diff_path, self.option("infile_path"),
                    self.option("region"), self.option("marktype"),
                    self.option("samples"), self.output_dir, self.option('region_type'))
        self.logger.info("cmd:{}".format(cmd1))
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "cnv_diff_{}.sh".format(now_time)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd1)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path) )
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始进行cnv_diff")
        command1 = self.add_command("cnv_diff", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("cnv_diff完成！")
        else:
            self.set_error("cnv_diff出错！")
        os.system('rm {}'.format(file_path))

    def multiple_analysis(self):
        """
        如果是多条件比较的时候，要将各自比较的结果进行求交集
        :return:
        """
        cnv = defaultdict(list)
        gene = defaultdict(list)
        sample_title = []
        samples = self.get_samples_num()
        analysis_times = len(self.option("marktype").split(','))
        marktype = self.option("marktype").split(',')[0]
        sample_1, sample_2 = self.option("samples").split(',')[0].split('|')
        first_file = os.path.join(self.output_dir, '_VS_'.join([sample_1, sample_2]) + '_{}.xls'.format(marktype))
        out_file = os.path.join(self.output_dir, "multiple_same.xls")
        if len(open(first_file, 'r').readlines()) < 2:
            self.set_error("比较分析结果为空，不进行后面的分析，请设置筛选条件重新运行！")
        self.logger.info("第一个文件的路径：{}".format(first_file))
        with open(first_file, 'r') as r:
            for line in r:
                temp = line.strip().split('\t')
                if line.startswith('#'):
                    sample_title.extend([temp[5], temp[6], temp[7], temp[8]])
                else:
                    cnv['\t'.join([temp[0], temp[1], temp[2], temp[3], temp[4]])] = \
                        [temp[5], temp[6], temp[7], temp[8]]
                    gene['\t'.join([temp[0], temp[1], temp[2], temp[3], temp[4]])] = [temp[9], temp[10]]
        for m in os.listdir(self.output_dir):
            if re.match(".*same\.xls|.*diff\.xls", m):
                if m in ['_'.join([sample_1, sample_2, marktype]) + '.xls', "multiple_same.xls"]:
                    continue
                else:
                    with open(os.path.join(self.output_dir, m), 'r') as r:
                        for line in r:
                            temp = line.strip().split('\t')
                            if line.startswith('#'):
                                sample_title.extend([temp[5], temp[6], temp[7], temp[8]])
                            else:
                                cnv_ = '\t'.join([temp[0], temp[1], temp[2], temp[3], temp[4]])
                                if cnv_ in cnv.keys():
                                    cnv[cnv_].extend([temp[5], temp[6], temp[7], temp[8]])
        # self.logger.info(sample_title)
        # self.logger.info(cnv)
        with open(out_file, 'w') as w:
            new_title = []
            for j in samples:
                new_title.extend([j + "_genotype", j + "_pvalue"])
            w.write('#chr\tstart\tend\tlength\ttype\t{}\tgene_num\tgene\n'.format('\t'.join(new_title)))
            for key in sorted(cnv.keys()):
                if len(cnv[key]) == analysis_times * 4:
                    ge_pv = []
                    for j in samples:
                        ge = sample_title.index(j + "_genotype")
                        pv = sample_title.index(j + "_pvalue")
                        ge_pv.extend([cnv[key][ge], cnv[key][pv]])
                    w.write("{}\t{}\t{}\n".format(key, '\t'.join(ge_pv), "\t".join(gene[key])))
        self.cnv_stat(self.output_dir + "/multiple_same.xls")
        for m in os.listdir(self.output_dir):
            if m not in ["multiple_same.stat.xls", "multiple_same.xls"]:
                os.remove(os.path.join(self.output_dir, m))

    def cnv_stat(self, file_path):
        """
        对cnv结果进行差异统计
        :return:
        """
        chr_info = defaultdict(list)
        out = self.output_dir + '/multiple_same.stat.xls'
        with open(file_path, 'r') as r, open(out, 'w') as w:
            w.write("#chrid\tdeletion\tduplication\tgene\n")
            for line in r:
                if not re.match('^#.*', line):
                    temp = line.strip().split('\t')
                    if temp[0] in chr_info.keys():
                        if temp[4] == 'deletion':
                            chr_info[temp[0]][0] += 1
                        if temp[4] == 'duplication':
                            chr_info[temp[0]][1] += 1
                        # noinspection PyBroadException
                        try:
                            chr_info[temp[0]][2] += int(temp[9])
                        except:
                            chr_info[temp[0]][2] += 0
                    else:
                        chr_info[temp[0]] = [0, 0, 0]
            for key in sorted(chr_info.keys()):
                w.write('{}\t{}\t{}\t{}\n'.format(key, chr_info[key][0], chr_info[key][1], chr_info[key][2]))

    def get_samples_num(self):
        samples = []
        for m in self.option("samples").split(','):
            for n in m.split("|"):
                if n not in samples:
                    samples.append(n)
        return samples

    def run(self):
        super(CnvCompareTool, self).run()
        self.cnv_diff()
        if self.option("analysis_model") == 'multiple' and len(self.option("marktype").split(',')) != 1:
            self.multiple_analysis()
        self.end()
