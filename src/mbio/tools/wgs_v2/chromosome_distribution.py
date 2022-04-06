# -*- coding: utf-8 -*-
# __author__: HONGDONG
# modified: 20190318

import re
import os
import time
import functools
from biocluster.agent import Agent
from biocluster.tool import Tool
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


def time_count(func):
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run ' + func_name + ' at ' + start_time)
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End ' + func_name + ' at ' + end_time)
        print("{}函数执行时间约为{}s".format(func.__name__, end - start))
    return wrapper


class ChromosomeDistributionAgent(Agent):
    """
    染色体标记分布图
    """
    def __init__(self, parent):
        super(ChromosomeDistributionAgent, self).__init__(parent)
        options = [
            {"name": "pos_file", "type": "infile", 'format': "bsa.vcf", "required": True},  # 各种情况的位置信息文件
            {"name": "marker_type", "type": "string", "required": True},
            {"name": "data_type", "type": "string", "required": True},
            {"name": "win_step", "type": "int", 'default': 1000000, "required": True},
            {"name": "analysis_object", "type": "string", "required": True},
            {"name": "graphic_style", "type": "string", "required": True},
            {"name": "target_path", 'type': "string"},
            {"name": "update_info", "type": "string", "required": True},
            {"name": "main_id", "type": "string", "required": True},
            {"name": "project_type", 'type': "string"},
            {"name": "task_id", "type": "string"},
            {"name": "project_sn", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("marker_type") not in ['Gene', 'SNP+InDel', "SNP", "InDel", "CNV", "SV", "SSR"]:
            raise OptionError("marker_type：%s类型不合法" % self.option("marker_type"))
        if self.option("graphic_style") not in ['heatmap', 'area']:
            raise OptionError("图形样式参数%s不正确" % self.option("graphic_style"))
        if self.option("data_type") not in ['before', 'after', 'ref']:
            raise OptionError("分析数据来源参数%s不正确" % self.option("data_type"))
        if self.option("marker_type") != 'SNP+InDel' and self.option("analysis_object") == 'all':
            raise OptionError("标记类型不为'SNP+InDel'时选择对象参数不能为all")
        if self.option("data_type") == 'ref' and self.option("marker_type") not in ["Gene", "SSR"]:
            raise OptionError("当数据类型是ref的时候，标记类型必须为Gene 或者SSR！")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(ChromosomeDistributionAgent, self).end()


class ChromosomeDistributionTool(Tool):
    def __init__(self, config):
        super(ChromosomeDistributionTool, self).__init__(config)
        self.window_sliding = self.config.PACKAGE_DIR + "/wgs_v2/window_sliding.py"
        self.python_path = 'miniconda2/bin/python '
        self.pos_file = ""

    def run_window_sliding(self):
        cmd = '{} {} -i {} -s {} -o {}'.format(self.python_path, self.window_sliding, self.pos_file,
                                               self.option("win_step"), self.output_dir)
        self.logger.info("开始进行window_sliding")
        command = self.add_command("window_sliding", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("window_sliding完成！")
        else:
            self.set_error("window_sliding出错！")

    @time_count
    def set_db(self):
        self.logger.info("开始进行导表")
        api = self.api.api("wgs_v2.chromosome_distribution")
        if self.option("project_type") and self.option("project_type") != 'dna_wgs_v2':
            api._project_type = self.option("project_type")
        for m in os.listdir(self.output_dir):
            if self.option("graphic_style") == 'area':
                api.add_sg_area_detail(self.option("main_id"), self.output_dir + "/{}".format(m),
                                       self.option("project_sn"), self.option("task_id"),
                                       self.work_dir + "/chr_start.txt")
            else:
                api.add_sg_distribution_detail(self.option("main_id"), self.output_dir + "/{}".format(m),
                                               self.option("project_sn"), self.option("task_id"),
                                               self.work_dir + "/chr_start.txt")
        self.logger.info("导表完成")

    def run(self):
        super(ChromosomeDistributionTool, self).run()
        if self.option("data_type") == 'before':
            self.make_file_before()
        elif self.option("data_type") == 'after':
            self.make_file_after()
        else:
            self.make_file_ref()
        self.run_window_sliding()
        if self.option("main_id"):
            self.set_db()
        self.end()

    def make_file_before(self):
        """
        设置比较分析前原始数据的位置信息
        :return:
        """
        if self.option("marker_type") in ['SNP+InDel', 'SNP', "InDel"]:
            self.one_sample_snp(self.option("pos_file").prop['path'], self.option("marker_type"),
                                self.option("analysis_object"))
        elif self.option("marker_type") == 'CNV':
            self.pos_file = self.option("pos_file").prop['path']
        elif self.option("marker_type") == 'SV':
            self.one_sample_sv(self.option("pos_file").prop['path'], self.option("analysis_object"))
        elif self.option("marker_type") == 'SSR':
            self.one_sample_ssr(self.option("pos_file").prop['path'], self.option("analysis_object"))
        pass

    def make_file_after(self):
        """
        设置比较分析后的位置信息
        :return:
        """
        if self.option("marker_type") in ['SNP+InDel', 'SNP', "InDel"]:  # variant_compare.detail
            self.pos_file = self.option("pos_file").prop['path']
        elif self.option("marker_type") == 'CNV':
            self.pos_file = self.option("pos_file").prop['path']
        elif self.option("marker_type") == 'SV':
            self.sv_after(self.option("pos_file").prop['path'])
        elif self.option("marker_type") == 'SSR':
            self.pos_file = self.option("pos_file").prop['path']

    def make_file_ref(self):
        """
        设置参考基因组的位置信息
        :return:
        """
        self.ref(self.option("pos_file").prop['path'], self.option("marker_type"))

    @time_count
    def sv_after(self, file_path):
        """
        Sv比较分析的结果文件/mnt/ilustre/users/sanger-dev/workspace/20190227/Single_sv_compare/S
        vCompare/output/SRR5739119_VS_SRR5739120.detail.xls
        :param file_path:
        :return:
        """
        self.logger.info("开始进行SV比较分析提取信息")
        self.pos_file = self.work_dir + '/chr_start.txt'
        with open(file_path, 'r') as r, open(self.pos_file, 'w') as w:
            for line in r:
                if not re.match("#.*", line):
                    temp = line.strip().split("\t")
                    w.write('{}\t{}\n'.format(temp[1], temp[2]))
        self.logger.info("开始进行SV比较分析提取信息成功")

    @time_count
    def one_sample_snp(self, file_path, type_, sample_name):
        """
        获取原始的pop.final.vcf中的某一个样本的chr_start信息
        :param file_path:
        :param type_:
        :param sample_name:
        :return:
        """
        self.logger.info("开始进行原始vcf提取信息")
        self.pos_file = self.work_dir + '/chr_start.txt'
        index_ = 0
        with open(file_path, 'r') as r, open(self.pos_file, 'w') as w:
            for line in r:
                if sample_name != 'all':
                    if re.match('#CHROM.*', line):
                        sample_list = line.strip().split("\t")
                        if sample_name in sample_list:
                            index_ = sample_list.index(sample_name)
                            self.logger.info('index:{}'.format(index_))
                        else:
                            raise Exception("样本不在 {} vcf {} 文件中!".format(sample_name, file_path))
                if not re.match("#.*", line):
                    temp = line.strip().split("\t")
                    if sample_name != 'all':
                        if type_ != 'SNP+InDel':
                            marker_type = self.check_snp_indel(','.join([temp[3], temp[4]]))
                            if marker_type != type_:
                                continue
                        if temp[6] not in ['PASS']:
                            continue
                        if temp[index_].split(":")[0] not in ['0/0']:
                            w.write('{}\t{}\n'.format(temp[0], temp[1]))
                    else:
                        w.write('{}\t{}\n'.format(temp[0], temp[1]))
        self.logger.info("开始进行原始vcf提取信息完成！")

    @time_count
    def one_sample_sv(self, file_path, sample_name):
        """
        获取sv的vcf文件中指定样本的chr_start   pop.sort.sv.vcf
        :param file_path:
        :param sample_name:
        :return:
        """
        self.logger.info("开始进行原始SV提取信息")
        self.pos_file = self.work_dir + '/chr_start.txt'
        with open(file_path, 'r') as r, open(self.pos_file, 'w') as w:
            for line in r:
                if re.match('#CHROM.*', line):
                    sample_list = line.strip().split("\t")
                    if sample_name in sample_list:
                        index_ = sample_list.index(sample_name)
                        self.logger.info('index:{}'.format(index_))
                    else:
                        raise Exception("样本不在 {} vcf {} 文件中!".format(sample_name, file_path))
                if not re.match("#.*", line):
                    temp = line.strip().split("\t")
                    if temp[index_].split(":")[0] not in ['0/0'] or temp[index_].split(":")[3] in ['PASS']:
                        w.write('{}\t{}\n'.format(temp[0], temp[1]))
        self.logger.info("进行原始SV提取信息成功")

    @time_count
    def one_sample_ssr(self, file_path, sample_name):
        """
        获取ssr的vcf文件中指定样本的chr_start   pop.sort.sv.vcf
        :param file_path:
        :param sample_name:
        :return:
        """
        self.logger.info("开始进行原始SSR提取信息")
        self.pos_file = self.work_dir + '/chr_start.txt'
        with open(file_path, 'r') as r, open(self.pos_file, 'w') as w:
            index_ = 0
            for line in r:
                if re.match('#CHROM.*', line):
                    sample_list = line.strip().split("\t")
                    if sample_name in sample_list:
                        index_ = sample_list.index(sample_name)
                        self.logger.info('index:{}'.format(index_))
                    else:
                        raise Exception("样本不在 {} vcf {} 文件中!".format(sample_name, file_path))
                if not re.match("#.*", line):
                    temp = line.strip().split("\t")
                    if temp[index_] not in ['-:-:-:-:-']:
                        w.write('{}\t{}\n'.format(temp[0], temp[1]))
        self.logger.info("进行原始SSR提取信息成功")

    @time_count
    def ref(self, file_path, type_):
        """
        获取参考基因组的位置信息，Gene的时候，ref.gtf， SSR的时候ssr.ref.result.xls
        :param file_path:
        :param type_:
        :return:
        """
        self.logger.info("开始进行ref提取信息")
        self.pos_file = self.work_dir + '/chr_start.txt'
        with open(file_path, 'r') as r, open(self.pos_file, 'w') as w:
            for line in r:
                if not re.match("#.*", line):
                    temp = line.strip().split("\t")
                    if type_ == 'Gene' and temp[2] in ['gene', "CDS"]:
                        w.write('{}\t{}\n'.format(temp[0], temp[3]))
                    elif type_ == 'SSR':
                        w.write('{}\t{}\n'.format(temp[0], temp[5]))
                    else:
                        pass
        self.logger.info("进行ref提取信息成功")

    def check_snp_indel(self, data):
        """
        判断该点是snp还是indel
        :param data:
        :return:
        """
        marker_type = "SNP"
        for n in data.split(','):
            if len(n) > 1:
                marker_type = "InDel"
        return marker_type
