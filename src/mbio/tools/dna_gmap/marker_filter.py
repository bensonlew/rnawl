# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modify 20180623
# modify 20180707
# tool

from bson.objectid import ObjectId
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os


class MarkerFilterAgent(Agent):
    """
    先运行tool:vcf_marker.py，传给vcf PID MID父母本ID；
    work_flow内的pop.primary.marker文件为不过滤文件。过滤都在这个接口的tool里进行
    若：
        接口某个参数缺失，会携带之前那个主表内的params传进来。pl和tool内的默认值和起作用
    """
    def __init__(self, parent):
        super(MarkerFilterAgent, self).__init__(parent)
        options = [
            {"name": "vcf", "type": "infile", "format": "dna_gmap.vcf"},  # 传入vcf或vcf.gz
            {"name": "detail_info", "type": "infile", "format": "dna_gmap.marker"},     # 主表里detail_info_path
            {"name": "type", "type": "string"},     # SNP/INDEL/ALL in controller
            {"name": "pdep", "type": "string"},     # 亲本深度1_20 or _20 or 1
            {"name": "odep", "type": "string"},     # 子代深度1_20 or _20 or 1
            {"name": "popt", "type": "string"},     # 群体类型
            {"name": "miss_tatio", "type": "float"},   # 缺失率：controller将[0,100]/100传过来
            {"name": "signif", "type": "float"},  # P值：[0,1]
            {"name": "marker_upload", "type": "infile", "format": "dna_gmap.marker"},  # 客户上传marker文件;
            # {"name": "marker_upload", "type": "string"},
            {"name": "child_list", "type": "string"},      # 子代列表的objectID str
            {"name": "filter_marker", "type": "outfile", "format": "dna_gmap.marker"}      # outfile
        ]
        self.add_option(options)
        self.step.add_steps('MarkerFilter')
        self._cpu = ''
        self._memory = ''
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.MarkerFilter.start()
        self.step.update()

    def step_end(self):
        self.step.MarkerFilter.finish()
        self.step.update()

    def check_options(self):
        if not self.option("vcf"):
            raise OptionError("请设置vcf文件", code="34800701")
        if not self.option("detail_info"):
            raise OptionError("请设置detail_info文件", code="34800702")
        if self.option("type") not in ["SNP", "INDEL", "ALL"]:
            raise OptionError("群体类型%s不属于SNP, INDEL, ALL", variables=(self.option("type")), code="34800703")
        if not self.option("pdep"):
            raise OptionError("请设置亲本平均测序深度,如10_20 or _20 or 10", code="34800704")
        if not self.option("odep"):
            raise OptionError("请设置子代平均测序深度,如1 or 1_30 or _30", code="34800705")
        if re.match(r'ri\d+', self.option("popt"), re.I):
            pass
        else:
            if self.option("popt").upper() not in ["BC", "DH", "CP", "F2", "F1"]:
                raise OptionError("群体类型%s不属于BC, DH, CP, F2, RIL, Ri\d+", variables=(self.option("popt")), code="34800706")
        if self.option("miss_tatio") != '':     # 参数设置页面选填，controller里判断
            if self.option("miss_tatio") < 0 or self.option("miss_tatio") > 100:
                raise OptionError("%s设置错误，需在[0,100]范围内", variables=(self.option("miss_tatio")), code="34800707")
        if self.option("signif") != "":     # 同上
            if self.option("signif") < 0 or self.option("signif") > 1:
                raise OptionError("%s设置错误，需在[0,1]范围内", variables=(self.option("signif")), code="34800708")

    def set_resource(self):
        """
        申请的运行所需资源。
        """
        self._cpu = 2
        self._memory = '40G'

    def end(self):
        super(MarkerFilterAgent, self).end()


class MarkerFilterTool(Tool):
    def __init__(self, config):
        super(MarkerFilterTool, self).__init__(config)
        self.perl_path = 'miniconda2/bin/perl '
        self.mkerfilter_path = self.config.PACKAGE_DIR + "/dna_gmap/markerfilter.pl"
        self.markerinfo_path = self.config.PACKAGE_DIR + "/dna_gmap/markerinfo_tian.pl"
        self.sort_markinfo_path = self.config.PACKAGE_DIR + "/dna_gmap/sort_markinfo.pl"

    # def get_child_list(self):
    #     """
    #     return子代列表的样本id:self.sample_child
    #     """
    #     marker_filter_api = self.api.api("dna_gmap.marker_filter")
    #     if self.option('child_list'):
    #         self.logger.info("^^^^^ child_list:{}存在".format(self.option('child_list')))
    #         child_list_api = marker_filter_api.get_child_list(ObjectId(self.option('child_list')))
    #         child_list = child_list_api["spcimen_ids"]
    #         child_list = list(set(child_list))
    #         child_list.sort()
    #         self.child_list = ','.join(child_list)

    def check_sample_upload(self):
        """
        check传入文件的sample id
        """
        if self.option('marker_upload').is_set:
            print(self.option('marker_upload').prop['path'])
            self.sample_vcf = []
            self.sample_upload = []     # 无marker_upload文件，该list为空
            sample_notvcf = []
            with open(self.option('vcf').prop['path'], "r") as f:
                self.logger.info("vcf读入ok")
                line = f.readline()
                while line:
                    if re.match('#CHROM', line):
                        item = line.strip().split("\t")
                        self.sample_vcf = item[9:]
                        break
                    line = f.readline()
            self.logger.info("^^^^^2 vcf样品为{}".format(self.sample_vcf))
            if len(self.sample_vcf) <= 2:
                self.set_error("遗传图谱的项目分析样本数必须多于三个", code="34800701")
            with open(self.option('marker_upload').prop['path'], "r") as f:
                lines = f.readlines()
                header = lines[0].strip().split("\t")
                self.sample_upload = header[2:]      # 用子代列表筛完上传的样品id,剩下的是否为空
                if len(self.sample_upload) != 0:
                    for i in self.sample_upload:
                        if i not in self.sample_vcf:
                            sample_notvcf.append(i)
                        self.logger.info("^^^^^ 上传文件内id：{}不在vcf,从而被剔除".format(sample_notvcf))
                else:
                    self.set_error("样本为空%s", variables=(self.sample_upload), code="34800702")

    def get_popt(self, popt):
        """
        f1输入F1
        :return 非F1代结果
        """
        if re.match(r'ri\d+', popt, re.I):
            # return 'Ri10'
            return popt.upper()
        else:
            if popt.upper() in ["BC", "DH", "CP", "F2", 'F1']:
                if popt.upper() == 'F1':
                    popt = 'CP'
                return popt.upper()

    def run_markerfilter(self):
        """
        统计页面1.1 1.2 1.3 1.5数据
        markerfilter.pl
        self.option('miss_tatio')传入的是 3e-05 要转成 '0.00003'
        """
        detail_info_path = self.option("detail_info").prop['path']
        popt = self.get_popt(self.option('popt'))
        cmd = "{} {} -info {}".format(self.perl_path, self.mkerfilter_path, detail_info_path)
        cmd += " -out {} -pop {}".format(self.output_dir + "/pop.filtered", popt)
        cmd += " -seg {} -mis {}".format(self.option('signif'), '{:.5f}'.format(self.option('miss_tatio')))
        cmd += " -vtype {}".format(self.option('type'))
        for ii in ['pdep', 'odep']:
            range_ = self.option(ii).strip().split('_')
            min_ = range_[0] if range_[0] != '' else '0'
            if len(range_) == 2:
                max_ = range_[1] if range_[1] != '' else 1000000
            else:
                max_ = 1000000
            min__ = 0.001 if min_ == '0' else int(min_)
            if ii == 'pdep':
                cmd += ' -minPMdep {} -maxPMdep {}'.format(min__, int(max_))
            else:
                cmd += ' -minOdep {} -maxOdep {}'.format(min__, int(max_))
        if self.option('child_list'):
            cmd += ' -indilist {}'.format(self.option('child_list'))
        if self.option('marker_upload').is_set:
            cmd += ' -newmarker {}'.format(self.option('marker_upload').prop['path'])
        self.run_cmd(cmd, "markerfilter")

    def run_markerinfo(self):
        """
        1.4
        cmd：markerinfo.pl
        :input vcf pop.filtered.marker
        :return  marker_info.xls
        生成标记数据统计表
        """
        self.check_exists(self.output_dir + "/pop.filtered.marker")
        # self.check_file_lines(self.output_dir + "/pop.filtered.marker")     # 核查文件是否只有抬头
        self.option("filter_marker", self.output_dir + "/pop.filtered.marker")
        cmd = "{} {} -vcf {} ".format(self.perl_path, self.markerinfo_path, self.option('vcf').prop['path'])
        cmd += "-mark {} -out {}".format(self.output_dir + "/pop.filtered.marker", self.work_dir + "/marker_info.xls")
        self.run_cmd(cmd, "markerinfo")

    def run_markerinfo_sort(self):
        """
        将流程生成的markinfo文件排序
        """
        self.check_exists(self.work_dir + "/marker_info.xls")
        cmd = "{} {} -i {}".format(self.perl_path, self.sort_markinfo_path, self.work_dir + "/marker_info.xls")
        cmd += " -o {}".format(self.output_dir + "/marker_info.xls")
        self.run_cmd(cmd, "markerinfo_sort")

    def check_exists(self, file_path):
        """
        用于检查文件及文件夹是否存在
        """
        if not os.path.exists(file_path):
            self.set_error("文件或文件夹%s不存在！", variables=(file_path), code="34800703")
            self.set_error("文件或文件夹%s不存在！", variables=(file_path), code="34800709")

    # def check_file_lines(self, mark_path):
    #     with open(mark_path, "r") as f:
    #         lines = f.readlines()
    #         num = len(lines)
    #         if num <= 1:
    #             self.set_error("文件:{}不得少于两行，请检查".format(mark_path))
    #     return True

    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("%s运行失败", variables=(cmd_name), code="34800704")
            self.set_error("%s运行失败", variables=(cmd_name), code="34800710")

    def run(self):
        """
        1 得到子代列表,分割
        2 检查上传样品的id,扔掉不在vcf的样品，筛完后样本为空报错
        3 程序开始过滤
        1.4 统计标记信息
        """
        super(MarkerFilterTool, self).run()
        # self.get_child_list()       # 1
        self.check_sample_upload()  # 2
        self.run_markerfilter()     # 3
        self.run_markerinfo()       # 1.4
        self.run_markerinfo_sort()  # 1.4排序
        self.end()
