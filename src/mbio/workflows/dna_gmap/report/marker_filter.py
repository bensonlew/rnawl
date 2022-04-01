# -*- coding: utf-8 -*-
# __author__ = 'qing_mei'
# modified 20180709
# workflow

import re
import os
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class MarkerFilterWorkflow(Workflow):
    """
    标记筛选的接口
    sg_marker_id
    F1流程输入内需要输入CP
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MarkerFilterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "vcf", "type": "infile", "format": "dna_gmap.vcf"},  # 传入vcf或vcf.gz
            {"name": "detail_info", "type": "infile", "format": "dna_gmap.marker"},     # 主表里detail_info_path
            {"name": "type", "type": "string"},     # sNP/INDEL/ALL in controller
            {"name": "pdep", "type": "string"},     # 亲本深度1_20 or _20 or 1
            {"name": "odep", "type": "string"},     # 子代深度1_20 or _20 or 1
            {"name": "popt", "type": "string"},     # 群体类型
            {"name": "miss_tatio", "type": "float"},   # 缺失率：controller将[0,100]/100传过来
            {"name": "signif", "type": "float"},  # P值：[0,1]
            {"name": "marker_upload", "type": "infile", "format": "dna_gmap.marker"},  # 客户上传marker文件;
            # {"name": "marker_upload", "type": "string"},
            {"name": "child_list", "type": "string"},      # 子代列表的objectID
            {"name": "update_info", "type": "string"},
            {"name": "ref_chrlist_path", "type": "string"},
            {"name": "filter_marker", "type": "outfile", "format": "dna_gmap.marker"},      # outfile
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.marker_filter = self.add_tool('dna_gmap.marker_filter')
        self.target_dir_ = ""

    def check_options(self):
        if not self.option("main_id"):
            raise OptionError("请设置main_id", code="14800401")
        if not self.option("task_id"):
            raise OptionError("请设置task_id", code="14800402")
        if not self.option("vcf"):
            raise OptionError("请设置vcf文件", code="14800403")
        if not self.option("detail_info").is_set:   # infile类型判断时必须写生is_set
            raise OptionError("请设置detail_info文件", code="14800404")
        if self.option("type") not in ["SNP", "INDEL", "ALL"]:
            raise OptionError("群体类型%s不属于SNP, INDEL, ALL", variables=(self.option("type")), code="14800405")
        if not self.option("pdep"):
            raise OptionError("请设置亲本平均测序深度,如10_20 or _20 or 10", code="14800406")
        if not self.option("odep"):
            raise OptionError("请设置子代平均测序深度,如1 or 1_30 or _30", code="14800407")
        if re.match(r'ri\d+', self.option("popt"), re.I):
            pass
        else:
            if self.option("popt").upper() not in ["BC", "DH", "CP", "F2", "F1"]:
                raise OptionError("群体类型%s不属于BC, DH, CP, F2, F1, RIL, Ri\d+", variables=(self.option("popt")), code="14800408")
        if self.option("miss_tatio") != '':     # 因参数设置页面选填，controller里判断是否为None
            if self.option("miss_tatio") < 0 or self.option("miss_tatio") > 1:
                raise OptionError("%s设置错误，需在[0,100]范围内", variables=(self.option("miss_tatio")), code="14800409")
        if self.option("signif") != "":     # 同上
            if self.option("signif") < 0 or self.option("signif") > 1:
                raise OptionError("%s设置错误，需在[0,1]范围内", variables=(self.option("signif")), code="14800410")
        if not self.option("ref_chrlist_path"):
            raise OptionError("请核查sg_task内的ref_chrlist_path", code="14800411")
        return True

    def get_file_path(self):
        """
        集中处理
        detail_info表结构内的path,rere
        vcf_path也要变更路径后然后传给tool的option
        """
        self.detail_info = self.option("detail_info").prop['path']  # infile后需要加prop['path']
        self.vcf_path = self.option('vcf').prop['path']  # infile后需要加prop['path']
        # self.chr_list = []  # 无用
        # with open(self.option("ref_chrlist_path"), 'r') as f:
        #     lines = f.readlines()
        #     for i in lines:
        #         chr = i.strip().split('\t')
        #         if chr[0] not in self.chr_list:
        #             self.chr_list.append(chr[0])
        self.filtered_marker_path = self.output_dir + "/pop.filtered.marker"

    def marker_filter_run(self):
        """
        缺失率和P值为什么是必填选项？
        首先如果参数页面不输入的话，代表params里2参数为"，
        但是markerfilter.pl会默认两个参数为空时分别赋值0.3 0.05。
        因此和表内参数为""矛盾，故在controller那里就判断参数为"时赋值0.3 0.05，
        并存进去parmas内，供页面展示。
        然后再传给workflow,因此为必须参数。
        """
        opt = {
            "vcf": self.vcf_path,
            "detail_info": self.detail_info,
            "type": self.option('type'),
            "pdep": self.option('pdep'),
            "odep": self.option('odep'),
            "popt": self.option('popt').upper(),
            "miss_tatio": self.option('miss_tatio'),
            "signif": self.option('signif'),
        }
        if self.option('marker_upload').is_set:
            opt.update({
                "marker_upload": self.option('marker_upload').prop['path'],
            })
        if self.option('child_list'):
            opt.update({
                "child_list": self.option('child_list'),    # str
            })
        self.marker_filter.set_options(opt)
        self.marker_filter.on('end', self.set_output, 'marker_filter')
        self.marker_filter.run()

    # def get_chrlist(self, chrlist, filtered_marker_path):
    #     """
    #     将chrlist的id排在前面，这个chrid前提存在于marker文件内
    #     """
    #     marker_chrlist = []
    #     with open(filtered_marker_path, 'r') as f:
    #         lines = f.readlines()
    #         for line in lines[1:]:
    #             line_list = line.strip().split('_')
    #             if line_list[0] not in marker_chrlist:
    #                 marker_chrlist.append(line_list[0])
    #     return marker_chrlist

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'marker_filter':
            self.linkdir(obj.output_dir, self.output_dir)
        self.option("filter_marker", self.output_dir + "/pop.filtered.marker")
        # self.chr_list_ = self.get_chrlist(self.chr_list, self.output_dir + "/pop.filtered.marker")
        self.set_db()
        # self.end()  # 导表后在导表内end
        pass

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        """
        更新主表的chr_list和path:marker.filtered.marker
        self.filtered_marker_path
        """
        self.get_target_dir()
        marker_id = ObjectId(self.option("main_id"))
        self.logger.info("开始marker_filter的导表！")
        api = self.api.api("dna_gmap.marker_filter")
        task_id = self.option('task_id')
        api.get_update_chrlist(marker_id, self.target_dir_ + "/pop.filtered.marker",
                               self.target_dir_ + "/pop.filtered.detail.info")
        api.add_sg_marker_detail(task_id, marker_id, self.output_dir + "/pop.filtered.gtype.stat", self.option('type'))
        api.add_sg_distribution(task_id, marker_id, self.filtered_marker_path)
        # api.add_sg_heatmap(task_id, marker_id, self.output_dir + "/pop.filtered.matrix",
        #                    path=self.target_dir_ + "/pop.filtered.detail.info")  # 变re上传？？
        api.add_sg_marker_stat(task_id, marker_id, self.output_dir + "/marker_info.xls")
        api.add_sg_marker_child_stat(task_id, marker_id, self.output_dir + "/pop.filtered.indi.stat")
        self.logger.info("设置marker_filter的导表成功！")
        self.end()

    def run(self):
        self.get_file_path()
        self.marker_filter_run()
        super(MarkerFilterWorkflow, self).run()

    def end(self):
        """
        这里后面要重新定义下文件名字
        :return:
        """
        # self.set_output_file()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MarkerFilterWorkflow, self).end()

    def get_target_dir(self):
        """
        获取远程磁盘的路径
        :return:
        """
        # self.target_dir_ = self._sheet.output.split('://')[1]
        self.target_dir_ = self._sheet.output.rstrip('/')
