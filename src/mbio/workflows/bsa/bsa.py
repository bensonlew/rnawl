# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'

"""Bsa基础工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import re
import json
import time
import shutil


class BsaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        version = 1.0.0
        lasted modifed by HONGDONG 20180307
        """
        self._sheet = wsheet_object
        super(BsaWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'bsa_path', 'type': 'infile', 'format': 'bsa.bsa_path'},  # 产品线上传的数据的基础路径，该路径下还有其他的子路径，要进行下判断
            {"name": "marker_type", "type": "string", "default": "ALL"},  # 标记类型 默认为all（SNP+INDEL），或者只有SNP
            {"name": "group_type", "type": "string", "default": "F2"},  # 群体类型 默认F2（仅qtlseq.pl使用）
            {"name": "s_deep", "type": "int", "default": 10},  # 亲本测序深度，默认10X
            {"name": "b_deep", "type": "int", "default": 10},  # 混池测序深度，默认10X
            {"name": "slidingwin_type", "type": "string", "default": "distance"},  # 滑窗方式，variant_num/distance
            {"name": "slidingwin_value", "type": "string", "default": "1-5"},  # 滑窗策略，1M-5K
            {"name": "threshold", "type": "string", "default": "quantile"},  # 阈值确定， quantile or index 默认是quantile
            {"name": "threshold_value", "type": "string", "default": "0.999"},
            {"name": "sample_info", "type": "string"}  # 样本信息，{"p1":"","p2":"","b1":"","b2":""}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.base_info = self.api.api("bsa.api_base")
        self.index_calc = self.add_tool('bsa.index_calc')
        self.import_baseinfo = ""  # 导入基础信息的tool
        self.slidingwin_analysis = self.add_module("bsa.slidingwin_analysis")  # 滑窗计算的module
        self.region_analysis = self.add_module("bsa.region_analysis")  # 区域分析模块
        self.slidingwin_filter = self.add_tool("bsa.slidingwin_filter")  # 区域定位，按照分位值与index值进行过滤
        self.go_anno = self.add_tool("bsa.anno_analysis")
        self.kegg_anno = self.add_tool("bsa.anno_analysis")
        self.eggnog_anno = self.add_tool("bsa.anno_analysis")
        self.circos = self.add_tool("bsa.draw_circos")
        self.step.add_steps("baseinfo", "slidingwin", "filter", "region", "anno")  # 添加分析步骤
        self.final_vcf = ""
        self.pop_summary = ""
        self.ref_path = ""
        self.slide_type = ""
        self.wp = ''   # 亲本野生
        self.mp = ''   # 亲本突变
        self.wb = ''  # 混池野生
        self.mb = ''  # 混池突变
        self.win_ = ''   # the window size
        self.step_ = ''   # the step size
        self.anno_tools = []
        self.target_dir = ""  # 文件上传到磁盘对应的路径，bsa整个流程中，所有的不导表但是需要用于后
        # 面计算都保存磁盘的路径
        self.slidingwin_id = ""
        self.ref_chrlist = ""
        self.region_id = ""
        self.circos_path = ""  # 用于存储circos的图片
        self.index_path = ""    # 用于存储0x 0x的index的结果文件
        self.is_old = True
        self.target_vcf = ""
        self.target_fa = ""
        self.target_pop = ''
        self.target_chrlist = ''
        self.region = ''

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option("bsa_path"):
            raise OptionError("必须要有bsa_path参数", code="11500101")
        if not self.option("marker_type"):
            raise OptionError("必须设置marker_type", code="11500102")
        elif self.option("marker_type") not in ['ALL', 'SNP']:
            raise OptionError("标记类型%s不合法！", variables=(self.option("marker_type")), code="11500103")
        if self.option("slidingwin_type") not in ['variant_num', 'distance']:
            raise OptionError("滑窗方式%s不合法！", variables=(self.option("slidingwin_type")), code="11500104")
        if self.option("threshold") not in ['quantile', 'index']:
            raise OptionError("阈值确定%s不合法！", variables=(self.option("threshold")), code="11500105")
        return True

    def run_index_calc(self):
        """
        计算0x，0x的vcf文件的index值，这里的数据 用于基因详情页的展示
        :return:
        """
        self.index_calc.set_options({
            "pop_vcf": self.final_vcf,
            "wp": self.wp,
            "mp": self.mp,
            "wb": self.wb,
            "mb": self.mb,
            "pdep": 0,
            "bdep": 0,
            "popt": self.option("group_type"),
            "variant_type": self.option("marker_type")
        })
        self.index_calc.on("end", self.set_output, "index_calc")
        self.index_calc.run()

    def run_import_baseinfo(self):
        """
        导入产品线手动上传的文件
        :return:
        """
        if self.is_old:
            self.import_baseinfo = self.add_tool("bsa.import_baseinfo_old")
        else:
            self.import_baseinfo = self.add_tool("bsa.import_baseinfo")
        self.import_baseinfo.set_options({
            "bsa_path": self.option("bsa_path"),
            "task_id": self._sheet.id,
            "project_sn": self._sheet.project_sn,
            "member_id": self._sheet.member_id
        })
        self.import_baseinfo.on("start", self.set_step, {'start': self.step.baseinfo})
        self.import_baseinfo.on("end", self.set_step, {'end': self.step.baseinfo})
        self.import_baseinfo.on("end", self.run_slidingwin)
        self.import_baseinfo.run()

    def run_slidingwin(self):
        """
        进行划窗分析
        :return:
        """
        self.slidingwin_analysis.set_options({
            "pop_vcf":  self.final_vcf,
            "wp": self.wp,
            "mp": self.mp,
            "wb": self.wb,
            "mb": self.mb,
            "pdep": self.option("s_deep"),
            "bdep": self.option("b_deep"),
            "popt": self.option("group_type"),
            "method": 'bp' if self.option("slidingwin_type") == 'distance' else 'num',
            "win": self.win_,
            "step": self.step_,
            "variant_type": self.option("marker_type"),
            "abs": "0" if self.wp or self.mp else "1"   # 两个亲本都没有的话为1
        })
        self.slidingwin_analysis.on("start", self.set_step, {'start': self.step.slidingwin})
        self.slidingwin_analysis.on("end", self.set_output, "slidingwin")
        self.slidingwin_analysis.on("end", self.set_step, {'end': self.step.slidingwin})
        self.slidingwin_analysis.run()

    def run_slidingwinfilter(self):
        """
        进行过滤，进行区域定位
        :return:
        """
        self.slidingwin_filter.set_options({
            "slidingwin_file": self.slidingwin_analysis.output_dir + "/sliding-win.result",
            "region_type": self.option("threshold"),  # 阈值类型，分位值quantile/index值index
            "region_value": self.option("threshold_value"),  # 阈值
            "slid_file": self.slidingwin_analysis.output_dir + "/sliding-win.slid.result"
        })
        self.slidingwin_filter.on("start", self.set_step, {'start': self.step.filter})
        self.slidingwin_filter.on("end", self.set_output, "slidingwin_filter")
        self.slidingwin_filter.on("end", self.set_step, {'end': self.step.filter})
        self.slidingwin_filter.run()

    def run_region_analysis(self):
        """
        对定位后的区域进行分析，包含了variant，gene，vcf分析
        :return:
        """
        self.region_analysis.set_options({
            "i_c_result": self.slidingwin_analysis.output_dir + "/index-calc.result.index",
            "pop_summary": self.pop_summary,
            "p_f_vcf": self.final_vcf,
            "s_w_select": self.slidingwin_filter.output_dir + "/sliding-win.threshold.select",
            "wp": self.wp,  # 野生型亲本名称
            "mp": self.mp,  # 突变型亲本名称
            "wb": self.wb,  # 野生型混池名称
            "mb": self.mb  # 突变型混池名称
            # "step": self.step_  # 滑窗策略较小的数值
        })
        self.region_analysis.on("start", self.set_step, {'start': self.step.region})
        self.region_analysis.on("end", self.set_output, "region_analysis")
        self.region_analysis.on("end", self.set_step, {'end': self.step.region})
        self.region_analysis.run()

    def run_go_anno(self):
        """
        对定位后的区域进行基因注释分析 go anno
        :return:
        """
        self.go_anno.set_options({
            "anno_stat": self.region_analysis.output_dir + "/region_gene/region.threshold.gene.go.stat",
            "anno_type": "go"
        })
        self.go_anno.on("end", self.set_output, "go_anno")
        self.go_anno.run()

    def run_kegg_anno(self):
        """
        对定位后的区域进行基因注释分析 kegg anno
        :return:
        """
        self.kegg_anno.set_options({
            "anno_stat": self.region_analysis.output_dir + "/region_gene/region.threshold.gene.kegg.stat",
            "anno_type": "kegg"
        })
        self.kegg_anno.on("end", self.set_output, "kegg_anno")
        self.kegg_anno.on("start", self.set_step, {'start': self.step.anno})
        self.kegg_anno.on("end", self.set_step, {'end': self.step.anno})
        self.kegg_anno.run()

    def run_eggnog_anno(self):
        """
        对定位后的区域进行基因注释分析 eggnog anno
        :return:
        """
        self.eggnog_anno.set_options({
            "anno_stat": self.region_analysis.output_dir + "/region_gene/region.threshold.gene.eggnog.stat",
            "anno_type": "eggnog"
        })
        self.eggnog_anno.on("end", self.set_output, "eggnog_anno")
        self.eggnog_anno.run()

    def run_draw_circos(self):
        """
        用于生成circos的结果文件
        :return:
        """
        self.circos.set_options({
            "p_v_vcf": self.final_vcf,
            "s_w_result": self.slidingwin_analysis.output_dir + "/sliding-win.result",
            "chrlist": self.ref_chrlist,
            "pop_summary": self.pop_summary
        })
        self.circos.on("end", self.set_output, "circos")
        self.circos.run()

    def check_base_file(self):
        """
        用于检查并初始化pop.final.vcf与pop.summary
        :return:
        """
        vcf = os.path.join(self.option("bsa_path").prop['path'], "05.annovar/combine_variants/pop.final.vcf")
        summary = os.path.join(self.option("bsa_path").prop['path'], "05.annovar/anno_count/pop.summary")
        ref = os.path.join(self.option("bsa_path").prop['path'], "02.reference/ref.fa")
        if os.path.exists(vcf) and os.path.exists(summary) and os.path.exists(ref):
            self.is_old = False
        else:
            vcf = os.path.join(self.option("bsa_path").prop['path'], "05.annovar/pop.final.vcf.gz")
            summary = os.path.join(self.option("bsa_path").prop['path'], "05.annovar/pop.summary")
            ref = os.path.join(self.option("bsa_path").prop['path'], "02.reference/ref.fa.gz")
        ref_chrlist = os.path.join(self.option("bsa_path").prop['path'], "02.reference/ref.chrlist")
        self.logger.info(vcf)
        self.logger.info(summary)
        self.logger.info(ref)
        if not os.path.isfile(vcf) or not os.path.isfile(summary) or not os.path.isfile(ref):
            self.set_error("文件%s或%s或%s不存在！请检查上传的原始文件是否完整！", variables=(vcf, summary, ref), code="11500101")
            self.set_error("文件%s或%s或%s不存在！请检查上传的原始文件是否完整！", variables=(vcf, summary, ref), code="11500118")
        else:
            self.final_vcf, self.pop_summary, self.ref_path, self.ref_chrlist = \
                self.link_to_output(vcf, summary, ref, ref_chrlist)

    def link_to_output(self, vcf, summary, ref, ref_chrlist):
        """
        将接口后续需要的文件link到workflow的output中
        :return:
        """
        if not os.path.exists(self.output_dir + '/temp'):
            os.mkdir(self.output_dir + '/temp')
        if not self.is_old:
            new_vcf = self.output_dir + "/temp/pop.final.vcf"
            new_ref = self.output_dir + "/temp/ref.fa"
        else:
            new_vcf = self.output_dir + "/temp/pop.final.vcf.gz"
            new_ref = self.output_dir + "/temp/ref.fa.gz"
        new_summary = self.output_dir + "/temp/pop.summary"
        new_chrlist = self.output_dir + "/temp/ref.chrlist"
        if os.path.exists(new_vcf):
            os.remove(new_vcf)
        os.link(vcf, new_vcf)
        if os.path.exists(new_summary):
            os.remove(new_summary)
        os.link(summary, new_summary)
        if os.path.exists(new_ref):
            os.remove(new_ref)
        os.link(ref, new_ref)
        if os.path.exists(new_chrlist):
            os.remove(new_chrlist)
        os.link(ref_chrlist, new_chrlist)
        self.logger.info("链接vcf，summary，ref，ref_chrlist到output目录成功！")
        return new_vcf, new_summary, new_ref, new_chrlist

    def get_pid_bid(self):
        """
        用于重组亲本列表与混池列表  sample_info = {"p1":"","p2":"","b1":"","b2":""}
        p1,p2是亲本（p2是突变型），b1,b2是混池（b2是突变型）
        self.wp = ''   # 亲本野生
        self.mp = ''   # 亲本突变
        self.wb = ''  # 混池野生
        self.mb = ''  # 混池突变
        :return:
        """
        sample_info = json.loads(self.option("sample_info"))
        if sample_info['p1']:
            self.wp = sample_info['p1']
        if sample_info['p2']:
            self.mp = sample_info['p2']
        if sample_info['b1']:
            self.wb = sample_info['b1']
        if sample_info['b2']:
            self.mb = sample_info['b2']

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def move2outputdir(self, olddir, newname):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        :param olddir: 初始路径
        :param newname: 目标路径，可以自定义
        :return:
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code="11500119")
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        # self.logger.info(newfiles)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        end = time.time()
        duration = end - start
        self.logger.info("文件夹{}到{}移动耗时{}s".format(olddir, newdir, duration))

    def move_file(self, old_file, new_file):
        """
        递归移动文件或者文件到指定路径
        :param old_file: 初始路径
        :param new_file: 目的路径
        :return:
        """
        if os.path.isfile(old_file):
            os.link(old_file, new_file)
        else:
            os.mkdir(new_file)
            for file_ in os.listdir(old_file):
                file_path = os.path.join(old_file, file_)
                new_path = os.path.join(new_file, file_)
                self.move_file(file_path, new_path)

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] == "index_calc":
            self.move2outputdir(obj.output_dir, self.output_dir + "/temp/index_calc")
        if event['data'] == "slidingwin":
            self.move2outputdir(obj.output_dir, self.output_dir + "/slidingwin")
        if event['data'] == "slidingwin_filter":
            self.move2outputdir(obj.output_dir, self.output_dir + "/slidingwin_filter")
        if event['data'] == "region_analysis":
            self.move2outputdir(obj.output_dir, self.output_dir + "/region_analysis")
        if event['data'] == "kegg_anno":
            self.move2outputdir(obj.output_dir, self.output_dir + "/gene_annotation/kegg_anno")
        if event['data'] == "go_anno":
            self.move2outputdir(obj.output_dir, self.output_dir + "/gene_annotation/go_anno")
        if event['data'] == "eggnog_anno":
            self.move2outputdir(obj.output_dir, self.output_dir + "/gene_annotation/eggnog_anno")
        if event['data'] == "circos":
            self.move2outputdir(obj.output_dir, self.output_dir + "/temp/circos")

    def run(self):
        task_info = self.api.api('bsa.base_info')
        self.check_base_file()  # 设置vcf，pop.summary文件路径
        self.get_pid_bid()  # 初始化一下 获取亲本pid，与混池bid
        self.get_target_dir()  # 初始化一下 获取到远程磁盘的路径，用于保存对应文件的路径到主表中
        self.check_win_step()  # 初始化下 win 与step值
        task_info.add_sg_task(self._sheet.member_id, self._sheet.member_type, self._sheet.cmd_id,
                              self._sheet.project_sn, self._sheet.id, self.target_vcf,  self.target_pop, self.target_fa,
                              self.target_chrlist, self.region)
        self.logger.info("wp:{}, mp:{}, wb:{}, mb:{}".format(self.wp, self.mp, self.wb, self.mb))
        self.logger.info("该流程是运行的is_old值{}".format(self.is_old))
        self.slidingwin_analysis.on("end", self.run_slidingwinfilter)
        self.slidingwin_analysis.on("end", self.run_draw_circos)
        self.slidingwin_filter.on("end", self.run_region_analysis)
        self.region_analysis.on("end", self.run_go_anno)
        self.region_analysis.on("end", self.run_kegg_anno)
        self.region_analysis.on("end", self.run_eggnog_anno)
        self.on_rely([self.go_anno, self.kegg_anno, self.eggnog_anno, self.index_calc, self.circos], self.end)
        self.run_index_calc()
        self.run_import_baseinfo()
        super(BsaWorkflow, self).run()

    def run_api(self):
        """
        对所有结果进行导表操作
        :return:
        """
        self.import_specimen_type()
        self.import_index()
        self.import_slidingwin()
        self.import_slidingwinfilteranalysis()
        self.import_manhattan()
        self.import_circos()
        pass

    def import_specimen_type(self):
        """
        更新参与分析的样本，增加type为1，便于关联区域详情筛选功能样本展示正确
        """
        task_info = self.api.api('bsa.base_info')
        task_info.update_specimen_type(task_id=self._sheet.id, wp=self.wp, mp=self.mp, wb=self.wb, mb=self.mb)
        self.logger.info("更新sg_specimen表中参与分析的样本类型结束")

    def import_index(self):
        """
        更新0X 0X的index/variant文件到sg_task
        :return:
        """
        self.base_info.update_db_record("sg_task", {"task_id": self._sheet.id},
                                        {"index_path": self.index_path + "/index-calc.result.index",
                                        "variant_path": self.index_path + "/index-calc.result.variant"})

    def import_slidingwin(self):
        """
        导入滑窗的结果文件
        :return:
        """
        slidingwin = self.api.api("bsa.slidingwin")
        slidingwin_params = {"s_deep": self.option("s_deep"), "b_deep": self.option("b_deep"),
                             "sliding_type": self.option("slidingwin_type"),
                             "sliding_strategy": self.option("slidingwin_value"), "submit_location": "slidingwin",
                             "task_type": 2}
        self.slidingwin_id = self.base_info.add_main_table("sg_slidingwin", self._sheet.id, self._sheet.project_sn,
                                                           slidingwin_params, "origin_slidingwin", "标记筛选和分析主表",
                                                           self._sheet.member_id, 'false')
        calc_index_path = self.target_dir + "/slidingwin/index-calc.result.index"
        calc_variant_path = self.target_dir + "/slidingwin/index-calc.result.variant"
        slidingwin_result_path = self.target_dir + "/slidingwin/sliding-win.result"
        slidingwin_slid_result_path = self.target_dir + "/slidingwin/sliding-win.slid.result"
        slidingwin_stat_path = self.slidingwin_analysis.output_dir + "/index-calc.result.final.stat"
        self.base_info.update_db_record("sg_slidingwin", {"_id": self.slidingwin_id},
                                        {"wp": self.wp, "mp": self.mp, "wb": self.wb, "mb": self.mb,
                                         "calc_index_path": calc_index_path, "calc_variant_path": calc_variant_path,
                                         "slidingwin_result_path": slidingwin_result_path,
                                         "variant_type": ("ALL" if self.option("marker_type") == "all" else 'SNP'),
                                         "slid_result_path": slidingwin_slid_result_path, "status": "end"})
        slidingwin.add_sg_slidingwin_stat(slidingwin_stat_path, self.slidingwin_id)

    def import_slidingwinfilteranalysis(self):
        """
        导入到过滤后的滑窗结果,以及过滤后的定位区域的分析结果，以及注释结果
        :return:
        """
        pathway_dir = self._sheet.output.rstrip('/') + "/gene_annotation/kegg_anno/pathway_dir"
        bsa_region = self.api.api("bsa.bsa_region")
        bsa_anno = self.api.api("bsa.bsa_region_anno")
        region_variant_path = self.region_analysis.output_dir + "/region_variant/region.threshold.variant.total"
        region_gene_path = self.region_analysis.output_dir + "/region_gene/region.threshold.gene.total"
        region_vcf_path = self.region_analysis.output_dir + "/region_vcf/region.threshold.vcf.total"
        # 导入到过滤后的滑窗结果， 关联区域过滤主表
        region_params = {"slidingwin_id": str(self.slidingwin_id), "threshold": self.option("threshold"),
                         "threshold_value": self.option("threshold_value"), "submit_location": "region", "task_type": 2}
        self.region_id = bsa_region.add_sg_region(self._sheet.id, self._sheet.project_sn, self.slidingwin_id,
                                                  self.wp, self.mp, self.wb, self.mb, "origin_filter_analysis",
                                                  region_params)
        # 关联区域变异位点统计表
        start = time.time()
        bsa_region.add_sg_region_variant(self.region_id, region_gene_path, region_variant_path)
        # 关联区域变异位点详情表
        end2 = time.time()
        self.logger.info("导sg_region_variant花费时间:{}".format(end2 - start))
        bsa_region.add_sg_region_vcf(self.region_id, region_vcf_path)
        end = time.time()
        self.logger.info("导sg_region_vcf花费时间:{}".format(end - end2))
        # 关联区域基因型频率详情表
        bsa_region.add_sg_region_index(self.region_id, region_variant_path)
        end1 = time.time()
        self.logger.info("导sg_region_index花费时间:{}".format(end1 - end))
        # 关联区域基因注释主表
        anno_id = bsa_anno.add_sg_region_anno(self.region_id)
        # 所有的候选区域的基因注释信息，包含nr，uniport, go, kegg, eggnog
        bsa_anno.sg_region_anno_detail(anno_id, region_gene_path)
        bsa_anno.sg_region_anno_go_stat(anno_id, self.go_anno.output_dir + "/region.threshold.gene.go.final.stat.detail")
        bsa_anno.sg_region_anno_kegg_stat(anno_id,
                                          self.kegg_anno.output_dir + "/region.threshold.gene.kegg.final.stat.detail",
                                          pathway_dir)
        bsa_anno.sg_region_anno_eggnog_stat(
            anno_id, self.eggnog_anno.output_dir + "/region.threshold.gene.eggnog.final.stat.detail")

    def import_circos(self):
        """
        将circos的图导入到表中
        :return:
        """
        start = time.time()
        # index = self.circos_path.find("/rerewrweset")
        # circos_path_ = self.circos_path[index:]
        slidingwin = self.api.api("bsa.slidingwin")
        slidingwin.add_sg_circos(self.slidingwin_id, self.circos_path + "/gene.num.csv",
                                 self.circos_path + "/circos.chrlist", self.output_dir + "/temp/circos/circos.chrlist",
                                 self.circos_path + "/snp.win.csv", self.circos_path + "/indel.win.csv",
                                 self.circos_path + "/sliding.win.csv")
        end = time.time()
        self.logger.info("导sg_circos花费时间:{}".format(end - start))

    def import_manhattan(self):
        """
        导入曼哈顿图
        :return:
        """
        start = time.time()
        slidingwin = self.api.api("bsa.slidingwin")
        slidingwin_result_path = self.slidingwin_analysis.output_dir + "/sliding-win.result"
        slidingwin.add_sg_manhattan(self.slidingwin_id, slidingwin_result_path, self.ref_chrlist)
        end = time.time()
        self.logger.info("导sg_manhattan花费时间:{}".format(end - start))
        bsa_region = self.api.api("bsa.bsa_region")
        bsa_region.update_sg_region_quantile(self.region_id, self.slidingwin_filter.output_dir + "/quantile.index")

    def get_target_dir(self):
        """
        获取远程磁盘的路径
        :return:
        """
        # self.target_dir = self._sheet.output.strip().split(':')[1]
        self.target_dir = self._sheet.output.rstrip('/')
        self.region = self._sheet.output.strip().split('://')[0]
        self.circos_path = self._sheet.output.rstrip('/') + "/temp/circos"
        self.index_path = self.target_dir + "/temp/index_calc"
        if not self.is_old:
            self.target_vcf = self.target_dir + "/temp/pop.final.vcf"
            self.target_fa = self.target_dir + "/temp/ref.fa"
        else:
            self.target_vcf = self.target_dir + "/temp/pop.final.vcf.gz"
            self.target_fa = self.target_dir + "/temp/ref.fa.gz"
        self.target_pop = self.target_dir + "/temp/pop.summary"
        self.target_chrlist = self.target_dir + "/temp/ref.chrlist"

    def check_win_step(self):
        """
        这一步放到这里来检查，之前放在check option中 self.win_与self.step_赋值不了，暂时不知道原因
        :return:
        """
        if self.option("slidingwin_type") == "distance":
            # self.slide_type = "bp"
            m = re.match(r"(.*)-(.*)", self.option("slidingwin_value"))
            if m:
                self.win_ = float(m.group(1)) * 1000000
                self.step_ = float(m.group(2)) * 1000
                if self.win_ < self.step_:
                    self.set_error("滑窗策略M前数值%s需大于K前数值%s!", variables=(self.win_, self.step_), code="11500102")
                    self.set_error("滑窗策略M前数值%s需大于K前数值%s!" , variables=(self.win_, self.step_), code="11500120")
                self.win_ = str(self.win_)
                self.step_ = str(self.step_)
            else:
                self.set_error("滑窗策略%s必须是'数值M-数值K'的形式!", variables=(self.option("slidingwin_value")), code="11500103")
                self.set_error("滑窗策略%s必须是'数值M-数值K'的形式!" , variables=( self.option("slidingwin_value")), code="11500121")
        else:
            # self.slide_type = "num"
            m = re.match(r"(\d+)-(\d+)", self.option("slidingwin_value"))
            if m:
                self.win_ = int(m.group(1))
                self.step_ = int(m.group(2))
                if self.win_ < self.step_:
                    self.set_error("滑窗策略中划线前数值%s需大于中划线后数值%s!", variables=(self.win_, self.step_), code="11500104")
                    self.set_error("滑窗策略中划线前数值%s需大于中划线后数值%s!" , variables=(self.win_, self.step_), code="11500122")
                self.win_ = str(self.win_)
                self.step_ = str(self.step_)
            else:
                self.set_error("滑窗策略%s必须是'正整数-正整数'的形式!", variables=(self.option("slidingwin_value")), code="11500105")
                self.set_error("滑窗策略%s必须是'正整数-正整数'的形式!" , variables=( self.option("slidingwin_value")), code="11500123")

    def remove_file(self):
        """
        删除不需要上传到磁盘的文件
        :return:
        """
        rm_list = list()
        rm_list.append(self.output_dir + "/index_calc/index-calc.result.final.stat")
        rm_list.append(self.output_dir + "/index_calc/index-calc.result.variant")
        rm_list.append(self.output_dir + "/slidingwin/index-calc.result.variant")
        rm_list.append(self.output_dir + "/slidingwin_filter/quantile.index")
        for files in rm_list:
            if os.path.isfile(files):
                os.remove(files)
                self.logger.info("删除文件{}成功！".format(files))
            else:
                self.logger.info("文件{}不存在，不用删除！".format(files))

    def send_files(self):
        self.remove_file()
        repaths = [
            [".", "", "基础分析结果文件夹",0,"310001"],
            ["temp", "", "临时文件目录",0,"310002"],
            ["temp/pop.final.vcf", "", "vcf合并文件",0,"310003"],
            ["temp/pop.summary", "", "基因功能注释表",0,"310004"],
            ["temp/ref.fa", "", "参考基因组fa序列",0,"310005"],
            ["temp/ref.chrlist", "", "参考基因组染色体列表文件",0,"310006"],
            ["temp/circos", "", "Circos图数据目录",0,"310007"],
            ["temp/index_calc", "", "Index数据目录",0,"310008"],
            ["gene_annotation", "", "关联区域基因注释文件目录",0,"310009"],
            ["gene_annotation/eggnog_anno", "", "EGGNOG注释文件目录",0,"310010"],
            ["gene_annotation/go_anno", "", "GO注释文件目录",0,"310011"],
            ["gene_annotation/kegg_anno", "", "KEGG注释文件目录",0,"310012"],
            ["gene_annotation/eggnog_anno/region.threshold.gene.eggnog.final.stat.detail", "", "EGGNOG分类统计文件",0,"310013"],
            ["gene_annotation/go_anno/region.threshold.gene.go.final.stat.detail", "", "GO分类统计文件",0,"310014"],
            ["gene_annotation/kegg_anno/region.threshold.gene.kegg.final.stat.detail", "", "KEGG分类统计文件",0,"310015"],
            ["gene_annotation/kegg_anno/pathway_dir", "", "KEGG  Pathway  文件目录",0,"310016"],
            ["region_analysis", "", "关联区域分析结果文件目录",0,"310017"],
            ["region_analysis/region_variant", "", "关联区域基因型频率文件目录",0,"310018"],
            ["region_analysis/region_variant/region.threshold.variant.total", "", "关联区域全部基因型频率文件",0,"310019"],
            ["region_analysis/region_variant/region.threshold.variant.eff", "", "关联区域功能区基因型频率文件",0,"310020"],
            ["region_analysis/region_vcf", "", "关联区域变异位点文件目录",0,"310021"],
            ["region_analysis/region_vcf/region.threshold.vcf.total", "", "关联区域全部变异位点文件",0,"310022"],
            ["region_analysis/region_vcf/region.threshold.vcf.eff", "", "关联区域功能区变异位点文件",0,"310023"],
            ["region_analysis/region_gene", "", "关联区域基因注释详情文件目录",0,"310024"],
            ["region_analysis/region_gene/region.threshold.gene.eff", "", "关联区域功能区基因注释详情文件",0,"310025"],
            ["region_analysis/region_gene/region.threshold.gene.eggnog.stat", "", "关联区域基因EGGNOG分类详情文件",0,"310026"],
            ["region_analysis/region_gene/region.threshold.gene.go.stat", "", "关联区域基因GO分类详情文件",0,"310027"],
            ["region_analysis/region_gene/region.threshold.gene.kegg.stat", "", "关联区域基因KEGG分类详情文件",0,"310028"],
            ["region_analysis/region_gene/region.threshold.gene.total", "", "关联区域全部基因注释详情文件",0,"310029"],
            ["slidingwin", "", "标记筛选和分析结果文件目录",0,"310030"],
            ["slidingwin/index-calc.result.final.stat", "", "标记筛选统计结果文件",0,"310031"],
            ["slidingwin/index-calc.result.index", "", "标记筛选基因型频率详情文件",0,"310032"],
            ["slidingwin/sliding-win.slid.result", "", "",0,"310033"],
            ["slidingwin/sliding-win.result", "", "全部标记基因型频率文件",0,"310034"],
            ["slidingwin_filter", "", "关联区域定位结果文件目录",0,"310035"],
            ["slidingwin_filter/sliding-win.threshold.select", "", "关联区域定位结果文件",0,"310036"]
        ]
        regexps = [
            [r"gene_annotatio/kegg_anno/pathway_dir/*\.png$", "", "Pathway图文件",0,"310037"],
            [r"gene_annotatio/kegg_anno/pathway_dir/*\.pdf$", "", "Pathway图文件",0, "310038"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        # for i in self.get_upload_files():
        #     self.logger.info('upload file:{}'.format(str(i)))

    def end(self):
        self.send_files()
        self.logger.info("上传文件目录完成")
        self.run_api()
        self.logger.info("导表完成")
        super(BsaWorkflow, self).end()
