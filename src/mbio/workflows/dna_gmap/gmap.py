# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'

"""GMAP基础工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId
import os
import re
import gevent
import time
import shutil
import json


class GmapWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        version = 1.0.0
        lasted modifed by HONGDONG 20180713
        """
        self._sheet = wsheet_object
        super(GmapWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'wgs_path', 'type': 'infile', 'format': 'bsa.bsa_path'},  # 输入wgs分析的路径
            {"name": "pid", "type": "string"},  # 页面选择亲本-父本
            {"name": "mid", "type": "string"},  # 页面选择亲本-母本
            {"name": 'popt', 'type': 'string'},  # 群体类型 F1/F2/DH/BC1/RIL
            {"name": "daishu", "type": "int", "default": 10},   # 当RIL群体的时候要设置子代数
            {"name": "chr_num", "type": "int"},  # 染色体个数
            {"name": "marker_table", "type": "infile", "format": "dna_gmap.marker"},  # 自定义标记列表
            {"name": "signif", "type": "float", "default": 0.05},  # 偏分离P值
            {"name": "miss_tatio", "type": "float", "default": 30},  # 缺失率
            {"name": "pdep", "type": "string", "default": "10_"},  # 亲本深度
            {"name": "odep", "type": "string", "default": "1_"},  # 子代深度
            {"name": "marker_type", "type": "string", "default": "snp,indel"},  # 标记类型选择
            {"name": "w_bin", "type": "int"},  # windows大小 2Mb   0.1~20MB
            {"name": "w_step", "type": 'int'},  # windows步长 100 kb
            {"name": "group_type", "type": "string", "default": "ref"},  # 分群方法 mlod+ref
            {"name": "lod_start", "type": "int", "default": 4},  # LOD 开始值
            {"name": "lod_end", "type": "int", "default": 20},   # LOD 结束值
            {"name": "min_group", "type": "int", "default": 20},  # 最小群MARKER数
            {"name": "max_group", "type": "int", "default": 500},   # 最大群Marker数
            {"name": "use_marker_table", "type": "string", "default": "false"}  # 上传的自定义标记列表做Bin
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.json_path = self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome/"  # 参考组配置文件
        self.wgs_api = self.api.api("wgs.api_base")
        self.api_base = self.api.api("dna_gmap.api_base")   # 导表的基础模块
        self.gmap_base = self.api.api("dna_gmap.gmap_base")
        self.linkage_groupping = self.add_module('dna_gmap.linkage_groupping')  # 连锁分群
        self.binner_calculate = self.add_tool("dna_gmap.binner_calculate")  # bin marker计算
        self.vcf_marker = self.add_tool("dna_gmap.vcf_marker")  # vcf文件转为marker
        self.marker_filters = []  # 遗传标记过滤
        self.mlod_tree = self.add_tool("dna_gmap.mlod_tree")  # 生成树图
        self.step.add_steps("vcf_marker", "marker_filter", "binner_calculate", "linkage_groupping")  # 添加分析步骤
        self.target_dir = ""  # 文件上传到磁盘对应的路径，所有的不导表但是需要用于后面计算都保存磁盘的路径
        self.wgs_path = ""
        self.genome_version_id = ''
        self.is_bin_marker = False
        self.has_wgs_id = False
        self.wgs_task_id = ''   # 后面导表的需要，所以要把输入的wgs的taskid给获取到
        self.software = ""  # 获取wgs项目中进行snp与indel使用的方法
        self.is_cnv = False  # 判断是否做了cnv
        self.is_sv = False  # 判断是否做了sv
        self.import_wgs = True   # 用于判断是不是要导入wgs的数据
        self.marker_id = ""
        self.is_total_lg = True  # 默认分群成功！
        self.has_child_list = True   # 默认是输入子代列表
        self.generationid = ""  # 子代样本id
        self.popt = ""
        self.total_chrlist = ""
        self.pop_final_vcf = ""  # 设置pop final vcf 绝对路径
        self.ref_chrlist = ""
        self.min_group = 0
        self.max_group = 0
        self.sample_num = 0
        self.region = ""

    def check_options(self):
        """
        检查参数设置--样本信息与分组方案要重新写下
        """
        if not self.option("wgs_path").is_set:
            raise OptionError("必须要输入wgs结算的结果文件路径！", code="14800101")
        if not self.option("pid") or not self.option("mid"):
            raise OptionError("必须要选择父本与母本", code="14800102")
        if self.option("popt"):
            if self.option('popt').lower() not in ["f1", "f2", "dh", "bc1", "ril"]:
                raise OptionError("群体类型%s不正确， 必须为F1/F2/DH/BC1/RIL！", variables=(self.option('popt')), code="14800103")
            if self.option('popt').lower() == 'ril' and self.option("daishu") <= 2:
                raise OptionError("RIL群体的代数必须要大于2", code="14800104")
        else:
            raise OptionError("必须要输入群体类型", code="14800105")
        if not self.option("chr_num"):
            raise OptionError("必须要输入染色体的个数！", code="14800106")
        if self.option("marker_type"):
            if self.option("marker_type").lower() not in ['indel', "snp", "snp,indel", "indel,snp"]:
                raise OptionError("输入的群体类型%s不合法！", variables=(self.option("marker_type")), code="14800107")
        else:
            raise OptionError("必须要输入群体类型参数！", code="14800108")
        if not self.option("group_type"):
            raise OptionError("必须输入group_type参数！", code="14800109")
        else:
            if self.option("group_type") not in ["ref", "mlod"]:
                raise OptionError("group_type参数必须为ref 或者 mlod!", code="14800110")
        if self.option("min_group") > self.option("max_group"):
            raise OptionError("最小makrer数必须要小于最大marker数", code="14800111")
        return True

    def run_vcf_marker(self):
        """
        vcf文件转换成marker
        :return:
        """
        self.vcf_marker.set_options({
            "vcf_path": self.pop_final_vcf,
            "pid": self.option("pid"),
            "mid": self.option("mid")
        })
        self.vcf_marker.on("end", self.set_output, "vcf_marker")
        self.vcf_marker.on("start", self.set_step, {'start': self.step.vcf_marker})
        self.vcf_marker.on("end", self.set_step, {'end': self.step.vcf_marker})
        self.vcf_marker.run()

    def run_marker_filter(self):
        if self.option('popt').lower() == 'ril':
            self.popt = "RIL{}".format(self.option('daishu'))
        else:
            self.popt = self.option('popt')
        ana_dict = [{"type": "ALL", "pdep": "0", "odep": "0", "miss_tatio": 1.0, "signif": 0.0},
                    {"type": "ALL" if self.option("marker_type").lower() == "snp,indel" else self.option("marker_type"),
                     "pdep": self.option("pdep"), "odep": self.option("odep"),
                     "miss_tatio": self.option("miss_tatio") / 100, "signif": self.option("signif")}]
        for m in ana_dict:
            self.logger.info("{}".format(m))
            marker = self.add_tool("dna_gmap.marker_filter")
            opt = {
                "vcf": self.pop_final_vcf,
                "detail_info": self.vcf_marker.option("primary_marker"),
                "type": m['type'],
                "pdep": m['pdep'],
                "odep": m['odep'],
                "popt": self.popt,
                "miss_tatio": m['miss_tatio'],
                "signif": m['signif']
            }
            if self.has_child_list:
                opt.update({
                    "child_list": self.generationid
                })
            if self.option("marker_table").is_set:
                opt.update({
                    "marker_upload": self.option("marker_table").prop['path']
                })
            marker.set_options(opt)
            self.logger.info("has_child_list:{}; generationid:{}".format(self.has_child_list, self.generationid))
            self.marker_filters.append(marker)
        self.on_rely(self.marker_filters, self.check_is_binmarker)
        self.marker_filters[0].on("start", self.set_step, {'start': self.step.marker_filter})
        self.marker_filters[0].on("end", self.set_output, "marker_filter_0_0")
        self.marker_filters[1].on("end", self.set_step, {'end': self.step.marker_filter})
        self.marker_filters[1].on("end", self.set_output, "marker_filter")
        for tool in self.marker_filters:
            gevent.sleep(1)
            tool.run()

    def check_is_binmarker(self):
        """
        用于自动检测，是否进行binmaker计算， 一条染色体的平均marker数是否是样品总个数的6倍
        :return:
        """
        file_path = self.marker_filters[1].option("filter_marker").prop['path']
        markers = self.get_file_len(file_path)
        if markers < 2000:
            raise Exception("过滤后的标记数量小于2000条，无法进行后续的分析，流程自动停止！")
        chrs = self.get_file_len(self.ref_chrlist)
        if markers / float(chrs) > self.sample_num * 6:
            self.logger.info("一条染色体的平均marker数是否是样品总个数的6倍，将进行bin marker计算！")
            self.is_bin_marker = True
            self.run_binner_calculate()
        else:
            self.linkage_groupping_run()

    def run_binner_calculate(self):
        opt = {
            "pop_final_vcf": self.pop_final_vcf,
            "genotype_matrix": self.marker_filters[1].option("filter_marker"),
            "pop_type": self.option("popt"),
            "window_size": 500,   # 0.5M
            "window_step": 100    # 100k
        }
        self.binner_calculate.set_options(opt)
        self.binner_calculate.on("end", self.set_output, "binner_calculate")
        self.binner_calculate.on("end", self.set_step, {'end': self.step.binner_calculate})
        self.binner_calculate.on("start", self.set_step, {'start': self.step.binner_calculate})
        self.binner_calculate.on("end", self.linkage_groupping_run)
        self.binner_calculate.run()

    def linkage_groupping_run(self):
        file = self.binner_calculate.option("total_bin_marker").prop['path'] if self.is_bin_marker else \
            self.marker_filters[1].option("filter_marker").prop['path']
        markers = self.get_file_len(file)
        if markers < 2000:
            raise Exception("过滤后的标记数量小于2000条，无法进行后续的分析，流程自动停止！")
        self.min_group, self.max_group = self.check_markers()
        self.logger.info("min_group:{}, max_group:{}".format(self.min_group, self.max_group))
        opt = {
            'marker': self.binner_calculate.option("total_bin_marker") if self.is_bin_marker else
            self.marker_filters[1].option("filter_marker"),
            'ref_chrlist': self.ref_chrlist,
            'poptype': self.option("popt"),
            'bin': "yes" if self.is_bin_marker else "no",
            'is_ref': "true",
            'group_type': self.option('group_type'),
            "pop_marker_detail": self.marker_filters[1].output_dir + "/pop.filtered.detail.info"
        }
        if self.option('group_type') == "mlod":
            opt.update({
                'chr_num': self.option('chr_num'),
                'start_lod': self.option('lod_start'),
                'end_lod': self.option('lod_end'),
                'min_group': self.min_group,
                'max_group': 500
                # 'max_group': int(self.max_group - 1 / float(5))
            })
        self.linkage_groupping.set_options(opt)
        self.linkage_groupping.on('end', self.set_output, 'linkage_groupping')
        self.linkage_groupping.on('end', self.run_groupping)
        self.linkage_groupping.on("start", self.set_step, {'start': self.step.linkage_groupping})
        self.linkage_groupping.run()

    def run_groupping(self):
        """
        start_lod=2， end_lod=50， min_group=20， max_group=实际marker个数
        :return:
        """
        self.mlod_tree.set_options({
            "total_mlod": self.linkage_groupping.option("total_mlod"),
            "marker": self.binner_calculate.option("total_bin_marker") if self.is_bin_marker else
            self.marker_filters[1].option("filter_marker"),
            'chr_num': self.option("chr_num"),
            "start_lod": 2,
            "end_lod": 50
        })
        self.mlod_tree.on("end", self.set_output, "tree_data")
        self.mlod_tree.on("end", self.set_step, {'end': self.step.linkage_groupping})
        # self.mlod_tree.on("end", self.end)
        self.mlod_tree.run()

    def check_markers(self):
        if self.is_bin_marker:
            file_path = self.binner_calculate.option("total_bin_marker").prop['path']
        else:
            file_path = self.marker_filters[1].option("filter_marker").prop['path']
        lens = len(open(file_path, 'rU').readlines())
        if lens < 22:
            return 4, lens
        else:
            return 20, lens

    def run(self):
        self.get_target_dir()  # 初始化一下 获取到远程磁盘的路径，用于保对应文件的路径到主表中
        self.get_wgs_taskid()   # 检查 输入的结果是不是wgs的项目
        self.check_chr_ref()  # 检查是不是全部sca水平，并选择ref进行分组
        self.has_child_list, self.generationid, self.sample_num = self.gmap_base.get_child_id(self._sheet.id)
        self.gmap_base.add_sg_task(self._sheet.member_id, self._sheet.member_type, self._sheet.cmd_id)
        self.vcf_marker.on('end', self.run_marker_filter)
        self.run_vcf_marker()
        super(GmapWorkflow, self).run()
        # self.start_listener()
        # self.fire("start")
        # self.end()

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
            raise Exception('需要移动到output目录的文件夹不存在。')
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
        if event['data'] == "vcf_marker":
            self.move2outputdir(obj.output_dir, self.output_dir + "/01.vcf_convert")
        if event['data'] == "marker_filter":
            self.move2outputdir(obj.output_dir, self.output_dir + "/02.marker_filter")
        if event['data'] == "binner_calculate":
            self.move2outputdir(obj.output_dir, self.output_dir + "/03.binner")
        if event['data'] == "linkage_groupping":
            self.move2outputdir(obj.output_dir, self.output_dir + "/04.grouping")
        if event['data'] == "tree_data":
            self.move2outputdir(obj.output_dir, self.output_dir + "/temp/tree_data")
            self.end()
        if event['data'] == "marker_filter_0_0":
            self.move2outputdir(obj.output_dir, self.output_dir + "/temp/marker_filter")

    def send_files(self):

        # self.remove_file()
        repaths = [
            [".", "", "基础分析结果文件夹",0,"330001"],
            ["01.vcf_convert", "", "vcf转化文件目录",0,"330002"],
            ["01.vcf_convert/pop.primary.marker", "", "",0,"330003"],
            ["01.vcf_convert/pop.primary.marker.log", "", "",0,"330004"],
            ["02.marker_filter/marker_info.xls", "", "marker过滤结果文件目录",0,"330005"],
            ["02.marker_filter/pop.filtered.detail.info", "", "",0,"330006"],
            ["02.marker_filter/pop.filtered.gtype.stat", "", "",0,"330007"],
            ["02.marker_filter/pop.filtered.indi.stat", "", "",0,"330008"],
            ["02.marker_filter/pop.filtered.marker", "", "",0,"330009"],
            ["02.marker_filter/pop.filtered.markerStat", "", "",0,"330010"],
            ["02.marker_filter/pop.filtered.matrix", "", "",0,"330011"],
            ["03.binner", "", "bin marker结果文件目录",0,"330012"],
            ["03.binner/bin_info.xls", "", "",0,"330013"],
            ["03.binner/bin_pos.xls", "", "",0,"330014"],
            ["03.binner/bin_stat.xls", "", "",0,"330015"],
            ["03.binner/Total.bin.marker", "", "",0,"330016"],
            ["04.grouping", "", "连锁分群结果文件目录",0,"330017"],
            ["04.grouping/evalutaion", "", "",0,"330018"],
            ["04.grouping/evalutaion/total.csv", "", "",0,"330019"],
            ["04.grouping/evalutaion/total.loc", "", "",0,"330020"],
            ["04.grouping/evalutaion/total.map", "", "",0,"330021"],
            ["04.grouping/evalutaion/total.mapstat", "", "",0,"330022"],
            ["04.grouping/evalutaion/total.marker", "", "",0,"330023"],
            ["04.grouping/evalutaion/total.marker.info", "", "",0,"330024"],
            ["04.grouping/evalutaion/total.marker.info.seg.region", "", "",0,"330025"],
            ["04.grouping/evalutaion/total.phy.spearman.xls", "", "",0,"330026"],
            ["04.grouping/groupping", "", "",0,"330027"],
            ["04.grouping/groupping/Total.lg", "", "",0,"330028"],
            ["04.grouping/marker_order", "", "",0,"330029"],
            ["04.grouping/marker_order/map_circle1", "", "",0,"330030"],
            ["04.grouping/marker_order/map_circle2", "", "",0,"330031"],
            ["04.grouping/marker_order/map_circle3", "", "",0,"330032"],
            ["04.grouping/mlod_calc", "", "",0,"330033"],
            ["04.grouping/mlod_calc/Total.mlod", "", "",0,"330034"],
            ["04.grouping/mlod_calc/mlod_calc_dir", "", "",0,"330035"],
            ["04.grouping/split_markers", "", "",0,"330036"]
        ]
        regexps = [
            [r"04.grouping/split_markers/.*\.loc", "", "correct.loc文件",0,"330037"],
            [r"04.grouping/split_markers/.*\.detail", "", "correct.loc详细信息文件",0,"330038"],
            [r"04.grouping/mlod_calc/mlod_calc_dir/.*\.mlod", "", "mlod文件",0,"330039"],
            [r"04.grouping/marker_order/map_circle1/.*\.marker$", "", "第一次排图marker结果文件",0,"330040"],
            [r"04.grouping/marker_order/map_circle1/.*\.marker\.detail", "", "第一次排图marker.detail结果文件",0,"330041"],
            [r"04.grouping/marker_order/map_circle1/.*\.out", "", "第一次排图out结果文件",0,"330042"],
            [r"04.grouping/marker_order/map_circle2/.*\.marker$", "", "第二次排图marker结果文件",0,"330043"],
            [r"04.grouping/marker_order/map_circle2/.*\.marker\.detail", "", "第二次排图marker.detail结果文件",0,"330044"],
            [r"04.grouping/marker_order/map_circle2/.*\.out", "", "第二次排图out结果文件",0,"330045"],
            [r"04.grouping/marker_order/map_circle3/.*\.marker$", "", "第三次排图marker结果文件",0,"330046"],
            [r"04.grouping/marker_order/map_circle3/.*\.marker\.detail", "", "第三次排图marker.detail结果文件",0,"330047"],
            [r"04.grouping/marker_order/map_circle3/.*\.out", "", "第三次排图out结果文件",0,"330048"],
            [r"04.grouping/evalutaion/.*\.mapstat", "", "map统计文件",0,"330049"],
            [r"04.grouping/evalutaion/total\..*\.info$", "", "图谱标记详情结果文件",0,"330050"],
            [r"04.grouping/evalutaion/total\..*\.info\.seg\.region", "", "图谱标记详情区域结果文件",0,"330051"],
            [r"04.grouping/evalutaion/total\..*\.phy\.spearman\.xls", "", "相关性结果文件",0,"330052"],
            [r"04.grouping/evalutaion/total\..*\.loc$", "", "图谱评估loc文件",0,"330053"],
            [r"04.grouping/evalutaion/total\..*\.map$", "", "图谱评估map文件",0,"330054"],
            [r"04.grouping/evalutaion/total\..*\.phase$", "", "图谱评估phase文件",0,"330055"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        # for i in self.get_upload_files():
        #     self.logger.info('upload file:{}'.format(str(i)))

    def run_api(self):
        """
        对所有结果进行导表操作
        :return:
        """
        self.stop_timeout_check()
        self.import_specimen_info()
        if self.import_wgs:
            self.import_mapping_results()
            self.import_snp_results()
            self.import_indel_results()
            # if self.is_cnv:
            #     self.import_cnv_results()
            # if self.is_sv:
            #     self.import_sv_results()
            self.import_annovar_results()
            # self.import_ssr_results()
            self.import_genomic_view()
            self.import_files_paths()
        self.add_all_marker_filter()
        if self.is_bin_marker:
            self.import_bin_marker()
        if os.path.exists(self.output_dir + "/04.grouping/groupping/Total.lg"):
            self.logger.info("生成了total.lg文件后面将进行连锁分群导表！")
            self.import_linkage_grouping()
        else:
            self.logger.info("没有生成total.lg文件后面不进行连锁分群导表！")

    def get_wgs_taskid(self):
        """
        根据输入的文件夹路径获取wgs的路径，files/m_188/188_5af3dc68d1df4/tsanger_30196/workflow_results
        如果能够获取到wgs task id则直接从wgs项目中拷贝相关的样本信息，否则就要给出样本信息列表，具体参考fastq.txt,
        :return:
        """
        task_id = os.path.basename(os.path.dirname(self.wgs_path.rstrip('/')))
        self.logger.info("WGS的task_id为：{}".format(task_id))
        # if not re.match("tsg_|tsanger_|sanger_|i-sanger_.*", task_id):
        #     raise Exception("task id is not correct!")
        result = self.wgs_api.col_find_one("sg_task", {"task_id": task_id})
        if result:
            self.has_wgs_id = True
            self.wgs_task_id = task_id
            self.genome_version_id = result['genome_version_id']
            params = self.wgs_api.col_find_one("sg_snp_call", {"task_id": task_id})
            if params:
                self.software = params['params']
            else:
                self.software = '{\"\":\"\"}'
            cnv = self.wgs_api.col_find_one("sg_cnv_call", {"task_id": task_id})
            if cnv:
                self.is_cnv = True
            sv = self.wgs_api.col_find_one("sg_sv_call", {"task_id": task_id})
            if sv:
                self.is_sv = True
            self.logger.info("设置软件名称成功！")
        else:
            self.logger.info("{}在wgs项目的sg_task中没有找到对应信息！")

    def import_specimen_info(self):
        self.logger.info("开始进行样本信息导表")
        wgs_base = self.api.api("wgs.wgs_base")
        if self.has_wgs_id:
            sample_info = wgs_base.get_specimen_other_info(self.wgs_task_id)
            self.gmap_base.add_sg_specimen(sample_info, self.option("pid"), self.option("mid"))
        else:
            if os.path.exists(self.option("wgs_path").prop['path'] + "/fastq.txt"):
                self.gmap_base.add_sg_specimen_other(self.option("wgs_path").prop['path'] + "/fastq.txt",
                                                     self._sheet.id)
                self.api.del_api("wgs.wgs_base")  # 删除api引用
                wgs_base = self.api.api("wgs.wgs_base")
                wgs_base.project_type = "dna_gmap"
                wgs_base.add_sg_specimen()
            else:
                self.import_wgs = False
                self.logger.info("没有wgs的样本信息，不进行前面的导表！")
                return
        self.logger.info("样本信息导表成功！")
        self.api.del_api("wgs.wgs_base")  # 删除api引用
        wgs_base = self.api.api("wgs.wgs_base")
        wgs_base.project_type = "dna_gmap"  # 切换为dna_gmap
        # wgs_base.update_clean_path(self.option("wgs_path").prop['path'] + "/01.fastq_qc/clean_data")
        self.logger.info("更新clean fastq成功！")
        self.logger.info("开始进行样本质控信息导表")
        wgs_base.add_sg_specimen_qc(self.option("wgs_path").prop['path'] + "/01.fastq_qc/rawdata_qc/qc.xls",
                                    self.option("wgs_path").prop['path'] + "/01.fastq_qc/cleandata_qc/qc.xls")
        self.logger.info("样本质控信息导表成功！")
        self.logger.info("开始进行碱基含量分布导表")
        wgs_base.add_qc_atgc_curve(self.option("wgs_path").prop['path'] + "/01.fastq_qc/rawdata_qc/atgc", "raw_reads")
        wgs_base.add_qc_atgc_curve(self.option("wgs_path").prop['path'] + "/01.fastq_qc/cleandata_qc/atgc",
                                   "clean_reads")
        self.logger.info("碱基含量分布导表成功！")
        self.logger.info("开始进行碱基错误率分布导表")
        wgs_base.add_qc_qual_curve(self.option("wgs_path").prop['path'] + "/01.fastq_qc/rawdata_qc/qual",
                                   "error_raw_reads")
        wgs_base.add_qc_qual_curve(self.option("wgs_path").prop['path'] + "/01.fastq_qc/cleandata_qc/qual",
                                   "error_clean_reads")
        self.logger.info("碱基错误率分布导表成功！")

    def import_mapping_results(self):
        self.logger.info("开始进行比对结果导表")
        self.api.del_api("wgs.wgs_base")  # 删除api引用
        wgs_base = self.api.api("wgs.wgs_base")
        wgs_base.project_type = "dna_gmap"
        mapping_id = wgs_base.add_sg_mapping()
        wgs_base.add_sg_mapping_detail(mapping_id, self.option("wgs_path").prop['path'] +
                                       "/03.map_stat/result.stat/Total.mapped.detail.xls")
        self.logger.info("比对结果导表成功！")
        self.logger.info("开始进行插入片段分布导表")
        wgs_base.insert_sg_mapping_curve(self.option("wgs_path").prop['path'] + "/03.map_stat/insert", mapping_id,
                                         "insert")
        self.logger.info("插入片段分布导表成功！")
        self.logger.info("开始进行测序深度分布导表")
        wgs_base.insert_sg_mapping_curve(self.option("wgs_path").prop['path'] + "/03.map_stat/depth", mapping_id, "dep")
        self.logger.info("测序深度分布导表成功！")
        self.logger.info("开始进行基因组覆盖度分布导表")
        wgs_base.add_sg_area_detail(mapping_id, self.option("wgs_path").prop['path'] + "/03.map_stat/coverage")
        self.logger.info("基因组覆盖度分布导表成功！")

    def import_snp_results(self):
        snp_api = self.api.api('wgs.snp')
        snp_api.project_type = "dna_gmap"
        self.logger.info("开始进行snp统计导表")
        call_id = snp_api.add_sg_snp_call(self._sheet.project_sn, self._sheet.id, params=self.software)
        snp_api.add_sg_snp_call_stat(call_id, self.option("wgs_path").prop['path'] +
                                     "/04.snp_indel/variant_stat/snp.stat")
        self.logger.info("snp统计导表成功！")
        self.logger.info("开始进行snp质量评估导表")
        snp_api.add_snp_qc_curve(self._sheet.id, call_id,
                                 self.option("wgs_path").prop['path'] + "/04.snp_indel/variant_stat/snp.GQ", "snp_qc",
                                 "snp_qc")
        snp_api.add_snp_qc_curve(self._sheet.id, call_id,
                                 self.option("wgs_path").prop['path'] + "/04.snp_indel/variant_stat/snp.depth",
                                 "snp_depth", "snp_depth")
        self.logger.info("snp质量评估导表成功！")
        self.logger.info("开始进行snp功能注释导表")
        anno_id = snp_api.add_sg_snp_anno(self._sheet.project_sn, self._sheet.id, params='{\"method\":\"SNPEff\"}')
        # snp功能统计表
        snp_api.add_sg_snp_anno_stat(anno_id, self.option("wgs_path").prop['path'] +
                                     "/04.snp_indel/anno_stat/snp.stat", "annotion")
        # snp功效统计表
        snp_api.add_sg_snp_anno_stat(anno_id, self.option("wgs_path").prop['path'] +
                                     "/04.snp_indel/anno_stat/snp.stat", "effect")
        # snp功效与功能累加图与直方图
        snp_api.add_sg_snp_anno_bar(self._sheet.project_sn, self._sheet.id, anno_id,
                                    self.option("wgs_path").prop['path'] + "/04.snp_indel/anno_stat/snp.stat")
        self.logger.info("snp功能注释导表成功！")

    def import_indel_results(self):
        indel_api = self.api.api('wgs.indel')
        indel_api.project_type = "dna_gmap"
        self.logger.info("开始进行indel统计导表")
        call_id = indel_api.add_sg_indel_call(self._sheet.project_sn, self._sheet.id, params=self.software)
        indel_api.add_sg_indel_call_stat(call_id, self.option("wgs_path").prop['path'] +
                                         "/04.snp_indel/variant_stat/indel.stat")
        self.logger.info("indel统计表导表成功")
        self.logger.info("开始进行indel长度分布导表")
        indel_api.add_indel_length_bar(self._sheet.id, call_id,
                                       self.option("wgs_path").prop['path'] + "/04.snp_indel/variant_stat/indel.len")
        self.logger.info("indel长度分布导表成功！")
        self.logger.info("开始进行indel质量评估导表")
        indel_api.add_indel_qc_curve(self._sheet.id, call_id,
                                     self.option("wgs_path").prop['path'] + "/04.snp_indel/variant_stat/indel.GQ",
                                     "indel_qc")
        indel_api.add_indel_qc_curve(self._sheet.id, call_id,
                                     self.option("wgs_path").prop['path'] + "/04.snp_indel/variant_stat/indel.depth",
                                     "indel_depth")
        self.logger.info("indel质量评估导表成功！")
        self.logger.info("开始进行indel功能注释导表")
        anno_id = indel_api.add_sg_indel_anno(self._sheet.project_sn, self._sheet.id, params='{\"method\":\"SNPEff\"}')
        # indel功能统计表
        indel_api.add_sg_indel_anno_stat(anno_id, self.option("wgs_path").prop['path'] +
                                         "/04.snp_indel/anno_stat/indel.stat", "annotion")
        # indel功效统计表
        indel_api.add_sg_indel_anno_stat(anno_id, self.option("wgs_path").prop['path'] +
                                         "/04.snp_indel/anno_stat/indel.stat", "effect")
        # indel功效与功能累加图与直方图
        indel_api.add_sg_indel_anno_bar(self._sheet.project_sn, self._sheet.id, anno_id,
                                        self.option("wgs_path").prop['path'] + "/04.snp_indel/anno_stat/indel.stat")
        self.logger.info("indel功能注释导表成功！")

    def import_cnv_results(self):
        cnv_api = self.api.api('wgs.cnv')
        cnv_api.project_type = "dna_gmap"
        self.logger.info("开始进行cnv统计导表")
        call_id = cnv_api.add_sg_cnv_call(self._sheet.project_sn, self._sheet.id, params='{\"method\":\"CNVnator\"}')
        cnv_api.add_sg_cnv_call_stat(call_id, self.option("wgs_path").prop['path'] + "/06.cnv/cnv.stat.xls")
        self.logger.info("cnv统计表导表成功！")
        self.logger.info("开始进行cnv长度统计导表")
        cnv_api.add_cnv_length_bar(self._sheet.id, call_id,
                                   self.option("wgs_path").prop['path'] + "/06.cnv/length")   # 导入length文件夹
        self.logger.info("cnv长度统计导表成功！")

    def import_sv_results(self):
        sv_api = self.api.api('wgs.sv')
        sv_api.project_type = "dna_gmap"
        self.logger.info("开始进行sv统计导表")
        call_id = sv_api.add_sg_sv_call(self._sheet.project_sn, self._sheet.id, params='{\"method\":\"BreakDancer\"}')
        sv_api.add_sg_sv_call_stat(call_id, self.option("wgs_path").prop['path'] + "/07.sv/sv.stat.xls")
        self.logger.info("sv统计表导表成功！")
        self.logger.info("开始进行sv长度统计导表")
        sv_api.add_sv_length_bar(self._sheet.id, call_id, self.option("wgs_path").prop['path'] + "/07.sv/length")
        # 导入length文件夹
        self.logger.info("sv长度统计导表成功！")

    def import_annovar_results(self):
        anno_api = self.api.api('wgs.anno')
        anno_api.project_type = "dna_gmap"
        self.logger.info("开始进行基因功能注释导表")
        anno_id = anno_api.add_sg_anno(self._sheet.id, self._sheet.project_sn)
        anno_api.add_sg_anno_detail(self.option("wgs_path").prop['path'] + "/05.annovar/anno_count/pop.summary",
                                    anno_id, self._sheet.id, self.genome_version_id)
        self.logger.info("基因功能注释导表成功！")
        self.logger.info("开始进行GO/KEGG/EGGNOG统计图表导表")
        anno_api.sg_anno_go_stat(anno_id, self.option("wgs_path").prop['path'] +
                                 "/05.annovar/go_anno/pop.go.final.stat.detail")
        anno_api.sg_anno_kegg_stat(anno_id, self.option("wgs_path").prop['path'] +
                                   "/05.annovar/kegg_anno/pop.kegg.final.stat.detail",
                                   self.option("wgs_path").prop['path'] + "/05.annovar/kegg_anno/pathway_dir",
                                   self.wgs_path + "/05.annovar/kegg_anno/pathway_dir")
        anno_api.sg_anno_eggnog_stat(anno_id, self.option("wgs_path").prop['path'] +
                                     "/05.annovar/eggnog_anno/pop.eggnog.final.stat.detail")

    def import_ssr_results(self):
        ssr_api = self.api.api("wgs.ssr_analysis")
        ssr_api.project_type = "dna_gmap"
        self.logger.info("开始进行参考基因组ssr统计导表")
        ssr_id = ssr_api.add_sg_ssr(self._sheet.project_sn, self._sheet.id)
        ssr_api.add_sg_ssr_stat(ssr_id, self.option("wgs_path").prop['path'] + "/08.ssr/ssr.stat")
        ssr_api.add_sg_ssr_detail(ssr_id, self.option("wgs_path").prop['path'] + "/08.ssr/ssr.ref.result.xls")
        self.logger.info("参考基因组ssr统计导表成功！")

    def import_genomic_view(self):
        self.logger.info("开始进行基因组浏览器结果导表")
        self.api.del_api("wgs.wgs_base")  # 删除api引用
        wgs_base = self.api.api("wgs.wgs_base")
        wgs_base.project_type = "dna_gmap"
        wgs_base.add_sg_igv(self.wgs_path + "/03.map_stat/map_bam/sort_bams",
                            self.wgs_path + "/05.annovar/combine_variants/pop.final.vcf",
                            self.wgs_path + "/02.reference/ref.2bit")
        self.logger.info("基因组浏览器结果导表成功！")
        if self.is_cnv and self.is_sv:
            self.logger.info("开始进行circos结果导表")
            wgs_base.add_sg_circos(self.wgs_path + "/09.circos",
                                   self.option("wgs_path").prop['path'] + "/09.circos", self.total_chrlist,
                                   params='{\"chrs\": \"\"}')
            self.logger.info("circos结果导表成功！")

    def import_files_paths(self):
        file_paths = {
            # "clean_fastq_path": self.wgs_path + "/01.fastq_qc/clean_data/",
            "bam_path": self.wgs_path + "/03.map_stat/map_bam/sort_bams/",
            "indel_anno_vcf": self.wgs_path + "/04.snp_indel/eff/indel.anno.primary.vcf",
            "snp_anno_vcf": self.wgs_path + "/04.snp_indel/eff/snp.anno.primary.vcf",
            "pop_final_vcf": self.wgs_path + "/05.annovar/combine_variants/pop.final.vcf",
            "cnv_anno_path": self.wgs_path + "/06.cnv/anno/",
            "sv_anno_path": self.wgs_path + "/07.sv/anno/",
            "genome_version_id": ObjectId(self.genome_version_id),
            "gtree_hash": self.target_dir + "/temp/tree_data/Total.gTree.hash",
            "poptype": self.option("popt"),
            "bin": "yes" if self.is_bin_marker else "no",
            "ref_chrlist": self.wgs_path + '/02.reference/ref.chrlist',
            "pid": self.option("pid"),
            "mid": self.option('mid'),
            "region": self.region,
            "detail_info": self.target_dir + "/02.marker_filter/pop.filtered.detail.info",
            "daishu": self.option("daishu"),
            "pop_summary": self.wgs_path + "/05.annovar/anno_count/pop.summary"
        }
        self.api_base.update_db_record("sg_task", {"task_id": self._sheet.id}, file_paths)
        pass

    def add_all_marker_filter(self):
        """
        marker filter运行两次，导表两次
        :return:
        """
        gmap_base = self.api.api("dna_gmap.gmap_base")
        chrs = gmap_base.get_chrs(self.ref_chrlist)
        ana_dict = [{"type": "snp,indel", "pdep": "0", "odep": "0", "miss_tatio": 0, "signif": 0,
                     "gtype_stat": self.output_dir + "/temp/marker_filter/pop.filtered.gtype.stat",
                     "distribution": self.output_dir + "/temp/marker_filter/pop.filtered.marker",
                     "heatmap": self.output_dir + "/temp/marker_filter/pop.filtered.matrix",
                     "marker_stat": self.output_dir + "/temp/marker_filter/marker_info.xls",
                     "child_stat": self.output_dir + "/temp/marker_filter/pop.filtered.indi.stat",
                     "target_marker": self.target_dir + "/temp/marker_filter/pop.filtered.marker",
                     "target_info": self.target_dir + "/temp/marker_filter/pop.filtered.detail.info",
                     "data_type": "1"},
                    {"type": "snp,indel" if self.option("marker_type").lower() == "snp,indel" else self.option("marker_type"),
                     "pdep": str(self.option("pdep")), "odep": str(self.option("odep")),
                     "miss_tatio": self.option("miss_tatio"), "signif": self.option("signif"),
                     "gtype_stat": self.output_dir + "/02.marker_filter/pop.filtered.gtype.stat",
                     "distribution": self.output_dir + "/02.marker_filter/pop.filtered.marker",
                     "heatmap": self.output_dir + "/02.marker_filter/pop.filtered.matrix",
                     "marker_stat": self.output_dir + "/02.marker_filter/marker_info.xls",
                     "child_stat": self.output_dir + "/02.marker_filter/pop.filtered.indi.stat",
                     "target_marker": self.target_dir + "/02.marker_filter/pop.filtered.marker",
                     "target_info": self.target_dir + "/02.marker_filter/pop.filtered.detail.info", "data_type": "2"}]
        for m in ana_dict:
            self.import_marker_filter(m['pdep'], m['odep'], m['type'], m['miss_tatio'], m['signif'], m['gtype_stat'],
                                      m['distribution'], m['heatmap'], m['marker_stat'], m['child_stat'], chrs,
                                      m['target_marker'], m['target_info'], m['data_type'])

    def import_marker_filter(self, pdep, odep, type_, miss_tatio, signif, gtype_stat, distribution, heatmap,
                             marker_stat, child_stat, chr_list, target_marker, target_info, data_type):
        """
        一次运行结果导入marker过滤的结果表
        "type": "ALL", "pdep": "0", "odep": "0", "miss_tatio": 0, "signif": 0
        :return:
        """
        self.logger.info("开始导入marker filter数据")
        marker_filter = self.api.api("dna_gmap.marker_filter")
        if data_type == '1':
            marker_filter.add_sg_subtype_matrix(self._sheet.id, target_info)
        else:
            params = {"matrix_id": "", "pid": self.option("pid"), "mid": self.option("mid"),
                      "type": type_, "pdep": pdep, "odep": odep, "popt": self.option("popt"),
                      "miss_tatio": miss_tatio, "signif": signif, "submit_location": "marker", "task_type": 2}
            self.marker_id = marker_filter.add_sg_marker(self._sheet.project_sn, self._sheet.id, type_, params,
                                                         detail_info_path=target_info, filter_marker_path=target_marker,
                                                         chr_list=chr_list)
            marker_filter.add_sg_marker_detail(self._sheet.id, self.marker_id, gtype_stat, type_)
            marker_filter.add_sg_distribution(self._sheet.id, self.marker_id, distribution)
            # marker_filter.add_sg_heatmap(self._sheet.id, self.marker_id, heatmap)
            marker_filter.add_sg_marker_stat(self._sheet.id, self.marker_id, marker_stat)
            marker_filter.add_sg_marker_child_stat(self._sheet.id, self.marker_id, child_stat)
        self.logger.info("导入marker filter数据成功")

    def import_bin_marker(self):
        """
        导入bin marker的结果文件
        :return:
        """
        self.logger.info("开始导入bin_marker数据")
        bin_marker = self.api.api("dna_gmap.binner_calculate")
        params = {"matrix_id": str(self.marker_id), "window_size": 0.5, "window_step": 100,
                  "submit_location": "binmarker", "task_type": 2, "upload_marker": False}
        self.marker_id = bin_marker.add_sg_binmarker(self._sheet.project_sn, self._sheet.id, params)
        bin_marker.add_sg_binmarker_bin(self.marker_id, self.output_dir + "/03.binner/bin_stat.xls")
        bin_marker.add_sg_binmarker_var(self.marker_id, self.output_dir + "/03.binner/bin_info.xls")
        bin_marker.add_sg_binmarker_var_detail(self.marker_id, self.output_dir + "/03.binner/bin_pos.xls")
        self.api_base.update_db_record("sg_binmarker", {"_id": self.marker_id},
                                       {"filtered_marker_path": self.target_dir + "/03.binner/Total.bin.marker",
                                        "detail_info_path":
                                            self.target_dir + "/02.marker_filter/pop.filtered.detail.info"})
        self.logger.info("导入bin_marker数据成功")

    def import_linkage_grouping(self):
        """
        导入连锁分群的结果表格
        :return:
        """
        self.logger.info("开始导入连锁分群结果数据")
        if self.is_bin_marker:
            marker_ = self.output_dir + "/03.binner/Total.bin.marker"
        else:
            marker_ = self.output_dir + "/02.marker_filter/pop.filtered.marker"
        linkage = self.api.api("dna_gmap.linkage_groupping")
        params = {"marker_id": str(self.marker_id), "group_type": self.option("group_type"),
                  "submit_location": "lg", "task_type": 2}
        if self.option("group_type") == "mlod":
            params.update({"chr_num": self.option("chr_num"), "start_lod": self.option("lod_start"),
                           "end_lod": self.option("lod_end"), "min_group": self.option("min_group"),
                           "max_group": self.option("max_group")})
        lg_id = linkage.add_sg_lg(self._sheet.id, self._sheet.project_sn, params=params)
        evalutaion_id = linkage.add_sg_evalutaion(lg_id)
        if self.option("popt").lower() not in ["f1", 'cp']:
            linkage.add_sg_evalutaion_stat(self.output_dir + "/04.grouping/evalutaion/total.mapstat", evalutaion_id)
            linkage.add_yichuantupu(self._sheet.id, evalutaion_id,
                                    self.output_dir + "/04.grouping/evalutaion/total.map", "total", True)
            linkage.collinearity(self._sheet.id, self.output_dir + "/04.grouping/evalutaion/total.map", evalutaion_id,
                                 self.output_dir + "/04.grouping/evalutaion/total.phy.spearman.xls", "total")
            linkage.add_sg_evalutaion_detail(self.output_dir + "/04.grouping/evalutaion/total.marker.info",
                                             evalutaion_id)
            marker_info_path = self.output_dir + "/04.grouping/evalutaion/total.marker.info"
            linkage.update_sg_lg(lg_id, self.target_dir + "/04.grouping/evalutaion/total.loc",
                                 self.target_dir + "/04.grouping/evalutaion/total.map", "nocp",
                                 self.target_dir + "/04.grouping/evalutaion/total.csv")
            linkage.add_sg_lg_detail(self.output_dir + "/04.grouping/groupping/Total.lg",
                                     self.output_dir + "/04.grouping/evalutaion/total.map", marker_, lg_id)
            linkage.add_sg_reorganization_heatmap(self.output_dir + "/04.grouping/evalutaion/fig",
                                                  self.target_dir + "/04.grouping/evalutaion/fig", evalutaion_id)
            # linkage.add_genetype_heatmap(self.output_dir + "/04.grouping/evalutaion/total.csv", "total",
            #                              self._sheet.id, evalutaion_id)
            self.api_base.update_db_record("sg_evalutaion", {"_id": evalutaion_id}, {"parent_source": ['total']})
        else:
            marker_info_path = self.output_dir + "/04.grouping/evalutaion/total.sexAver.info"
            linkage.add_sg_evalutaion_stat(self.output_dir + "/04.grouping/evalutaion/male.mapstat", evalutaion_id,
                                           "male")
            linkage.add_sg_evalutaion_stat(self.output_dir + "/04.grouping/evalutaion/sexAver.mapstat", evalutaion_id,
                                           "sexaver")
            linkage.add_sg_evalutaion_stat(self.output_dir + "/04.grouping/evalutaion/female.mapstat", evalutaion_id,
                                           "female")
            linkage.add_yichuantupu(self._sheet.id, evalutaion_id,
                                    self.output_dir + "/04.grouping/evalutaion/total.female.map", "female", True)
            linkage.add_yichuantupu(self._sheet.id, evalutaion_id,
                                    self.output_dir + "/04.grouping/evalutaion/total.male.map", "male", True)
            linkage.add_yichuantupu(self._sheet.id, evalutaion_id,
                                    self.output_dir + "/04.grouping/evalutaion/total.sexAver.map", "sexaver", True)
            linkage.collinearity(self._sheet.id, self.output_dir + "/04.grouping/evalutaion/total.female.map",
                                 evalutaion_id,
                                 self.output_dir + "/04.grouping/evalutaion/total.female.phy.spearman.xls", "female")
            linkage.collinearity(self._sheet.id, self.output_dir + "/04.grouping/evalutaion/total.male.map",
                                 evalutaion_id, self.output_dir + "/04.grouping/evalutaion/total.male.phy.spearman.xls",
                                 "male")
            linkage.collinearity(self._sheet.id, self.output_dir + "/04.grouping/evalutaion/total.sexAver.map",
                                 evalutaion_id,
                                 self.output_dir + "/04.grouping/evalutaion/total.sexaver.phy.spearman.xls", "sexaver")
            linkage.add_sg_evalutaion_detail(self.output_dir + "/04.grouping/evalutaion/total.male.info",
                                             evalutaion_id, "male")
            linkage.add_sg_evalutaion_detail(self.output_dir + "/04.grouping/evalutaion/total.female.info",
                                             evalutaion_id, "female")
            linkage.add_sg_evalutaion_detail(self.output_dir + "/04.grouping/evalutaion/total.sexAver.info",
                                             evalutaion_id, "sexaver")
            linkage.update_sg_lg(lg_id, self.target_dir + "/04.grouping/total.sexAver.loc",
                                 self.target_dir + "/04.grouping/total.sexAver.map", "cp")
            linkage.add_sg_lg_detail(self.output_dir + "/04.grouping/groupping/Total.lg",
                                     self.output_dir + "/04.grouping/evalutaion/total.sexAver.map", marker_, lg_id)
            linkage.add_sg_reorganization_heatmap(self.output_dir + "/04.grouping/evalutaion/fig",
                                                  self.target_dir + "/04.grouping/evalutaion/fig", evalutaion_id, "cp")
            # linkage.add_genetype_heatmap(self.output_dir + "/04.grouping/evalutaion/total.female.phase", "female",
            #                              self._sheet.id, evalutaion_id)
            # linkage.add_genetype_heatmap(self.output_dir + "/04.grouping/evalutaion/total.male.phase", "male",
            #                              self._sheet.id, evalutaion_id)
            # linkage.add_genetype_heatmap(self.output_dir + "/04.grouping/evalutaion/total.sexAver.phase", "sexaver",
            #                              self._sheet.id, evalutaion_id)
            self.api_base.update_db_record("sg_evalutaion", {"_id": evalutaion_id}, {"parent_source": ["female", "male",
                                                                                                       "sexaver"]})
        linkage.add_sg_lg_stat(self.output_dir + "/04.grouping/marker_stat/total.marker.stat.xls", lg_id)
        linkage.add_sg_tree(self.output_dir + "/temp/tree_data/Total.gTree.hash.dumper",
                            self._sheet.id, "grouping_tree", marker_info_path)
        self.api_base.update_db_record("sg_task", {"task_id": self._sheet.id}, {"origin_lg": lg_id})
        self.api_base.update_db_record("sg_lg", {"_id": lg_id},
                                       {"total_lg": self.target_dir + "/04.grouping/groupping/Total.lg",
                                        "marker_info_path": self.target_dir + "/04.grouping/evalutaion/",
                                        "marker_id": self.marker_id,
                                        "marker_type": "binmarker" if self.is_bin_marker else "marker"})
        self.logger.info("导入连锁分群结果数据成功")

    def end(self):
        self.send_files()
        self.logger.info("上传文件目录完成！")
        self.run_api()
        self.logger.info("全部导表完成！")
        super(GmapWorkflow, self).end()

    def get_target_dir(self):
        """
        获取远程磁盘的路径
        :return:
        """
        self.target_dir = self._sheet.output.rstrip("/")
        self.region = self._sheet.output.strip().split('://')[0]
        path = self.get_wgs_path()
        # m = re.match(r'(.*):(.*)', path)
        if not path.startswith("/mnt"):
            self.pop_final_vcf = os.path.join(self.option('wgs_path').prop['path'].rstrip('/') +
                                              "/05.annovar/combine_variants/pop.final.vcf")
            self.ref_chrlist = os.path.join(self.option('wgs_path').prop['path'].rstrip('/') +
                                            "/02.reference/ref.chrlist")
            self.wgs_path = path.rstrip('/')
        else:
            self.wgs_path = path.rstrip('/')
            self.pop_final_vcf = path.rstrip('/') + "/05.annovar/combine_variants/pop.final.vcf"
            self.ref_chrlist = path.rstrip('/') + "/02.reference/ref.chrlist"
        self.logger.info("wgs_path:{}".format(self.wgs_path))
        self.logger.info("pop_final_vcf:{}".format(self.pop_final_vcf))

    def check_chr_ref(self):
        """
        检查基因组 是不是都是sca水平，如果sca水平的话 就不能使用参考染色体去分群
        :return:
        """
        wgs_base = self.api.api("wgs.wgs_base")
        self.total_chrlist = wgs_base.total_chrlist(self.json_path, self.genome_version_id, self.option("group_type"))

    def get_file_len(self, file_path):
        """
        获取文件行数
        :return:
        """
        return len(open(file_path, 'rU').readlines())

    def get_wgs_path(self):
        """
        从mapping_file.txt中获取wgs的路径，用于后面导表使用
        :return:
        """
        file_path = self.work_dir + "/remote_input/wgs_path/mapping_file.txt"
        if not os.path.exists(file_path):
            raise Exception("file:{} is not exists!")
        with open(file_path, 'r') as r:
            data = r.readlines()[0]
            json_data = json.loads(data)
        temp_path = json_data['wgs_path'][0]["file_path"].split("workflow_results")[0].rstrip('/') + '/workflow_results'
        return temp_path


