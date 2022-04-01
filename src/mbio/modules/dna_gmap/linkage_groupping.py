#!/usr/bin/env python
# -*- coding: utf-8 -*-

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os


class LinkageGrouppingModule(Module):
    """
    标记连锁分群-的module，包含：Module：mlod_calc， Tool：groupping， Moudle：split_markers， Module：marker_order，
     Tool：evalutaion, Tool:hand_groupping, Tool:get_markers
    version 1.0
    author: HONGDONG
    last_modify: 20180619
    """
    def __init__(self, work_id):
        super(LinkageGrouppingModule, self).__init__(work_id)
        options = [
            {"name": "marker", "type": "infile", "format": "dna_gmap.marker"},
            # 输入marker文件， pop.filtered.marker或者Total.bin.marker
            {"name": "is_ref", "type": "string", "default": "true"},  # 判断是否是有参
            {"name": "group_type", "type": "string", "default": 'ref'},  # mlod or ref
            {"name": "key_file", "type": "string", "defult": "Total"},  # 输入Key值
            {"name": "scaf_marker", "type": "string", "default": "no"},  # yes就添加sca 否则不添加
            {"name": "chr_num", "type": "int"},
            {"name": "poptype", "type": "string"},  # 传入poptype参数，例：F1,F2等
            {"name": "bin", "type": "string", "default": "no"},  # 用于判断是不是进行了bin计算
            {"name": "ref_chrlist", "type": "string"},
            {"name": "start_lod", "type": "int", "default": 3},
            {"name": "end_lod", "type": "int", "default": 20},
            {"name": "stepsize_lod", "type": "int", "default": 1},
            {"name": "min_group", "type": "int", "default": 20},
            {"name": "max_group", "type": "int", "default": 500},
            {"name": "gtree_hash", "type": "infile", "format": "dna_gmap.gtree_hash"},
            {"name": "tree_nodes", "type": "string"},
            {"name": "total_lg", "type": "infile", "format": "dna_gmap.lg"},
            {"name": "miss_ratio_start", "type": "float", "default": 30},
            {"name": "miss_ratio_end", "type": "float", "default": 100},
            {"name": "signif_start", "type": "float", "default": 0.05},
            {"name": "signif_end", "type": "float", "default": 1},
            {"name": "marker_info_path", "type": "string"},
            {"name": "useless_marker", "type": "string"},
            {"name": "analysis_type", "type": "string", "default": "1"}, # 用于区分生成连锁分群的方式，有1,2,3三种方式
            # 1是根据程序去生成连锁分群，当为2为手动去进行连锁分，3是图谱评估结果中生成新的连锁分群
            {"name": "total_mlod", "type": "outfile", "format": "dna_gmap.mlod"},
            {"name": "pop_marker_detail", "type": "infile", "format": "dna_gmap.marker"}  # pop.filterd.detail.info
        ]
        self.add_option(options)
        self.mlod_calc = self.add_module('dna_gmap.mlod_calc')
        self.groupping = self.add_tool('dna_gmap.groupping')
        self.hand_groupping = self.add_tool('dna_gmap.hand_groupping')
        self.split_markers = self.add_module('dna_gmap.split_markers')
        self.marker_order = self.add_module('dna_gmap.marker_order')
        self.evalutaion = self.add_tool('dna_gmap.map_evaluation')
        self.get_markers = self.add_tool('dna_gmap.get_markers')
        self.group_result = self.add_tool('dna_gmap.get_grouping_result')
        self.total_lg_file = ""

    def check_options(self):
        """
        检查参数
        """
        if not self.option('marker'):
            raise OptionError('必须提供marker结果表', code="24800101")
        if self.option('analysis_type'):
            if self.option('analysis_type') not in ["1", "2", "3"]:
                raise OptionError("分析方法类型不合法！必须为1 or 2 or 3", code="24800102")
        else:
            raise OptionError('必须提供analysis_type结果表', code="24800103")
        if not self.option("pop_marker_detail"):
            raise OptionError("必须要输入pop_marker_detail参数！", code="24800104")

    def mlod_calc_run(self):
        self.mlod_calc.set_options({
            'marker': self.option('marker').prop['path'],
            'only': 200
        })
        self.mlod_calc.on('end', self.set_output, 'mlod_calc')
        self.mlod_calc.run()

    def groupping_run(self):
        self.groupping.set_options({
            'group_type': self.option('group_type'),
            'total_mlod': self.mlod_calc.option("total_mlod"),
            'marker': self.option('marker'),
            'scaf_marker': self.option('scaf_marker'),
            'chr_num': self.option('chr_num'),
            'start_lod': self.option('start_lod'),
            'end_lod': self.option('end_lod'),
            'stepsize_lod': self.option('stepsize_lod'),
            'min_group': self.option('min_group'),
            'max_group': self.option('max_group')
        })
        self.groupping.on('end', self.set_output, 'groupping')
        self.groupping.run()

    def hand_groupping_run(self):
        """
        手动进行分群
        :return:
        """
        self.hand_groupping.set_options({
            'gtree_hash': self.option('gtree_hash').prop["path"],
            'tree_nodes': self.option('tree_nodes')
        })
        self.hand_groupping.on('end', self.set_output, 'hand_groupping')
        self.hand_groupping.run()

    def get_markers_run(self):
        """图谱评估中重新运行"""
        self.get_markers.set_options({
            'total_lg': self.option('total_lg').prop["path"],
            'miss_ratio_start': self.option('miss_ratio_start'),
            'miss_ratio_end': self.option('miss_ratio_end'),
            'signif_start': self.option('signif_start'),
            'signif_end': self.option('signif_end'),
            'marker_info_path': self.option('marker_info_path'),
            'useless_marker': self.option('useless_marker')
        })
        self.get_markers.on('end', self.set_output, 'get_markers')
        self.get_markers.run()

    def split_markers_run(self):
        if self.option("analysis_type") == "1" and not self.groupping.option("total_lg").is_set:
            self.logger.info("total.lg为空")
            self.end()
        else:
            opt = {
                "bin": self.option("bin"),  # 判断是否为bin，传入参数为yes/no
                "poptype": "CP" if self.option('poptype').lower() in ['f1', 'cp'] else self.option('poptype'),
                # 传入poptype参数，例：CP,F2等
                "marker_path": self.option('marker').prop['path'],  # 传入marker文件，没bin：pop.filtered.marker；
                #  有bin：Total.bin.marker
                "ref": "yes" if self.option("is_ref") == "true" else 'no'  # 判断是否为ref，传入参数为yes/no
            }
            if self.option("analysis_type") == '1':
                opt.update({"lg_path": self.groupping.option("total_lg")})
            elif self.option("analysis_type") == '2':
                opt.update({"lg_path": self.hand_groupping.option("total_lg")})
            else:
                opt.update({"lg_path": self.get_markers.option("filter_lg")})
            self.split_markers.set_options(opt)
            self.split_markers.on('end', self.set_output, 'split_markers')
            self.split_markers.on('end', self.marker_order_run)
            self.split_markers.on('end', self.group_result_run)
            self.split_markers.run()

    def marker_order_run(self):
        self.marker_order.set_options({
            'marker_list': self.split_markers.output_dir + "/ref.marker.list",
            'poptype': "CP" if self.option('poptype').lower() in ['f1', 'cp'] else self.option('poptype'),
            "ref": "yes" if self.option("is_ref") == "true" else 'no',  # 判断是否为ref，传入参数为yes/no
            "pri_marker_path": self.split_markers.option("pri_marker_path").prop['path'],
            "bin": self.option("bin")
        })
        self.marker_order.on('end', self.set_output, 'marker_order')
        self.marker_order.on('end', self.evalutaion_run)
        self.marker_order.run()

    def evalutaion_run(self):
        option = {
            'map_cycle_dir': self.marker_order.output_dir + "/map_circle3",
            'pop_type': "F1" if self.option('poptype').lower() in ['f1', 'cp'] else self.option('poptype'),
            'ref_chrlist': self.option("ref_chrlist"),
            'is_ref': self.option('is_ref')
        }
        if self.option('poptype').lower() in ['f1', 'cp']:
            option.update({"filtered_marker": self.option('marker')})
        self.evalutaion.set_options(option)
        self.evalutaion.on('end', self.set_output, 'evalutaion')
        self.evalutaion.run()

    def group_result_run(self):
        option = {
            'pop_marker_detail': self.option("pop_marker_detail").prop['path']
        }
        if self.option("analysis_type") == '1':
            option.update({"total_lg": self.groupping.option("total_lg")})
        elif self.option("analysis_type") == '2':
            option.update({"total_lg": self.hand_groupping.option("total_lg")})
        else:
            option.update({"total_lg": self.get_markers.option("filter_lg")})
        self.group_result.set_options(option)
        self.group_result.on('end', self.set_output, 'marker_stat')
        self.group_result.run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
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

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'mlod_calc':
            self.linkdir(obj.output_dir, 'mlod_calc')
        elif event['data'] == 'groupping':
            self.linkdir(obj.output_dir, 'groupping')
        elif event['data'] == 'hand_groupping':
            self.linkdir(obj.output_dir, 'hand_groupping')
        elif event['data'] == 'get_markers':
            self.linkdir(obj.output_dir, 'get_markers')
        elif event['data'] == 'split_markers':
            self.linkdir(obj.output_dir, 'split_markers')
        elif event['data'] == 'marker_order':
            self.linkdir(obj.output_dir, 'marker_order')
        elif event['data'] == 'evalutaion':
            self.linkdir(obj.output_dir, 'evalutaion')
        elif event['data'] == 'marker_stat':
            self.linkdir(obj.output_dir, 'marker_stat')
        else:
            pass

    def run(self):
        super(LinkageGrouppingModule, self).run()
        self.on_rely([self.group_result, self.evalutaion], self.end)
        if self.option("analysis_type") == '1':
            self.logger.info("执行ref_mlod分析接口")
            self.mlod_calc.on('end', self.groupping_run)
            self.groupping.on('end', self.split_markers_run)
            self.mlod_calc_run()
        elif self.option("analysis_type") == '2':
            self.logger.info("执行图谱手动分群接口")
            self.hand_groupping.on('end', self.split_markers_run)
            self.hand_groupping_run()
        elif self.option("analysis_type") == '3':
            self.logger.info("执行图谱评估结果筛选接口")
            self.get_markers.on('end', self.split_markers_run)
            self.get_markers_run()

    def end(self):
        if self.option("analysis_type") == '1':  # 只有默认的ref或者mlod的时候才会有total.mlod文件
            self.option("total_mlod").set_path(self.mlod_calc.option("total_mlod").prop['path'])
        super(LinkageGrouppingModule, self).end()
