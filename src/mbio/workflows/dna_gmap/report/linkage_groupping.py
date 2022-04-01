# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180701

import re
import os
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class LinkageGrouppingWorkflow(Workflow):
    """
    标记连锁分群-的接口workflow，里面主要就是linkage_groupping这个module，然后就是导表了
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LinkageGrouppingWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "marker", "type": "infile", "format": "dna_gmap.marker"},
            # 输入marker文件， pop.filtered.marker或者Total.bin.marker
            {"name": "is_ref", "type": "string", "default": "true"},  # 判断是否是有参
            {"name": "group_type", "type": "string", "default": 'ref'},  # mlod or ref
            {"name": "key_file", "type": "string", "defult": "Total"},  # 输入Key值
            {"name": "scaf_marker", "type": "string", "default": "no"},  # yes就添加sca 否则不添加
            {"name": "chr_num", "type": "int"},
            {"name": "poptype", "type": "string"},  # 传入poptype参数，例：CP,F2等
            {"name": "bin", "type": "string", "default": "no"},  # 用于判断是不是进行了bin计算
            {"name": "ref_chrlist", "type": "infile", 'format': 'bsa.vcf'},
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
            {"name": "marker_info_path", "type": "infile", 'format': 'bsa.vcf'},
            {"name": "useless_marker", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "analysis_type", "type": "string", "default": "1"},  # 用于区分生成连锁分群的方式，有1,2,3三种方式
            # 1是根据程序去生成连锁分群，当为2为手动去进行连锁分，3是图谱评估结果中生成新的连锁分群
            {"name": "pop_marker_detail", "type": "infile", "format": "dna_gmap.marker"},  # pop.filterd.detail.info
            {"name": "marker_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.linkage_groupping = self.add_module('dna_gmap.linkage_groupping')
        self.target_dir = ''

    def check_options(self):
        if not self.option('marker'):
            raise OptionError('必须提供marker结果表', code="14800301")
        if self.option('analysis_type'):
            if self.option('analysis_type') not in ["1", "2", "3"]:
                raise OptionError("分析方法类型不合法！必须为1 or 2 or 3", code="14800302")
        else:
            raise OptionError('必须提供analysis_type结果表', code="14800303")
        if not self.option("pop_marker_detail"):
            raise OptionError("必须要提供pop_marker_detail", code="14800304")
        return True

    def linkage_groupping_run(self):
        opt = {
            'marker': self.option('marker').prop['path'],
            "pop_marker_detail": self.option("pop_marker_detail").prop['path'],
            'ref_chrlist': self.option("ref_chrlist").prop['path'],
            'poptype': self.option("poptype"),
            'bin': self.option("bin"),
            'is_ref': self.option("is_ref"),
            "analysis_type": self.option("analysis_type")
        }
        if self.option("analysis_type") == "1":
            opt.update({
                'group_type': self.option('group_type')
            })
            if self.option('group_type') == "mlod":
                opt.update({
                    'scaf_marker': self.option('scaf_marker'),
                    'chr_num': self.option('chr_num'),
                    'start_lod': self.option('start_lod'),
                    'end_lod': self.option('end_lod'),
                    'stepsize_lod': self.option('stepsize_lod'),
                    'min_group': self.option('min_group'),
                    'max_group': self.option('max_group')
                })
        elif self.option("analysis_type") == "2":
            opt.update({
                "gtree_hash": self.option("gtree_hash").prop['path'],
                "tree_nodes": self.option("tree_nodes")
            })
        elif self.option("analysis_type") == "3":
            opt.update({
                'total_lg': self.option('total_lg').prop["path"],
                'miss_ratio_start': self.option('miss_ratio_start'),
                'miss_ratio_end': self.option('miss_ratio_end'),
                'signif_start': self.option('signif_start'),
                'signif_end': self.option('signif_end'),
                'marker_info_path': self.option('marker_info_path').prop['path'],
                'useless_marker': self.option('useless_marker')
            })
        self.linkage_groupping.set_options(opt)
        self.linkage_groupping.on('end', self.set_output, 'linkage_groupping')
        self.linkage_groupping.run()

    def check_markers(self):
        file_path = self.option('marker').prop['path']
        lens = len(open(file_path, 'rU').readlines())
        if lens < 22:
            return 4, lens
        else:
            return 20, lens

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'linkage_groupping':
            self.linkdir(obj.output_dir, self.output_dir)
        self.set_db()
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
        if os.path.exists(self.output_dir + "/groupping/Total.lg") or os.path.exists(
                        self.output_dir + "/hand_groupping/Total.lg") or os.path.exists(
                    self.output_dir + "/get_markers/total.lg"):
            self.logger.info("设置linkage_groupping的导表！")
            api = self.api.api("dna_gmap.linkage_groupping")
            base_api = self.api.api("dna_gmap.api_base")
            # 遗传图谱评估
            evalutaion_id = api.add_sg_evalutaion(self.option("main_id"))
            if self.option("poptype").lower() in ['cp', "f1"]:
                api.add_sg_evalutaion_stat(self.output_dir + "/evalutaion/male.mapstat", evalutaion_id, "male")
                api.add_sg_evalutaion_stat(self.output_dir + "/evalutaion/sexAver.mapstat", evalutaion_id, "sexaver")
                api.add_sg_evalutaion_stat(self.output_dir + "/evalutaion/female.mapstat", evalutaion_id, "female")
                api.add_yichuantupu(self.option("task_id"), evalutaion_id,
                                    self.output_dir + "/evalutaion/total.female.map", "female", True)
                api.add_yichuantupu(self.option("task_id"), evalutaion_id,
                                    self.output_dir + "/evalutaion/total.male.map", "male", True)
                api.add_yichuantupu(self.option("task_id"), evalutaion_id,
                                    self.output_dir + "/evalutaion/total.sexAver.map", "sexaver", True)
                api.collinearity(self.option("task_id"), self.output_dir + "/evalutaion/total.female.map",
                                 evalutaion_id,
                                 self.output_dir + "/evalutaion/total.female.phy.spearman.xls", "female")
                api.collinearity(self.option("task_id"), self.output_dir + "/evalutaion/total.male.map", evalutaion_id,
                                 self.output_dir + "/evalutaion/total.male.phy.spearman.xls", "male")
                api.collinearity(self.option("task_id"), self.output_dir + "/evalutaion/total.sexAver.map",
                                 evalutaion_id,
                                 self.output_dir + "/evalutaion/total.sexaver.phy.spearman.xls", "sexaver")
                api.add_sg_evalutaion_detail(self.output_dir + "/evalutaion/total.male.info", evalutaion_id, "male")
                api.add_sg_evalutaion_detail(self.output_dir + "/evalutaion/total.female.info", evalutaion_id, "female")
                api.add_sg_evalutaion_detail(self.output_dir + "/evalutaion/total.sexAver.info", evalutaion_id,
                                             "sexaver")
                api.update_sg_lg(self.option("main_id"), self.target_dir + "/evalutaion/total.sexAver.loc",
                                 self.target_dir + "/evalutaion/total.sexAver.map", "cp")
                api.add_sg_lg_detail(self.output_dir + "/groupping/Total.lg",
                                     self.output_dir + "/evalutaion/total.sexAver.map",
                                     self.option("marker").prop['path'], ObjectId(self.option("main_id")))
                api.add_sg_reorganization_heatmap(self.output_dir + "/evalutaion/fig",
                                                  self.target_dir + "/evalutaion/fig", evalutaion_id, "cp")

                # api.add_genetype_heatmap(self.output_dir + "/evalutaion/total.female.phase", "female",
                #                          self.option("task_id"), evalutaion_id)
                # api.add_genetype_heatmap(self.output_dir + "/evalutaion/total.male.phase", "male",
                #                          self.option("task_id"), evalutaion_id)
                # api.add_genetype_heatmap(self.output_dir + "/evalutaion/total.sexAver.phase", "sexaver",
                #                          self.option("task_id"), evalutaion_id)
                base_api.update_db_record("sg_evalutaion", {"_id": evalutaion_id}, {"parent_source": ["female", "male",
                                                                                                      "sexaver"]})
            else:
                api.add_sg_evalutaion_stat(self.output_dir + "/evalutaion/total.mapstat", evalutaion_id)
                api.add_yichuantupu(self.option("task_id"), evalutaion_id,
                                    self.output_dir + "/evalutaion/total.map", "total", True)
                api.collinearity(self.option("task_id"), self.output_dir + "/evalutaion/total.map", evalutaion_id,
                                 self.output_dir + "/evalutaion/total.phy.spearman.xls", "total")
                api.add_sg_evalutaion_detail(self.output_dir + "/evalutaion/total.marker.info", evalutaion_id)
                api.update_sg_lg(self.option("main_id"), self.target_dir + "/evalutaion/total.loc",
                                 self.target_dir + "/evalutaion/total.map", "nocp",
                                 self.target_dir + "/evalutaion/total.csv")
                api.add_sg_lg_detail(self.output_dir + "/groupping/Total.lg",
                                     self.output_dir + "/evalutaion/total.map", self.option("marker").prop['path'],
                                     ObjectId(self.option("main_id")))
                api.add_sg_reorganization_heatmap(self.output_dir + "/evalutaion/fig",
                                                  self.target_dir + "/evalutaion/fig", evalutaion_id)
                base_api.update_db_record("sg_evalutaion", {"_id": evalutaion_id}, {"parent_source": ['total']})
            api.add_sg_lg_stat(self.output_dir + "/marker_stat/total.marker.stat.xls", ObjectId(self.option("main_id")))
            base_api.update_db_record("sg_lg", {"_id": ObjectId(self.option("main_id"))},
                                      {"total_lg": self.target_dir + "/groupping/Total.lg",
                                       "marker_id": ObjectId(self.option('marker_id')),
                                       "marker_type": "binmaker" if self.option("bin") == "yes" else "marker",
                                       "marker_info_path": self.target_dir + "/evalutaion/"})
            self.logger.info("设置linkage_groupping的导表成功！")
        else:
            self.logger.info("分群没有成功不进行导表！")
        self.end()

    def run(self):
        self.get_target_dir()
        self.linkage_groupping_run()
        super(LinkageGrouppingWorkflow, self).run()

    def end(self):
        """
        这里后面要重新定义下文件名字
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(LinkageGrouppingWorkflow, self).end()

    def get_target_dir(self):
        """
        获取远程磁盘的路径
        :return:
        """
        # self.target_dir = self._sheet.output.strip().split('://')[1]
        self.target_dir = self._sheet.output.rstrip('/')
