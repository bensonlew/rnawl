# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# last_modify:20180620

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import time


class MarkerOrderModule(Module):
    """
    排图，运行三次。
    lasted modified by hongdong@20180731 添加对hkxhk标记的判断
    lasted modified by hongdong@20190103 添加对f1群体排图结果进行判断，如果有map文件遗传距离都是0，就直接用第一次排图结果
    """
    def __init__(self, work_id):
        super(MarkerOrderModule, self).__init__(work_id)
        options = [
            {"name": "marker_list", "type": "infile", "format": "dna_gmap.marker_list"},  # pri.marker.list/ref.marker.list
            {"name": "poptype", "type": "string"},  # 传入poptype参数，例：CP,F2等
            {"name": "ref", "type": "string"},  # 判断是否为ref，传入参数为yes/no
            {"name": "min_lod", "type": "string", "default": "-1"},
            {"name": "randomise_order", "type": "string", "default": "1"},
            {"name": "knn", "type": "string", "default": "25"},
            {"name": "bin", "type": "string", "default": "no"},
            {"name": "pri_marker_path", "type": "string"},
            {"name": "map_cycle_dir", "type": "outfile", "format": "dna_gmap.map_cycle_dir"}  # 输出文件夹map_circle3
        ]
        self.add_option(options)
        self.marker_order1_tools = []
        self.marker_order2_tools = []
        self.marker_order3_tools = []

    def check_options(self):
        if not self.option("marker_list"):
            raise OptionError("缺少marker_list参数", code="24800201")
        if not self.option("poptype"):
            raise OptionError("缺少poptype参数", code="24800202")
        if not self.option("ref"):
            raise OptionError("缺少ref参数", code="24800203")
        if not self.option("min_lod"):
            raise OptionError("缺少min_lod参数", code="24800204")
        if not self.option("randomise_order"):
            raise OptionError("缺少randomise_order参数", code="24800205")
        if not self.option("knn"):
            raise OptionError("缺少knn参数", code="24800206")
        return True

    def get_list(self, file_path, file_name=None):
        """
        生成下一步要用的list
        file:module output中的文件夹名称
        """
        if file_name:
            file_name_ = os.path.join(file_path, (file_name + "/marker.list"))
            file_path_ = os.path.join(file_path, file_name)
        else:
            file_name_ = file_path + "/marker.list"
            file_path_ = file_path
        if os.path.exists(file_name_):
            os.remove(file_name_)
        with open(file_name_, "w") as fg:
            self.logger.info("file_path: {}".format(file_name_))
            for root, dirs, files in os.walk(file_path_):
                for name in files:
                    if re.match(".*\.correct\.marker$", name) or re.match(".*\.correct\.loc$", name):
                        self.logger.info('{} has file：{}'.format(file_path_, name))
                        name_split = name.strip().split(".")
                        chr_name = name_split[0]
                        if file_name:
                            loc_path = os.path.join(os.path.join(file_path, file_name), name)
                        else:
                            loc_path = os.path.join(file_path, name)
                        fg.write(chr_name + "\t" + loc_path + "\n")

    def map_circle1_run(self):
        with open(self.option("marker_list").prop["path"], "r")as f1:
            lines = f1.readlines()
            for line in lines:
                data_split = line.strip().split("\t")
                correct_loc = data_split[1]
                map_circle = self.add_tool("dna_gmap.marker_order")
                map_circle.set_options({
                    "correct_loc": correct_loc,
                    "poptype": self.option("poptype"),
                    "ref": self.option("ref"),
                    "min_lod": self.option("min_lod"),
                    "randomise_order": self.option("randomise_order"),
                    "knn": self.option("knn"),
                    "bin": self.option('bin'),
                    "pri_marker_path": self.option('pri_marker_path'),
                    "is_first": "true"
                })
                self.marker_order1_tools.append(map_circle)
        for j in range(len(self.marker_order1_tools)):
            self.marker_order1_tools[j].on("end", self.set_output, "map_circle1")
        if self.marker_order1_tools:
            if len(self.marker_order1_tools) > 1:
                self.on_rely(self.marker_order1_tools, self.map_circle2_run)
            elif len(self.marker_order1_tools) == 1:
                self.marker_order1_tools[0].on('end', self.map_circle2_run)
            for tool in self.marker_order1_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("marker_order1_tools列表为空！", code="24800207")

    def map_circle2_run(self):
        self.logger.info("第一次排图结束，开始第二次排图")
        if self.check_hkxhk():
            self.end()
        else:
            if len(self.marker_order1_tools) == 1:
                self.get_list(self.marker_order1_tools[0].output_dir)
                marker_path = os.path.join(self.marker_order1_tools[0].output_dir, "marker.list")
            else:
                self.get_list(self.output_dir, "map_circle1")
                marker_path = os.path.join(self.output_dir, "map_circle1/marker.list")
            self.logger.info("marker_path: {}".format(marker_path))
            with open(marker_path, "r")as f1:
                lines = f1.readlines()
                for line in lines:
                    data_split = line.strip().split("\t")
                    correct_loc = data_split[1]
                    map_circle = self.add_tool("dna_gmap.marker_order")
                    map_circle.set_options({
                        "correct_loc": correct_loc,
                        "poptype": self.option("poptype"),
                        # "ref": self.option("ref"),
                        "ref": 'no',
                        "min_lod": self.option("min_lod"),
                        "randomise_order": self.option("randomise_order"),
                        "knn": self.option("knn"),
                        "bin": self.option('bin'),
                        "pri_marker_path": self.option('pri_marker_path')
                    })
                    self.marker_order2_tools.append(map_circle)
            for j in range(len(self.marker_order2_tools)):
                self.marker_order2_tools[j].on("end", self.set_output, "map_circle2")
            if self.marker_order2_tools:
                if len(self.marker_order2_tools) > 1:
                    self.on_rely(self.marker_order2_tools, self.map_circle3_run)
                elif len(self.marker_order2_tools) == 1:
                    self.marker_order2_tools[0].on('end', self.map_circle3_run)
                for tool in self.marker_order2_tools:
                    gevent.sleep(1)
                    tool.run()
            else:
                self.set_error("marker_order2_tools列表为空！", code="24800208")

    def map_circle3_run(self):
        """
        添加对第二次排图的结果进行检查，如果所有的map文件中有文件遗传距离都是0的时候，就直接用第一次的结果
        :return:
        """
        self.logger.info("22222:{}".format(self.marker_order2_tools[0].output_dir))
        if self.check_map_distance():
            self.end()
        else:
            if len(self.marker_order2_tools) == 1:
                self.get_list(self.marker_order2_tools[0].output_dir)
                marker_path = os.path.join(self.marker_order2_tools[0].output_dir, "marker.list")
            else:
                self.get_list(self.output_dir, "map_circle2")
                marker_path = os.path.join(self.output_dir, "map_circle2/marker.list")
            with open(marker_path, "r")as f1:
                lines = f1.readlines()
                for line in lines:
                    data_split = line.strip().split("\t")
                    correct_loc = data_split[1]
                    map_circle = self.add_tool("dna_gmap.marker_order")
                    map_circle.set_options({
                        "correct_loc": correct_loc,
                        "poptype": self.option("poptype"),
                        # "ref": self.option("ref"),
                        "ref": 'no',
                        "min_lod": self.option("min_lod"),
                        "randomise_order": self.option("randomise_order"),
                        "knn": self.option("knn"),
                        "bin": self.option('bin'),
                        "pri_marker_path": self.option('pri_marker_path')
                    })
                    self.marker_order3_tools.append(map_circle)
            for j in range(len(self.marker_order3_tools)):
                self.marker_order3_tools[j].on("end", self.set_output, "map_circle3")
            if self.marker_order3_tools:
                if len(self.marker_order3_tools) > 1:
                    self.on_rely(self.marker_order3_tools, self.end)
                elif len(self.marker_order3_tools) == 1:
                    self.marker_order3_tools[0].on('end', self.end)
                for tool in self.marker_order3_tools:
                    gevent.sleep(1)
                    tool.run()
            else:
                self.set_error("marker_order3_tools列表为空！", code="24800209")

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
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'map_circle1':
            self.linkdir(obj.output_dir, "map_circle1")
        if event['data'] == 'map_circle2':
            self.linkdir(obj.output_dir, "map_circle2")
        if event['data'] == 'map_circle3':
            self.linkdir(obj.output_dir, "map_circle3")

    def run(self):
        super(MarkerOrderModule, self).run()
        self.map_circle1_run()

    def end(self):
        if len(self.marker_order3_tools) == 1:
            self.option("map_cycle_dir").set_path(self.marker_order3_tools[0].output_dir)
        else:
            self.option("map_cycle_dir").set_path(self.output_dir + "/map_circle3")
        super(MarkerOrderModule, self).end()

    def check_hkxhk(self):
        """
        用于检查run_crosslink_group运行结果中的*.000.loc 是不是都是hkxhk 如果都是这个标记的话 就跳过run_smooth_cp
        :return:
        """
        is_hkxhk = False
        if self.marker_order1_tools == 1:
            file_path = self.marker_order1_tools[0].output_dir
        else:
            file_path = self.output_dir + "/map_circle1"
        for file_ in os.listdir(file_path):
            if re.match('.*\.correct\.loc$', file_):
                hkxhk_num = 0
                with open(os.path.join(file_path, file_), "r") as r:
                    data = r.readlines()
                    for line in data:
                        # temp = line.strip().split('\t')
                        temp = re.split('[ \t]', line)
                        if temp[1] == "<hkxhk>":
                            hkxhk_num += 1
                    if hkxhk_num == len(data):
                        self.logger.info("{}中标记全部为hkxhk！只进行一次排图计算！")
                        is_hkxhk = True
                        break
        if is_hkxhk:
            self.linkdir(file_path, "map_circle3")
        return is_hkxhk

    def check_map_distance(self):
        map_distance_is_zero = False
        if self.marker_order2_tools == 1:
            file_path = self.marker_order2_tools[0].output_dir
        else:
            file_path = self.output_dir + "/map_circle2"
        for map_file in os.listdir(file_path):
            if re.match('.*\.female\.map$', map_file) or re.match('.*\.male\.'
                                                                  'map$', map_file) or re.match('.*\.sexAver\.'
                                                                                                'map$', map_file):
                with open(os.path.join(file_path, map_file), 'r') as r:
                    is_zero = True
                    for line in r:
                        if re.match(r'^group.*', line):
                            pass
                        else:
                            temp = line.strip().split('\t')
                            try:
                                map_distance = float(temp[1])
                            except:
                                continue
                            if map_distance > 0:
                                is_zero = False
                    if is_zero:
                        self.logger.info("文件{}遗传距离全部为0".format(map_file))
                        map_distance_is_zero = True
                        break
        if map_distance_is_zero:
            self.logger.info("有map文件的遗传距离为0，所以直接将第一次排图结果作为第三次结果用于后面计算!")
            self.linkdir(self.output_dir + "/map_circle1", "map_circle3")
        return map_distance_is_zero

