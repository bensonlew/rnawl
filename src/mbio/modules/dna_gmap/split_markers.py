# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# last_modify:20180614

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os


class SplitMarkersModule(Module):
    """
    用于包splitbylg.py和smooth.py两个tool。
    """
    def __init__(self, work_id):
        super(SplitMarkersModule, self).__init__(work_id)
        options = [
            {"name": "bin", "type": "string"},  # 判断是否为bin，传入参数为yes/no
            {"name": "poptype", "type": "string"},  # 传入poptype参数，例：CP,F2等
            {"name": "marker_path", "type": "infile", "format": "dna_gmap.marker"},
            # 传入marker文件，没bin：pop.filtered.marker；有bin：Total.bin.marker
            {"name": "lg_path", "type": "infile", "format": "dna_gmap.lg"},  # 传入Total.lg
            {"name": "ref", "type": "string"},  # 判断是否为ref，传入参数为yes/no
            {"name": "min_lod", "type": "string", "default": "-1"},  # crosslink_group参数。
            {"name": "ignore_cxr", "type": "string", "default": "1"},  # crosslink_group参数。
            {"name": "randomise_order", "type": "string", "default": "0"},  # crosslink_group参数。
            {"name": "marker_list", "type": "outfile", "format": "dna_gmap.marker_list"},
            {'name': "pri_marker_path", "type": "outfile", "format": "dna_gmap.marker_path"}
        ]
        self.add_option(options)
        self.splitbylg = self.add_tool("dna_gmap.splitbylg")
        self.smooth_tools = []

    def check_options(self):
        if not self.option("bin"):
            raise OptionError("缺少bin参数", code="24800401")
        if not self.option("poptype"):
            raise OptionError("缺少poptype参数", code="24800402")
        if not self.option("marker_path"):
            raise OptionError("缺少marker_path参数", code="24800403")
        if not self.option("lg_path"):
            raise OptionError("缺少lg_path参数", code="24800404")
        if self.option("ref") not in ["yes", "no"]:
            raise OptionError("ref参数应该为yes或no", code="24800405")
        if not self.option("min_lod"):
            raise OptionError("缺少min_lod参数", code="24800406")
        if not self.option("ignore_cxr"):
            raise OptionError("缺少ignore_cxr参数", code="24800407")
        if not self.option("randomise_order"):
            raise OptionError("缺少randomise_order参数", code="24800408")
        return True

    def get_name(self):
        """
        获取结果文件中的染色体名称。
        生成name.list文件，后续重跑时方便查看顺序。
        """
        name_list = []
        for m in os.listdir(self.splitbylg.option("file_path").prop["path"]):
            if re.match(r".*\.pri\.map$", m):
                if m not in ['.pri.marker', '.pri.map']:
                    map_path = m.strip().split(".")
                    name_list.append(map_path[0])
        with open(self.work_dir + "/name.list", "w") as fw:
            for i in range(len(name_list)):
                fw.write(name_list[i] + "\n")

    def splitbylg_run(self):
        self.splitbylg.set_options({
            "bin": self.option("bin"),
            "poptype": self.option("poptype"),
            "marker_path": self.option("marker_path").prop["path"],
            "lg_path": self.option("lg_path").prop["path"]
        })
        if self.option("ref") == "yes":
            self.splitbylg.on('end', self.smooth_run)
            self.splitbylg.run()
        elif self.option("ref") == "no":
            self.splitbylg.on("end", self.set_output, "splitbylg")
            self.splitbylg.on("end", self.end)
            self.splitbylg.run()
            self.option("marker_list").set_path(self.output_dir + "/pri.marker.list")

    def smooth_run(self):
        self.get_name()
        name_list_path = os.path.join(self.work_dir, "name.list")
        with open(name_list_path, "r") as fr:
            lines = fr.readlines()
            for line in lines:
                data_split = line.strip().split("\t")
                chr_name = data_split[0]
                pri_marker = os.path.join(self.splitbylg.option("file_path").prop["path"], (chr_name + ".pri.marker"))
                pri_map = os.path.join(self.splitbylg.option("file_path").prop["path"], (chr_name + ".pri.map"))
                self.logger.info(pri_marker)
                self.logger.info(pri_map)
                smooth = self.add_tool("dna_gmap.smooth")
                smooth.set_options({
                    "poptype": self.option("poptype"),
                    "pri_marker": pri_marker,
                    "pri_map": pri_map,
                    "min_lod": self.option("min_lod"),
                    "ignore_cxr": self.option("ignore_cxr"),
                    "randomise_order": self.option("randomise_order")
                })
                self.smooth_tools.append(smooth)
        for j in range(len(self.smooth_tools)):
            self.smooth_tools[j].on("end", self.set_output, 'smooth')
        if self.smooth_tools:
            if len(self.smooth_tools) > 1:
                self.on_rely(self.smooth_tools, self.ref_list)
            elif len(self.smooth_tools) == 1:
                self.smooth_tools[0].on('end', self.ref_list)
            for tool in self.smooth_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("smooth_tools列表为空！", code="24800403")

    def ref_list(self):
        """
        生成ref.marker.list
        """
        if len(self.smooth_tools) == 1:
            out_path = self.smooth_tools[0].output_dir
        else:
            out_path = self.output_dir
        with open(self.output_dir + "/ref.marker.list", "w") as f3:
            for root, dirs, files in os.walk(out_path):
                for name in files:
                    if re.match('.*\.correct\.loc$', name):
                        if name not in ['.correct.loc']:
                            loc_path = os.path.join(self.output_dir, name)
                            name_split = name.strip().split(".")
                            chr_name = name_split[0]
                            f3.write(chr_name + "\t" + loc_path + '\n')
        self.option("marker_list").set_path(self.output_dir + "/ref.marker.list")
        self.option("pri_marker_path").set_path(self.splitbylg.option("file_path").prop["path"])
        self.end()

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
        if event['data'] == 'splitbylg':
            self.linkdir(obj.output_dir, self.output_dir)
        if event['data'] == 'smooth':
            self.linkdir(obj.output_dir, self.output_dir)
        else:
            pass

    def run(self):
        super(SplitMarkersModule, self).run()
        self.splitbylg_run()

    def end(self):
        super(SplitMarkersModule, self).end()
