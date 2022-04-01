# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.06.15

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class MarkerOrderAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(MarkerOrderAgent, self).__init__(parent)
        options = [
            {"name": "correct_loc", "type": "infile", "format": "dna_gmap.marker"},  # 传入split_markers跑出来的correct.loc文件/pri.marker文件。
            {"name": "poptype", "type": "string"},  # 传入poptype参数，例：CP,F2等
            {"name": "ref", "type": "string"},  # 判断是否为ref，传入参数为yes/no
            {"name": "min_lod", "type": "string"},
            {"name": "randomise_order", "type": "string"},
            {"name": "knn", "type": "string"},
            {"name": "bin", "type": "string", "default": "yes"},
            {"name": "pri_marker_path", 'type': 'string'},
            {"name": "is_first", "type":"string", "default": "false"},
            {"name": "ignore_cxr", "type": 'string', "default": "1"}
            # {"name": "marker", "type": "outfile", "format": "dna_gmap.marker"}  # 输出marker/loc文件。
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("correct_loc"):
            raise OptionError("必须输入correct_loc", code="34800801")
        if not self.option("poptype"):
            raise OptionError("必须输入poptype", code="34800802")
        if self.option("poptype") == "CP":
            if self.option("ref") not in ["yes", "no"]:
                raise OptionError("ref参数应该为yes或no", code="34800803")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(MarkerOrderAgent, self).end()


class MarkerOrderTool(Tool):
    def __init__(self, config):
        super(MarkerOrderTool, self).__init__(config)
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.map_gather_pl = self.config.PACKAGE_DIR + '/dna_gmap/map-gather.pl'
        self.smooth_CP_pl = self.config.PACKAGE_DIR + '/dna_gmap/smooth-CP.pl'
        self.smooth_NOCP_pl = self.config.PACKAGE_DIR + '/dna_gmap/smooth-NOCP.pl'
        self.crosslink_pos = 'bioinfo/gmap/crosslink_pos'
        self.crosslink_group = 'bioinfo/gmap/crosslink_group'
        self.crosslink_map = 'bioinfo/gmap/crosslink_map'
        self.MSTmap = 'bioinfo/gmap/MSTmap'
        self.chr_name = ''
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')

    def clean_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))

    def set_chr_name(self):
        pri_marker_list = self.option("correct_loc").prop["path"].strip().split('/')
        chr_name_list = pri_marker_list[len(pri_marker_list) - 1].strip().split('.')
        self.chr_name = chr_name_list[0]

    def run_crosslink_pos(self):
        # if self.option('is_first') == 'false':
        #     loc = os.path.join(self.output_dir, (self.chr_name + ".0.000.loc"))
        # else:
        loc = self.option("correct_loc").prop["path"]
        out_path = os.path.join(self.output_dir, (self.chr_name + ".000.map"))
        cmd = "{} --inp={} --out={}".format(self.crosslink_pos, loc, out_path)
        command = self.add_command("crosslink_pos", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("crosslink_pos完成")
        else:
            self.set_error("crosslink_pos失败", code="34800801")

    def run_crosslink_group(self):
        out_path = os.path.join(self.output_dir, (self.chr_name + ".0."))
        cmd = "{} --inp={} --outbase={} --mapbase={} --min_lod={} --randomise_order={} --knn={}"\
            .format(self.crosslink_group, self.option("correct_loc").prop["path"], out_path, out_path,
                    self.option("min_lod"), self.option("randomise_order"), self.option("knn"))
        command = self.add_command("crosslink_group", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("crosslink_group完成")
        else:
            self.set_error("crosslink_group失败", code="34800802")

    def crosslink_group_run(self):
        """
        这里只针对f1群体，同时是bin的时候
        ref:yes
        poptype:CP
        crosslink_group
        """
        loc_name = os.path.join(self.output_dir, (self.chr_name + ".0."))
        cmd = "{} --inp={} --outbase={} --min_lod={} --ignore_cxr={} --randomise_order={} --knn=3" \
            .format(self.crosslink_group, self.option("correct_loc").prop["path"], loc_name, self.option("min_lod"),
                    self.option("ignore_cxr"), self.option("randomise_order"))
        command = self.add_command("crosslink_group_1", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("crosslink_group完成")
        else:
            self.set_error("crosslink_group失败", code="34800810")

    def run_crosslink_map(self):
        inp = os.path.join(self.output_dir, (self.chr_name + ".0.000.loc"))
        out = os.path.join(self.output_dir, (self.chr_name + ".000.loc"))
        map_path = os.path.join(self.output_dir, (self.chr_name + ".000.map"))
        cmd = "{} --inp={} --out={} --map={}".format(self.crosslink_map, inp, out, map_path)
        command = self.add_command("crosslink_map", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("crosslink_map完成")
        else:
            self.set_error("crosslink_map失败", code="34800803")

    def run_map_gather(self):
        """
        crosslink_pos/crosslink_group中可共用。
        """
        map_path = os.path.join(self.output_dir, (self.chr_name + ".000.map"))
        out_path = os.path.join(self.output_dir, self.chr_name)
        cmd = "{} {} -i {} -o {}".format(self.perl_path, self.map_gather_pl, map_path, out_path)
        command = self.add_command("map_gather", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("map_gather完成")
        else:
            self.set_error("map_gather失败", code="34800804")

    def run_smooth_cp(self):
        """
        ref:yes/ref:no传入参数不同。 注意这里要是是无参的话要添加loc_path文件
        """
        map_path = os.path.join(self.output_dir, (self.chr_name + ".sexAver.map"))
        loc_path = os.path.join(self.output_dir, (self.chr_name + ".000.loc"))
        # if self.option('bin') == 'yes':
        #     ref_loc = os.path.join(self.option('pri_marker_path'), (self.chr_name + ".pri.marker"))
        # else:
        ref_loc = self.option("correct_loc").prop["path"]
        if self.option("ref") == "yes":
            cmd = "{} {} -m {} -l {} -k {} -d {}".format(self.perl_path, self.smooth_CP_pl, map_path,
                                                         ref_loc, self.chr_name, self.output_dir)
        else:
            cmd = "{} {} -m {} -l {} -k {} -d {}".format(self.perl_path, self.smooth_CP_pl, map_path,
                                                         loc_path, self.chr_name, self.output_dir)
        command = self.add_command("smooth_cp", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("smooth_cp完成")
        else:
            self.set_error("smooth_cp失败", code="34800805")

    def run_mstmap(self):
        out_path = os.path.join(self.output_dir, (self.chr_name + ".out"))
        cmd = "{} {} {}".format(self.MSTmap, self.option("correct_loc").prop["path"], out_path)
        command = self.add_command("mstmap", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("MSTmap完成")
        else:
            self.set_error("MSTmap失败", code="34800806")

    def run_smooth_nocp(self):
        out_path = os.path.join(self.output_dir, (self.chr_name + ".out"))
        correct_loc = os.path.join(self.output_dir, (self.chr_name + ".correct.marker"))
        cmd = "{} {} -i {} -m {} -o {}".format(self.perl_path, self.smooth_NOCP_pl,
                                               self.option("correct_loc").prop["path"], out_path, correct_loc)
        command = self.add_command("smooth_nocp", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("smooth_nocp完成")
        else:
            self.set_error("smooth_nocp失败", code="34800807")

    def run(self):
        super(MarkerOrderTool, self).run()
        self.clean_output()
        self.set_chr_name()
        if self.option('poptype') == 'CP':
            if self.option('ref') == 'yes':
                # if self.option('is_first') == 'false':
                #     self.crosslink_group_run()
                #     # self.run_crosslink_map()
                # else:
                #     pass
                self.run_crosslink_pos()
                self.run_map_gather()
                self.run_smooth_cp()
                self.end()
            else:
                self.run_crosslink_group()
                self.run_crosslink_map()
                self.run_map_gather()
                self.run_smooth_cp()
                self.end()
        else:
            self.run_mstmap()
            self.run_smooth_nocp()
            self.end()
