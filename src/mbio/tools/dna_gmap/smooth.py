# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.06.12

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class SmoothAgent(Agent):
    """
    lasted modified by hongdong@20180731 添加了对hkxhk标记的判断，如果全部为该标记的话就不进行smooth-cp计算
    lasted modified by hongdong@20190103 当出现21.000.loc  21.001.loc两个分群的结果的时候，这个时候判断后并取行数最多的
    一个文件用于后面计算
    """
    def __init__(self, parent):
        super(SmoothAgent, self).__init__(parent)
        options = [
            {"name": "poptype", "type": "string"},  # 传入poptype参数，例：CP,F2等
            {"name": "pri_marker", "type": "infile", "format": "dna_gmap.marker"},  # 传入 染色体名称.pri.marker 文件
            {"name": "pri_map", "type": "infile", "format": "dna_gmap.map"},  # 传入 染色体名称.pri.map 文件
            {"name": "min_lod", "type": "string"},  # crosslink_group参数。
            {"name": "ignore_cxr", "type": "string"},  # crosslink_group参数。
            {"name": "randomise_order", "type": "string"}  # crosslink_group参数。
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("poptype"):
            raise OptionError("必须输入poptype", code="34801201")
        if not self.option("pri_marker"):
            raise OptionError("必须输入pri_marker", code="34801202")
        if not self.option("pri_map"):
            raise OptionError("必须输入pri_map", code="34801203")
        if not self.option("min_lod"):
            raise OptionError("必须输入min_lod", code="34801204")
        if not self.option("ignore_cxr"):
            raise OptionError("必须输入ignore_cxr", code="34801205")
        if not self.option("randomise_order"):
            raise OptionError("必须输入randomise_order", code="34801206")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(SmoothAgent, self).end()


class SmoothTool(Tool):
    def __init__(self, config):
        super(SmoothTool, self).__init__(config)
        self.perl_path = 'miniconda2/bin/perl'
        self.smooth_CP_pl = self.config.PACKAGE_DIR + '/dna_gmap/smooth-CP.pl'
        self.smooth_NOCP_pl = self.config.PACKAGE_DIR + '/dna_gmap/smooth-NOCP.pl'
        self.crosslink_group = 'bioinfo/gmap/crosslink_group'
        self.chr_name = ''
        self.loc_file = ''

    def clean_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))

    def set_chr_name(self):
        pri_marker_list = self.option("pri_marker").prop["path"].strip().split('/')
        chr_name_list = pri_marker_list[len(pri_marker_list) - 1].strip().split('.')
        self.chr_name = chr_name_list[0]

    def run_smooth_nocp(self):
        """
        ref:yes
        poptype:NOCP
        smooth-NOCP
        """
        loc_name = os.path.join(self.output_dir, (self.chr_name + ".correct.loc"))
        cmd = "{}  {} -i {} -o {} -m {}".format(self.perl_path, self.smooth_NOCP_pl,
                                                self.option("pri_marker").prop["path"],
                                                loc_name, self.option("pri_map").prop["path"])
        command = self.add_command("smooth_nocp", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("smooth_nocp完成")
        else:
            self.set_error("smooth_nocp失败", code="34801201")

    def run_crosslink_group(self):
        """
        ref:yes
        poptype:CP
        crosslink_group
        """
        loc_name = os.path.join(self.output_dir, (self.chr_name + "."))
        cmd = "{} --inp={} --outbase={} --min_lod={} --ignore_cxr={} --randomise_order={} --knn=25"\
            .format(self.crosslink_group, self.option("pri_marker").prop["path"], loc_name, self.option("min_lod"),
                    self.option("ignore_cxr"), self.option("randomise_order"))
        command = self.add_command("crosslink_group", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("crosslink_group完成")
        else:
            self.set_error("crosslink_group失败", code="34801202")

    def run_smooth_cp(self):
        """
        ref:yes
        poptype:CP
        smooth-CP
        """
        # loc_name = os.path.join(self.output_dir, (self.chr_name + ".000.loc"))
        loc_name = self.loc_file
        cmd = "{} {} -l {} -m {} -k {} -d {}".format(self.perl_path, self.smooth_CP_pl, loc_name,
                                                     self.option("pri_map").prop["path"], self.chr_name, self.output_dir)
        command = self.add_command("smooth_cp", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("smooth_cp完成")
        else:
            self.set_error("smooth_cp失败", code="34801203")

    def run(self):
        super(SmoothTool, self).run()
        self.clean_output()
        self.set_chr_name()
        if self.option('poptype') == 'CP':
            self.run_crosslink_group()
            self.check_loc_file()
            if not self.check_hkxhk():  # 返回为True的时候说明都是hkxhk标记
                self.run_smooth_cp()
            self.end()
        else:
            self.run_smooth_nocp()
            self.end()

    def check_hkxhk(self):
        """
        用于检查run_crosslink_group运行结果中的*.000.loc 是不是都是hkxhk 如果都是这个标记的话 就跳过run_smooth_cp
        :return:
        """
        # file_name = os.path.join(self.output_dir, (self.chr_name + ".000.loc"))
        file_name = self.loc_file
        hkxhk_num = 0
        with open(file_name, "r") as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split(' ')
                if temp[1] == "<hkxhk>":
                    hkxhk_num += 1
            if hkxhk_num == len(data):
                self.logger.info("{}中标记全部为hkxhk！".format(file_name))
                os.system('cp {} {}'.format(file_name, os.path.join(self.output_dir, (self.chr_name + ".correct.loc"))))
                return True
            else:
                return False

    def get_file_len(self, file_path):
        """
        获取文件行数
        :return:
        """
        return len(open(file_path, 'rU').readlines())

    def check_loc_file(self):
        """
        用于检查run_crosslink_group运行后有多少个loc文件，并取里面标记数最多的文件用于后面smooth_cp计算
        :return:
        """
        self.loc_file = self.output_dir + '/{}.000.loc'.format(self.chr_name)
        lines = self.get_file_len(self.loc_file)
        for m in os.listdir(self.output_dir):
            if re.match('.*\.([0-9]{3,})\.loc$', m):   # 21.001.loc 21.000.loc
                if m != '{}.000.loc'.format(self.chr_name):
                    if self.get_file_len(self.output_dir + '/{}'.format(m)) > lines:
                        self.loc_file = self.output_dir + '/{}'.format(m)
        self.logger.info("最终的loc文件为：{}".format(self.loc_file))
