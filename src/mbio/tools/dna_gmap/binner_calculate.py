# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.11

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class BinnerCalculateAgent(Agent):
    """
    遗传图谱：BinMarker分析
    """
    def __init__(self, parent):
        super(BinnerCalculateAgent, self).__init__(parent)
        options = [
            {"name": "pop_final_vcf", "type": "infile", "format": "bsa.vcf"},  # pop.final.vcf.gz文件，用于染色体bin统计
            {"name": "genotype_matrix", "type": "infile", "format": "dna_gmap.marker"},  # 分型矩阵，此处为pop.filtered.marker,遗传标记筛选的结果文件
            {"name": "marker_upload", "type": "infile", "format": "dna_gmap.marker"},  # 客户上传的标记列表
            {"name": "pop_type", "type": "string"},  # 群体类型，F2，F1，BC，CP，RIL  有CP是F1，无CP是F2,
            {"name": "window_size", "type": "float", "default": 2000},  # window size(kb)
            {"name": "window_step", "type": "float", "default": 100},  # step size(kb)
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "total_bin_marker", "type": "outfile", "format": "dna_gmap.marker"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("pop_final_vcf").is_set:
            raise OptionError("请设置pop.final.vcf文件", code="34800101")
        if not self.option("genotype_matrix").is_set:
            raise OptionError("请设置分型矩阵genotype_matrix", code="34800102")
        if not self.option("pop_type"):
            raise OptionError("请设置群体类型", code="34800103")
        if self.option("pop_type") not in ["F1", "F2", "BC", "CP", "RIL", "DH"]:
            raise OptionError("群体类型只能在F1/F2/BC/CP/RIL/DH内", code="34800104")

    def set_resource(self):
        self._cpu = 3
        if self.option("genotype_matrix").prop["marker_num"] < 100000:
            self._memory = "20G"
        elif self.option("genotype_matrix").prop["marker_num"] < 1000000:
            self._memory = "60G"
        else:
            self._memory = "100G"

    def end(self):
        super(BinnerCalculateAgent, self).end()


class BinnerCalculateTool(Tool):
    def __init__(self, config):
        super(BinnerCalculateTool, self).__init__(config)
        self.python_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/python"
        self.perl_path = "miniconda2/bin/perl"
        self.perl = self.config.SOFTWARE_DIR + "/miniconda2/bin/perl"
        self.nocp_pesudo = self.config.PACKAGE_DIR + "/dna_gmap/binNOCP-pesudo.pl"
        self.cp_pesudo = self.config.PACKAGE_DIR + "/dna_gmap/binCP-pesudo.pl"
        self.binmap_wintype = self.config.PACKAGE_DIR + "/dna_gmap/binmap-wintype.pl"
        self.binmap_phase = self.config.PACKAGE_DIR + "/dna_gmap/binmap-phase.pl"
        self.nocp_merge = self.config.PACKAGE_DIR + "/dna_gmap/binNOCP-merge.pl"
        self.cp_merge = self.config.PACKAGE_DIR + "/dna_gmap/binCP-merge.pl"
        # self.binner_chr = self.config.PACKAGE_DIR + "/dna_gmap/binner_chr.pl"
        self.binner_chr = self.config.PACKAGE_DIR + "/dna_gmap/binner_chr.py"
        self.binner_pos = self.config.PACKAGE_DIR + "/dna_gmap/binner_pos.py"
        self.binner_info = self.config.PACKAGE_DIR + "/dna_gmap/binner_info1.pl"
        self.para_fly = 'program/parafly-r2013-01-21/bin/bin/ParaFly'  # 并行投递cmd

    def filter_marker(self):
        """
        根据客户上传的标记列表过滤分型矩阵
        """
        self.logger.info("根据上传的标记列表进行分分型矩阵的过滤")
        marker_ids = []
        with open(self.option("marker_upload").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                marker_ids.append(item[0])
        genotype_matrix = os.path.join(self.work_dir, "new_genotype_matrix.xls")
        with open(self.option("genotype_matrix").prop["path"], "r") as f, open(genotype_matrix, "w") as w:
            lines = f.readlines()
            w.write(lines[0])
            for line in lines[1:]:
                item = line.strip().split("\t")
                if item[0] not in marker_ids:
                    w.write(line)
        return genotype_matrix

    def run_bincp_pesudo(self):
        """
        binCP-pesudo.pl
        """
        cmd = "{} {} -i {}".format(self.perl_path, self.cp_pesudo, self.genotype_matrix)
        cmd += " -o {} -k {}".format(self.work_dir, "Total")
        self.run_cmd(cmd, "bin_cp_pesub")

    def run_binnocp_pesudo(self):
        """
        binNOCP-pesudo.pl
        """
        cmd = "{} {} -i {}".format(self.perl_path, self.nocp_pesudo, self.genotype_matrix)
        cmd += " -o {} -k {}".format(self.work_dir, "Total")
        self.run_cmd(cmd, "bin_nocp_pesub")

    def run_binmap_wintype(self):
        """
        male/female: binmap-wintype.pl
        """
        cmd_list = []
        male_matrix = os.path.join(self.work_dir, "Total.male.matrix")
        male_cmd = "{} {} -i {} -o {}".format(self.perl, self.binmap_wintype, male_matrix, self.work_dir)
        male_cmd += " -k {} -win {} -step {}".format("Total.male", self.option("window_size"), self.option("window_step"))
        cmd_list.append(male_cmd)
        female_matrix = os.path.join(self.work_dir, "Total.female.matrix")
        female_cmd = "{} {} -i {} -o {}".format(self.perl, self.binmap_wintype, female_matrix, self.work_dir)
        female_cmd += " -k {} -win {} -step {}".format("Total.female", self.option("window_size"), self.option("window_step"))
        cmd_list.append(female_cmd)
        self.run_cmd_more(cmd_list, "binmap_wintype")

    def run_binmap_phase(self):
        """
        male/female: binmap-phase.pl
        """
        cmd_list = []
        male_wintype = os.path.join(self.work_dir, "Total.male.wintype")
        male_cmd = "{} {} -i {} -o {} -k {}".format(self.perl, self.binmap_phase, male_wintype, self.work_dir, "Total.male")
        cmd_list.append(male_cmd)
        female_wintype = os.path.join(self.work_dir, "Total.female.wintype")
        female_cmd = "{} {} -i {} -o {} -k {}".format(self.perl, self.binmap_phase, female_wintype, self.work_dir, "Total.female")
        cmd_list.append(female_cmd)
        self.run_cmd_more(cmd_list, "binmap_phase")

    def run_binnocp_merge(self):
        """
        binNOCP-merge.pl
        """
        female_bin_phase = os.path.join(self.work_dir, "Total.female.bin.phase")
        male_bin_phase = os.path.join(self.work_dir, "Total.male.bin.phase")
        cmd = "{} {} -f {} -m {}".format(self.perl_path, self.nocp_merge, female_bin_phase, male_bin_phase)
        cmd += " -o {} -k {}".format(self.work_dir, "Total")
        self.run_cmd(cmd, "bin_nocp_merge")

    def run_bincp_merge(self):
        """
        binCP-merge.pl
        """
        female_bin_phase = os.path.join(self.work_dir, "Total.female.bin.phase")
        male_bin_phase = os.path.join(self.work_dir, "Total.male.bin.phase")
        cmd = "{} {} -f {} -m {}".format(self.perl_path, self.cp_merge, female_bin_phase, male_bin_phase)
        cmd += " -o {} -k {}".format(self.work_dir, "Total")
        self.run_cmd(cmd, "bin_cp_merge")

    def run_bin_info(self):
        """
        binner_info1.pl
        """
        bin_marker = os.path.join(self.work_dir, "Total.bin.marker")
        info_cmd = "{} {} -popt {} -bin {}".format(self.perl_path, self.binner_info, self.option("pop_type"), bin_marker)
        info_cmd += " -mark {} -out {}".format(self.genotype_matrix, self.output_dir + "/bin_info.xls")
        info_cmd += " -out_pos {}".format(self.work_dir + "/bin_detail.xls")
        self.run_cmd(info_cmd, "bin_info")

    def run_binner_stat(self):
        """
        binner_chr.py/binner_pos.py
        """
        cmd_list = []
        bin_marker = os.path.join(self.work_dir, "Total.bin.marker")
        out_bin = os.path.join(self.output_dir, "Total.bin.marker")
        if os.path.exists(out_bin):
            os.remove(out_bin)
        os.link(bin_marker, out_bin)
        self.option("total_bin_marker", out_bin)
        chr_cmd = "{} {} -vcf {} -bin {}".format(self.python_path, self.binner_chr, self.option("pop_final_vcf").prop["path"], bin_marker)
        chr_cmd += " -mark {} -o {}".format(self.genotype_matrix, self.output_dir + "/bin_stat.xls")
        cmd_list.append(chr_cmd)
        pos_cmd = "{} {} -bin {} -vcf {}".format(self.python_path, self.binner_pos, self.work_dir + "/bin_detail.xls", self.option("pop_final_vcf").prop["path"])
        pos_cmd += " -o {}".format(self.output_dir + "/bin_pos.xls")
        cmd_list.append(pos_cmd)
        self.run_cmd_more(cmd_list, "bin_stat")

    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("%s运行失败", variables=(cmd_name), code="34800101")
            self.set_error("%s运行失败", variables=(cmd_name), code="34800104")

    def run_cmd_more(self, cmd_list, cmd_name):
        """
        将多个cmd命令并行执行
        """
        cmd_file = os.path.join(self.work_dir, "cmd_list_{}.txt".format(cmd_name))
        wrong_cmd = os.path.join(self.work_dir, "failed_cmd_{}.txt".format(cmd_name))
        with open(cmd_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.para_fly, cmd_file, 2, wrong_cmd)
        self.run_cmd(cmd_more, "more_" + cmd_name)

    def set_db(self):
        """
        将结果导入mongo数据库
        """
        self.logger.info("将结果导入mongo数据库")
        binner_api = self.api.api("dna_gmap.binner_calculate")
        binmarker_id = self.option("main_id")
        bin_stat = os.path.join(self.output_dir, "bin_stat.xls")
        bin_info = os.path.join(self.output_dir, "bin_info.xls")
        bin_pos = os.path.join(self.output_dir, "bin_pos.xls")
        binner_api.add_sg_binmarker_bin(binmarker_id, bin_stat)
        binner_api.add_sg_binmarker_var(binmarker_id, bin_info)
        binner_api.add_sg_binmarker_var_detail(binmarker_id, bin_pos)

    def run(self):
        super(BinnerCalculateTool, self).run()
        if self.option("marker_upload").is_set:
            self.genotype_matrix = self.filter_marker()
        else:
            self.genotype_matrix = self.option("genotype_matrix").prop["path"]
        if self.option("pop_type") == "F1":
            self.run_bincp_pesudo()
        else:
            self.run_binnocp_pesudo()
        self.run_binmap_wintype()
        self.run_binmap_phase()
        if self.option("pop_type") == "F1":
            self.run_bincp_merge()
        else:
            self.run_binnocp_merge()
        self.run_bin_info()
        self.run_binner_stat()
        if self.option("main_id"):
            self.set_db()
        self.end()
