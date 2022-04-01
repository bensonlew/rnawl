# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 2018.06.14

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class MapEvaluationAgent(Agent):
    """
    遗传图谱：图谱评估分析
    lasted modified by hongdong @ 20180614
    """
    def __init__(self, parent):
        super(MapEvaluationAgent, self).__init__(parent)
        options = [
            {"name": "map_cycle_dir", "type": "infile", "format": "dna_gmap.map_cycle_dir"},  # 排图结果路径，这里后面加一个file文件
            {"name": "pop_type", "type": "string", "default": "F2"},  # 群体类型，F2，F1，BC，CP，RIL
            # 有CP是F1，无CP是F2
            {"name": "ref_chrlist", "type": "string"},  # 染色体列表
            {"name": "is_ref", "type": "string", "default": "false"},  # 是否是有参的
            {"name": "filtered_marker", "type": "infile", "format": "dna_gmap.marker"}  # 可以没有，
            # 上游cp是f1群体的时候才会有该参数
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("map_cycle_dir"):
            raise OptionError("请设置map_cycle_dir文件", code="34800601")
        if not self.option("ref_chrlist"):
            raise OptionError("请设置ref_chrlist文件", code="34800602")
        if self.option("pop_type").lower() not in ["f1", "f2", "bc", "cp", "ril", 'dh']:
            raise OptionError("群体类型只能在F1/F2/BC/CP/RIL内", code="34800603")
        if self.option("pop_type").lower() == 'f1' and not self.option("filtered_marker"):
            raise OptionError("群体类型是F1，必须要有filtered_marker参数！", code="34800604")

    def set_resource(self):
        self._cpu = 5
        self._memory = "50G"

    def end(self):
        super(MapEvaluationAgent, self).end()


class MapEvaluationTool(Tool):
    def __init__(self, config):
        super(MapEvaluationTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.perl = self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin/perl"
        self.mapmergenocp = self.config.PACKAGE_DIR + "/dna_gmap/MapMergeNOCP.pl"
        self.mapmergecp = self.config.PACKAGE_DIR + "/dna_gmap/MapMergeCP.pl"
        self.mapestimate = self.config.PACKAGE_DIR + "/dna_gmap/mapEstimate.pl"
        self.markerinfo = self.config.PACKAGE_DIR + "/dna_gmap/markerinfo.pl"
        self.merge_map_total = self.config.PACKAGE_DIR + "/dna_gmap/merge-map-total.pl"
        self.ralationmap = self.config.PACKAGE_DIR + "/dna_gmap/RalationMap.pl"
        self.para_fly = 'program/parafly-r2013-01-21/bin/bin/ParaFly'
        self.drawmap = self.config.PACKAGE_DIR + "/dna_gmap/drawmap.R"
        self.drawbincp = self.config.PACKAGE_DIR + "/dna_gmap/drawbinCP.R"
        self.drawbincpsexAver = self.config.PACKAGE_DIR + "/dna_gmap/drawbinCP-sexAver.R"
        self.drawbinnocp = self.config.PACKAGE_DIR + "/dna_gmap/drawbinNOCP.R"
        self.r_path = 'program/R-3.3.1/bin/Rscript'
        self.rgcc5_path = 'program/R-3.3.1_gcc5.1/bin/Rscript'
        self.result = ''
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"

    def run_mapmergenocp(self):
        """
        MapMergeNOCP.pl
        /mnt/ilustre/users/sanger-dev/app/program/perl/perls/perl-5.24.0/bin/perl
        /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/dna_gmap/MapMergeNOCP.pl
        -dmap ~/sg-users/xuanhongdong/gmap/08.map-cycle3_nobin/ -o nobin_result -adjust
        """
        cmd = "{} {} -dmap {} -o {} -adjust".format(self.perl_path, self.mapmergenocp,
                                                    self.option("map_cycle_dir").prop['path'], self.result)
        self.run_cmd(cmd, "map_merge_nocp")

    def run_mapmergecp(self):
        """
        MapMergeCP.pl
        /mnt/ilustre/users/sanger-dev/app/program/perl/perls/perl-5.24.0/bin/perl
        /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/dna_gmap/MapMergeCP.pl
        -dmap ~/sg-users/xuanhongdong/gmap/08.map-cycle3_nobin/ -o nobin_result -adjust
        """

        cmd = "{} {} -dmap {} -o {} -mark {} -adjust".format(self.perl_path, self.mapmergecp,
                                                             self.option("map_cycle_dir").prop['path'], self.result,
                                                             self.option("filtered_marker").prop['path'])
        self.run_cmd(cmd, "map_merge_cp")

    def run_mapestimate(self, pop_type):
        """
        mapEstimate.pl
        /mnt/ilustre/users/sanger-dev/app/program/perl/perls/perl-5.24.0/bin/perl
        /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/dna_gmap/mapEstimate.pl -i total.map -o total.mapstat
        """

        if pop_type.lower() not in ['f1', 'cp']:
            cmd = "{} {}".format(self.perl_path, self.mapestimate)
            cmd += " -i {} -o {}".format(self.result + '/total.map', self.result + "/total.mapstat")
            self.run_cmd(cmd, "mapestimate")
        else:
            cmd = "{} {}".format(self.perl, self.mapestimate)
            cmd_list = []
            sexaver = cmd + " -i {} -o {}".format(self.result + '/total.sexAver.map', self.result + '/sexAver.mapstat')
            cmd_list.append(sexaver)
            male = cmd + " -i {} -o {}".format(self.result + '/total.male.map', self.result + '/male.mapstat')
            cmd_list.append(male)
            female = cmd + " -i {} -o {}".format(self.result + '/total.female.map', self.result + '/female.mapstat')
            cmd_list.append(female)
            self.run_cmd_more(cmd_list, "mapestimate", 3)

    def run_markerinfo(self, pop_type):
        """
        markerinfo.pl
        /mnt/ilustre/users/sanger-dev/app/program/perl/perls/perl-5.24.0/bin/perl
        /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/dna_gmap/markerinfo.pl
        -map total.map -input total.marker --pop f2 -out total.marker.info
        当为cp的时候，input应该是total.loc但是这里生成的该文件为空，所以改成total.final.loc
        """
        if pop_type.lower() not in ['f1', 'cp']:
            cmd = "{} {} --pop {}".format(self.perl_path, self.markerinfo, pop_type.lower())
            cmd += " -map {} -input {} --out {}".format(self.result + '/total.map', self.result + '/total.marker',
                                                        self.result + '/total.marker.info')
            self.run_cmd(cmd, "markerinfo")
        else:
            cmd = "{} {} --pop {}".format(self.perl, self.markerinfo, pop_type.lower())
            cmd_list = []
            sexaver = cmd + " -map {} -input {} --out {}"\
                .format(self.result + '/total.sexAver.map', self.result + '/total.sexAver.loc',
                        self.result + '/total.sexAver.info')
            cmd_list.append(sexaver)
            male = cmd + " -map {} -input {} --out {}"\
                .format(self.result + '/total.male.map', self.result + '/total.sexAver.loc',
                        self.result + '/total.male.info')
            cmd_list.append(male)
            female = cmd + " -map {} -input {} --out {}"\
                .format(self.result + '/total.female.map', self.result + '/total.sexAver.loc',
                        self.result + '/total.female.info')
            cmd_list.append(female)
            self.run_cmd_more(cmd_list, "markerinfo", 3)

    def run_merge_map_total(self, pop_type):
        """
        merge-map-total.pl
        /mnt/ilustre/users/sanger-dev/app/program/perl/perls/perl-5.24.0/bin/perl
        /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/dna_gmap/merge-map-total.pl -list ref.chrlist
          -mark total.map -out 3-15.xls
        """

        if pop_type.lower() not in ['f1', 'cp']:
            cmd = "{} {} --list {}".format(self.perl_path, self.merge_map_total, self.option("ref_chrlist"))
            cmd += " -mark {} --out {}".format(self.result + '/total.map', self.result + '/total.map.result.xls')
            self.run_cmd(cmd, "merge_map_total")
        else:
            cmd = "{} {} --list {}".format(self.perl, self.merge_map_total, self.option("ref_chrlist"))
            cmd_list = []
            sexaver = cmd + " -mark {} --out {}".format(self.result + '/total.sexAver.map',
                                                        self.result + '/sexaver.map.result.xls')
            cmd_list.append(sexaver)
            male = cmd + " -mark {} --out {}".format(self.result + '/total.male.map',
                                                     self.result + '/male.map.result.xls')
            cmd_list.append(male)
            female = cmd + " -mark {} --out {}".format(self.result + '/total.female.map',
                                                       self.result + '/female.map.result.xls')
            cmd_list.append(female)
            self.run_cmd_more(cmd_list, "merge_map_total", 3)

    def run_ralationmap(self, pop_type):
        if pop_type.lower() not in ['f1', 'cp']:
            cmd = "{} {} -k {} -o {} -m {}".format(self.perl_path, self.ralationmap, "total.phy", self.output_dir,
                                                   self.result + '/total.map')
            self.run_cmd(cmd, "ralationmap")
        else:
            cmd = "{} {} -o {}".format(self.perl, self.ralationmap, self.output_dir)
            cmd_list = []
            sexaver = cmd + " -m {} -k {}".format(self.result + '/total.sexAver.map', "total.sexaver.phy")
            cmd_list.append(sexaver)
            male = cmd + " -m {} -k {}".format(self.result + '/total.male.map', "total.male.phy")
            cmd_list.append(male)
            female = cmd + " -m {} -k {}".format(self.result + '/total.female.map', "total.female.phy")
            cmd_list.append(female)
            self.run_cmd_more(cmd_list, "ralationmap", 3)

    def run_drawmap(self, pop_type):
        """
        Rscript drawmap.R --mark  /result/total --out /result/fig --pop cp
        Rscript $Bin/bin/drawmap.R --mark $out/total  --out $out/fig --pop $pop
        :param pop_type:
        :return:
        """
        if pop_type.lower() in ['f1', 'cp']:
            cmd = "{} {} --mark {} --out {} --pop cp".format(self.rgcc5_path, self.drawmap, self.result + "/total",
                                                             self.output_dir + "/fig")
            self.run_cmd(cmd, "drawmap")
        else:
            cmd = "{} {} --mark {} --out {} --pop {}".format(self.rgcc5_path, self.drawmap, self.result + "/total",
                                                             self.output_dir + "/fig", self.option("pop_type"))
            self.run_cmd(cmd, "drawmap")

    def run_drawbin_heatmap(self, pop_type):
        """
        Rscript drawbinCP-sexAver.R --mark total.sexAver.phase --out ./total.sexAver.bin
        Rscript drawbinCP.R --mark total.female.phase --out ./total.female.bin
        Rscript drawbinNOCP.R --mark total.csv --out ./total.bin
        :param pop_type:
        :return:
        """
        heatmap_path = self.output_dir + "/fig1"
        if not os.path.exists(heatmap_path):
            os.mkdir(heatmap_path)
        if pop_type.lower() in ['f1', 'cp']:
            sexaver = "{} {} --mark {} --out {}".format(self.rgcc5_path, self.drawbincpsexAver,
                                                        self.result + '/total.sexAver.phase',
                                                        heatmap_path + "/total.sexAver.bin")
            self.run_cmd(sexaver, "sexaver_map")
            female = "{} {} --mark {} --out {}".format(self.rgcc5_path, self.drawbincp, self.result + '/total.female.phase',
                                                       heatmap_path + "/total.female.bin")
            self.run_cmd(female, "female_map")
            male = "{} {} --mark {} --out {}".format(self.rgcc5_path, self.drawbincp, self.result + '/total.male.phase',
                                                     heatmap_path + "/total.male.bin")
            self.run_cmd(male, "male_map")
        else:
            total = "{} {} --mark {} --out {}".format(self.rgcc5_path, self.drawbinnocp, self.result + '/total.csv',
                                                      heatmap_path + "/total.bin")
            self.run_cmd(total, "male_map")
            self.make_pdf(heatmap_path)

    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("%s运行失败", variables=(cmd_name), code="34800601")

    def create_dir(self, path):
        """
        检查目录在不在，不在的时候，直接创建，存在的时候，要先删除后再创建
        :param path:
        :return:
        """
        if not os.path.exists(path):
            os.mkdir(path)
            self.logger.info("创建文件夹{}成功！".format(path))
        else:
            os.rmdir(path)
            os.mkdir(path)
            self.logger.info("删除文件夹，并创建{}成功！".format(path))
        return path

    def make_pdf(self, file_path):
        """png转pdf"""
        cmd_list = []
        for m in os.listdir(file_path):
            n = re.match('(.*)\.png$', m)
            if n:
                png_path = os.path.join(self.output_dir + "/fig1/{}".format(m))
                pdf_path = os.path.join(self.output_dir + "/fig1/{}.pdf".format(n.group(1)))
                cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white {} {}'\
                    .format(png_path, pdf_path)
                cmd_list.append(cmd)
        self.run_cmd_more(cmd_list, "pdf2png", 3)

    def run_cmd_more(self, cmd_list, cmd_name, cpu):
        """
        将多个cmd命令并行执行
        """
        cmd_file = os.path.join(self.work_dir, "cmd_list_{}.txt".format(cmd_name))
        wrong_cmd = os.path.join(self.work_dir, "failed_cmd_{}.txt".format(cmd_name))
        with open(cmd_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.para_fly, cmd_file, cpu, wrong_cmd)
        self.run_cmd(cmd_more, "more_" + cmd_name)

    def run(self):
        super(MapEvaluationTool, self).run()
        # self.result = self.create_dir(self.work_dir + '/map_result')
        self.result = self.output_dir
        if self.option("pop_type").lower() == "f1":
            self.run_mapmergecp()
        else:
            self.run_mapmergenocp()
        self.run_mapestimate(self.option("pop_type").lower())
        self.run_markerinfo(self.option("pop_type").lower())
        # self.run_merge_map_total(self.option("pop_type").lower())
        if self.option("is_ref") == "true":  # 有参的时候执行相关性
            self.run_ralationmap(self.option("pop_type").lower())
        self.run_drawmap(self.option("pop_type").lower())
        self.run_drawbin_heatmap(self.option("pop_type").lower())
        self.end()
