# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last_modify: 2018.08.15

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
from bson.objectid import ObjectId


class RegionAnnoAgent(Agent):
    """
    基因注释接口
    一个controller和一个tool,一个api
    对应脚本cmd: region-gene.pl cmd: eff-enrich.R
    region_path,前三列只要是chrid start end就可以，位置信息有重复也没关系
    """
    def __init__(self, parent):
        super(RegionAnnoAgent, self).__init__(parent)
        options = [
            {"name": "main_id", "type": "string"},  #
            {"name": "task_id", "type": "string"},  #
            {"name": "pop_summary_path", "type": "infile", "format": "dna_evolution.pop_summary"},
            # {"name": "pop_summary_path", "type": "string"},     # 已经是绝对路径，待写成infile核查
            {"name": "region_path", "type": "infile", "format": "dna_evolution.bed"},    # sg_anno_params的region_path
            {"name": "update_info", "type": "string"},
            {"name": "go_summary_path", "type": "infile", "format": "dna_evolution.pop_enrich"},
            {"name": "genome_version_id", "type": "string"},
            {"name": "pathway_path", "type": "string"}  # 对象存储路径
        ]
        self.add_option(options)
        self.step.add_steps("region_anno")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.region_anno.start()
        self.step.update()

    def stepfinish(self):
        self.step.region_anno.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('pop_summary_path'):
            raise OptionError('必须输入pop_summary_path')
        if not self.option('region_path'):
            raise OptionError('必须输入region_path')
        if not self.option('main_id'):
            raise OptionError('必须输入main_id')
        if not self.option('task_id'):
            raise OptionError('必须输入task_id')
        if not self.option('genome_version_id'):
            raise OptionError('必须输入genome_version_id')
        if not self.option("pathway_path"):
            raise OptionError("请设置pathway_path")

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "6G"


class RegionAnnoTool(Tool):
    def __init__(self, config):
        super(RegionAnnoTool, self).__init__(config)
        self._version = "v1.0.1"
        self.perl_path = 'miniconda2/bin/perl '
        self.R_path = 'program/R-3.3.1/bin/Rscript '
        self.region_gene_path = self.config.PACKAGE_DIR + '/dna_evolution/region-gene.pl'
        self.effR_path = self.config.PACKAGE_DIR + '/dna_evolution/eff-enrich.R'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.ko_path = self.config.SOFTWARE_DIR + '/database/KEGG/kegg_2017-05-01/kegg/pathway/ko/'
        self.parafly = "program/parafly-r2013-01-21/src/ParaFly"

    # def get_select_region(self):
    #     """
    #     把传进来的region还原成文件
    #     eg：chr1-123-234;chr2-345-567;
    #     """
    #     region_list = self.option('select_region').split(';')
    #     with open(self.work_dir + '/select.region', 'w') as fw:
    #         fw.write('#chr\tstart\tend\n')
    #         for i in region_list:
    #             line_list = i.split('-')
    #             fw.write('\t'.join(line_list) + '\n')
    #     self.logger.info('~~~~~~select.region文件生成成功')

    def run_region_gene(self):
        """
        输入pop.summary select(三列信息，chrid start end)
        运行eff-enrich.R,计算注释统计
        """
        cmd = "{} {}".format(self.perl_path, self.region_gene_path)
        cmd += " -a {} -i {}".format(self.option('pop_summary_path').prop['path'],
                                     self.option('region_path').prop['path'])
        cmd += " -o {}".format(self.work_dir + "/select_region")
        self.logger.info('开始运行region_gene提取')
        self.run_cmd(cmd, "region_gene")

    def run_eff_detail(self, type_):
        """
        """
        cmd = "{} {}".format(self.R_path, self.effR_path)
        cmd += " --input {}".format(self.work_dir + "/select_region." + type_ + ".stat")
        cmd += " --output {}".format(self.output_dir + "/select_region." + type_ + ".stat")
        if type_ == 'eggnog':
            cmd += " --eggnog"
        else:
            cmd += " --top 1"
        self.logger.info('开始运行effrich.R注释统计:' + type_)
        self.run_cmd(cmd, "eff_detail_" + type_)

    def run_get_picture(self):
        """
        根据ko号找到对应的所有ko的pdf和png文件-- 添加如果存在pathway_dir就将该文件夹删除 add by hongdong 20180301
        :return:
        """
        if os.path.isdir(self.output_dir + '/pathway_dir'):
            os.system('rm -r %s' % os.path.join(self.output_dir, '/pathway_dir'))
        self.num = 0
        cmd_list = []
        os.system("mkdir {}".format(self.output_dir + '/pathway_dir'))
        all_ko = os.listdir(self.ko_path)
        with open(self.work_dir + "/select_region.kegg.stat")as fr:
            for line in fr:
                tmp = line.strip().split('\t')
                ko_png = tmp[0] + '.png'
                if ko_png in all_ko:
                    ko_path = os.path.join(self.ko_path, ko_png)
                    os.link(ko_path, self.output_dir + '/pathway_dir/' + ko_png)
                    pdf_path = os.path.join(self.output_dir + '/pathway_dir', tmp[0] + '.pdf')
                    cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + ko_path + ' ' + pdf_path
                    cmd_list.append(cmd)
                else:
                    self.set_error("ko{} png文件在本地数据库{}中没有，请核实!".format(tmp[0], self.ko_path))
        with open(self.work_dir + "/cmd.list", "w") as fw:
            for i in range(len(cmd_list)):
                fw.write(cmd_list[i] + "\n")
        cmd = self.parafly + " -c {} -CPU 10".format(self.work_dir + "/cmd.list")
        cmd_obj = self.add_command("cmd_list", cmd).run()
        self.wait(cmd_obj)
        if cmd_obj.return_code == 0:
            self.logger.info("cmd list执行成功")
        else:
            self.set_error("cmd list运行出错!")

    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("{}运行失败".format(cmd_name))
            raise Exception("{}运行失败".format(cmd_name))

    def check_exists(self, file_path):
        """
        用于检查文件及文件夹是否存在
        """
        if not os.path.exists(file_path):
            raise Exception("文件或文件夹{}不存在！".format(file_path))

    def run(self):
        """
        运行
        :return:
        """
        super(RegionAnnoTool, self).run()
        # self.get_select_region()
        self.run_region_gene()
        self.run_eff_detail("kegg")
        self.run_eff_detail("go")
        self.run_eff_detail("eggnog")
        if os.path.exists(self.output_dir + "/select_region.kegg.stat.detail"):
            self.logger.info('kegg')
            # 若结果文件产生了，就证明kegg.stat文件不为空且存在。就运行picture
            self.run_get_picture()
        self.set_db()
        # self.end()

    def set_db(self):
        """
        更新主表的chr_list和path:marker.filtered.marker
        self.filtered_marker_path
        """
        region_id = ObjectId(self.option("main_id"))
        genome_version_id = ObjectId(self.option('genome_version_id'))
        self.logger.info("开始region_anno的导表！")
        api = self.api.api("dna_evolution.region_anno")
        task_id = self.option('task_id')
        # api.get_update_chrlist(marker_id, self.target_dir_ + "/pop.filtered.marker",
        #                        self.target_dir_ + "/pop.filtered.detail.info")
        api.add_sg_region_anno_detail(region_id, self.option('pop_summary_path').prop['path'], genome_version_id)
        api.sg_region_anno_go_stat(task_id, region_id, self.output_dir + "/select_region.go.stat.detail",
                                   self.option('go_summary_path').prop['path'])
        api.sg_region_anno_kegg_stat(task_id, region_id, self.output_dir + "/select_region.kegg.stat.detail",
                                     self.output_dir + "/pathway_dir",
                                     os.path.join(self.option("pathway_path"), "region_anno/pathway_dir"))
        api.sg_region_anno_eggnog_stat(task_id, region_id, self.output_dir + "/select_region.eggnog.stat.detail")
        self.logger.info("设置region_anno的导表成功！")
        self.end()
