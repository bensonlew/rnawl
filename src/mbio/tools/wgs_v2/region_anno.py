# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190326

import os
import json
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class RegionAnnoAgent(Agent):
    """
    功能注释筛选及统计
    """
    def __init__(self, parent):
        super(RegionAnnoAgent, self).__init__(parent)
        options = [
            {"name": "pop_summary", "type": "infile", "format": "dna_evolution.pop_summary", "required": True},
            {"name": "go_level2", "type": "infile", "format": "wgs_v2.vcf", "required": True},  # GO id第二层级文件
            {"name": "region_select", "type": "string", "required": True},  # 区域
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "species_version_id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("pop_summary").is_set:
            raise OptionError("请设置pop.summary")
        if not self.option("region_select"):
            raise OptionError("请设置进行筛选的区域")
        if self.option("region_select") != "all":
            try:
                region_select = json.loads(self.option("region_select"))
            except:
                raise OptionError("参数region_select的格式不正确，请检查")
        if not self.option("go_level2").is_set:
            raise OptionError("请设置GO第二层级文件，用于导表")

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(RegionAnnoAgent, self).end()


class RegionAnnoTool(Tool):
    def __init__(self, config):
        super(RegionAnnoTool, self).__init__(config)
        self.python_path = "program/Python/bin/python"
        self.region_anno_path = self.config.PACKAGE_DIR + "/wgs_v2/region_anno.py"
        self.ko_path = self.config.SOFTWARE_DIR + "/database/KEGG/kegg_2017-05-01/kegg/pathway/ko/"
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.parafly = "program/parafly-r2013-01-21/src/ParaFly"

    def run_region_anno(self):
        """
        进行注释的筛选
        """
        self.logger.info(type(self.option("region_select")))
        cmd = "{} {} -i {} ".format(self.python_path, self.region_anno_path, self.option("pop_summary").prop["path"])
        cmd += "-r '{}' -o {}".format(self.option("region_select"), self.output_dir)
        command = self.add_command("region_anno", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("region_anno运行完成")
        else:
            self.set_error("region_anno运行失败，请检查")

    def run_get_picture(self):
        """
        根据ko号找到对应的所有ko的pdf和png文件-- 添加如果存在pathway_dir就将该文件夹删除
        """
        pathway_dir = os.path.join(self.output_dir, "pathway_dir")
        if os.path.isdir(pathway_dir):
            os.system('rm -r %s' % pathway_dir)
        os.system("mkdir {}".format(pathway_dir))
        self.num = 0
        cmd_list = []
        all_ko = os.listdir(self.ko_path)
        with open(self.output_dir + "/region.kegg.stat")as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split('\t')
                ko_png = item[0] + '.png'
                if ko_png in all_ko:
                    ko_path = os.path.join(self.ko_path, ko_png)
                    os.link(ko_path, self.output_dir + '/pathway_dir/' + ko_png)
                    pdf_path = os.path.join(self.output_dir + '/pathway_dir', item[0] + '.pdf')
                    cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + ko_path + ' ' + pdf_path
                    cmd_list.append(cmd)
                else:
                    self.set_error("ko{} png文件在本地数据库{}中没有，请核实!".format(item[0], self.ko_path))
        with open(self.work_dir + "/kegg_png_pdf_cmd.list", "w") as fw:
            for i in range(len(cmd_list)):
                fw.write(cmd_list[i] + "\n")
        cmd = self.parafly + " -c {} -CPU 10".format(self.work_dir + "/kegg_png_pdf_cmd.list")
        command = self.add_command("cmd_list", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("kegg的png转pdf运行成功")
        else:
            self.set_error("kegg的png转pdf运行出错!")

    def set_db(self):
        self.logger.info("开始导表")
        anno_api = self.api.api("wgs_v2.region_anno")
        anno_id = self.option("main_id")
        data_json = os.path.dirname(self.work_dir) + "/data.json"
        s3_output_dir = json.loads(open(data_json).read())["output"]
        species_version_id = self.option("species_version_id")
        pop_summary_path = os.path.join(self.output_dir, "region.summary")
        go_summary_path = self.option("go_level2").prop["path"]
        go_stat_detail = os.path.join(self.output_dir, "region.go.stat")
        kegg_stat_detail = os.path.join(self.output_dir, "region.kegg.stat")
        pathway_dir = os.path.join(s3_output_dir, "pathway_dir")
        origin_path = os.path.join(self.output_dir, "pathway_dir")
        eggnog_stat_detail = os.path.join(self.output_dir, "region.go.stat")
        pfam_stat_detail = os.path.join(self.output_dir, "region.pfam.stat")
        anno_api.add_sg_region_anno_detail(anno_id, pop_summary_path, species_version_id)
        anno_api.sg_region_anno_go_stat(anno_id, go_stat_detail, go_summary_path)
        anno_api.sg_region_anno_kegg_stat(anno_id, kegg_stat_detail, pathway_dir, origin_path)
        anno_api.sg_region_anno_eggnog_stat(anno_id, eggnog_stat_detail)
        anno_api.sg_region_anno_pfam_stat(anno_id, pfam_stat_detail)
        pop_summary = os.path.join(s3_output_dir, "region.summary")
        anno_api.update_info(coll="sg_region_anno", main_id=anno_id, update_dict={"pop_summary": pop_summary})
        self.logger.info("注释导表完成")

    def run(self):
        super(RegionAnnoTool, self).run()
        self.run_region_anno()
        self.run_get_picture()
        if self.option("main_id"):
            self.set_db()
        self.end()
