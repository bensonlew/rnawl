# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20180703

import os
import re
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class QtlAnalysisModule(Module):
    """
    qtl定位及关联区域注释
    """
    def __init__(self, work_id):
        super(QtlAnalysisModule, self).__init__(work_id)
        options = [
            {"name": "trait_dir", "type": "infile", "format": "dna_gmap.trait_dir"},  # 性状文件夹
            {"name": "total_map", "type": "infile", "format": "dna_gmap.qtl_map_dir"},  # pop_type为F2时是total.csv，CP的时候是文件夹
            {"name": "total_phase", "type": "infile", "format": "dna_gmap.qtl_map_dir"},  # pop_type为CP时才有
            {"name": "pop_final_vcf", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "pop_summary", "type": "infile", "format": "dna_gmap.lg"},
            {"name": "pop_type", "type": "string", "default": "F2"},  # pop type, CP/F2
            {"name": "location_method", "type": "string", "default": "cim"},  # 定位方法,cim/scanone
            {"name": "pm_num", "type": "int", "default": 1000},  # 置换检验次数，[0,+∞]
            {"name": "p_value", "type": "float"},  # 阈值，p value,默认值0.05
            {"name": "lod_value", "type": "int"},  # 阈值，LOD,默认值3
            {"name": "pid", "type": "string"},  # 父本
            {"name": "mid", "type": "string"},  # 母本
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.end_times = 0
        self.new_trit = []

    def check_options(self):
        if not self.option("trait_dir").is_set:
            raise OptionError("请设置性状文件夹", code="24800501")
        if not self.option("total_map").is_set:
            raise OptionError("请设置map的结果文件", code="24800502")
        if self.option("pop_type") not in ["F1", "F2", "DH", "BC", "Ril"]:
            raise OptionError("pop_type类型只能为F1/F2/DH/BC/Ril,而不能是%s", variables=(self.option("pop_type")), code="24800503"
)
        if self.option("pop_type") == "F1":
            if not self.option("total_phase").is_set:
                raise OptionError("请设置map对应的phase", code="24800504")
        if self.option("location_method") not in ["cim", "scanone"]:
            raise OptionError("定位方法出错，只能为cim/scanone", code="24800505")
        if self.option("pm_num") < 0:
            raise OptionError("置换检验次数只能是大于等于0的正整数", code="24800506")
        if not self.option("p_value") and not self.option("lod_value"):
            raise OptionError("必须设置阈值P_value或LOD", code="24800507")
        if not self.option("pid") or not self.option("mid"):
            raise OptionError("请设置亲本和母本", code="24800508")
        if self.option("pop_type") == "F1" and self.option("location_method") == "cim":
            raise OptionError("群体类型为F1的时候定位方法不能是cim,请检查", code="24800509")

    def run_qtl(self):
        self.qtl_list = []
        self.trait_file = self.option("trait_dir").prop["trait_file"]
        self.trait_list = self.trait_file.keys()
        options = {
            "pop_type": self.option("pop_type"),
            "location_method": self.option("location_method"),
            "pm_num": self.option("pm_num"),
            "p_value": self.option("p_value"),
            "lod_value": self.option("lod_value")
        }
        if self.option("pop_type") == "F1":
            for sex in ["female", "male", "sexAver"]:
                options["total_map"] = os.path.join(self.option("total_map").prop["path"], "total." + sex + ".map")
                # options["total_phase"] = self.option("total_phase").prop["path"] + "/total." + sex + ".loc"
                options["total_phase"] = os.path.join(self.option("total_phase").prop["path"], "total.sexAver.loc")
                for trait in self.trait_list:
                    options["trait_file"] = self.trait_file[trait]
                    qtl = self.add_tool("dna_gmap.qtl")
                    qtl.set_options(options)
                    qtl.on("end", self.run_regin_anno, trait)
                    qtl.on("end", self.set_output, "qtl")
                    self.qtl_list.append(qtl)
        else:
            options["total_map"] = os.path.join(self.option("total_map").prop["path"], "total.csv")
            for trait in self.trait_list:
                options["trait_file"] = self.trait_file[trait]
                qtl = self.add_tool("dna_gmap.qtl")
                qtl.set_options(options)
                qtl.on("end", self.run_regin_anno, trait)
                qtl.on("end", self.set_output, "qtl")
                self.qtl_list.append(qtl)
        for tool in self.qtl_list:
            tool.run()

    def run_regin_anno(self, event):
        obj = event["bind_object"]
        qtl_csv = ""
        for f in os.listdir(obj.output_dir):
            if f.endswith(".qtl.csv"):
                self.new_trit.append(event["data"])
                qtl_csv = os.path.join(obj.output_dir, f)
        if qtl_csv:
            options = {
                "qtl_csv": qtl_csv,
                "pop_final_vcf": self.option("pop_final_vcf"),
                "pop_summary": self.option("pop_summary")
            }
            region_anno = self.add_tool("dna_gmap.region_anno")
            region_anno.set_options(options)
            region_anno.on("end", self.run_anno_analysis, event["data"])
            region_anno.on("end", self.run_region_info, event["data"])
            region_anno.on("end", self.set_output, "region_anno")
            region_anno.run()
        else:
            self.end_times += 2
            self.set_output("end")

    def run_anno_analysis(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            if f.endswith(".kegg.stat"):
                kegg_stat = os.path.join(obj.output_dir, f)
        options = {
            "anno_stat": kegg_stat,
            "anno_type": "kegg"
        }
        anno_analysis = self.add_tool("bsa.anno_analysis")
        anno_analysis.set_options(options)
        anno_analysis.on("end", self.set_output, "kegg_anno")
        anno_analysis.run()

    def run_region_info(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            if f.endswith(".vcf.total"):
                vcf_total = os.path.join(obj.output_dir, f)
        options = {
            "vcf_total": vcf_total,
            "trait_file": self.trait_file[event["data"]],
            "pop_type": self.option("pop_type"),
            "pid": self.option("pid"),
            "mid": self.option("mid")
        }
        region_info = self.add_tool("dna_gmap.region_info")
        region_info.set_options(options)
        region_info.on("end", self.set_output, "region_info")
        region_info.run()

    def get_total_qtl_result(self):
        """
        统计所有性状的qtl.result文件，进行汇总
        """
        with open(self.output_dir + "/qtl_result.xls", "w") as w:
            w.write("Trait\tChr\tPosition(cM)\tLOD\tR2(%)\tStart(cM)\tEnd(cM)\tMarker ID\tMarker start\tMarker end\n")
            for trait in self.new_trit:
                qtl_csv = os.path.join(self.output_dir, trait + ".qtl.result")
                with open(qtl_csv, "r") as f:
                    lines = f.readlines()
                    for line in lines[1:]:
                        item = line.strip().split("\t")
                        w.write(trait + "\t" + item[0].split("_")[0] + "\t" + item[2] + "\t" + item[3] + "\t")
                        w.write(item[4] + "\t" + item[7] + "\t" + item[8] + "\t" + item[0] + "\t" + item[9])
                        w.write("\t" + item[10] + "\n")

    def set_output(self, event):
        try:
            obj = event["bind_object"]
            if event["data"] == "qtl":
                for f in os.listdir(obj.output_dir):
                    if f.endswith(".scan.csv"):
                        self.link(obj.output_dir, f)
            elif event["data"] == "region_anno":
                for f in os.listdir(obj.output_dir):
                    if re.match(r".*.(?:eggnog|go.stat|total)", f):
                        self.link(obj.output_dir, f)
                    if f.endswith(".qtl.result"):
                        self.link(obj.output_dir, f)
            else:
                for f in os.listdir(obj.output_dir):
                    self.link(obj.output_dir, f)
        except:
            self.logger.info("qtl没有定位到结果")
        self.end_times += 1
        if self.end_times == len(self.qtl_list) * 4:
            self.get_total_qtl_result()
            if self.option("main_id"):
                self.set_db()
            self.end()

    def link(self, old_dir, f):
        old = os.path.join(old_dir, f)
        new = os.path.join(self.output_dir, f)
        if f == "pathway_dir":
            if not os.path.exists(new):
                os.mkdir(new)
            for f1 in os.listdir(old):
                old1 = os.path.join(old, f1)
                new1 = os.path.join(new, f1)
                if not os.path.exists(new1):
                    os.link(old1, new1)
        else:
            if os.path.exists(new):
                os.remove(new)
            os.link(old, new)

    def set_db(self):
        """
        """
        self.logger.info("将结果保存到mongo")
        qtl_api = self.api.api("dna_gmap.qtl")
        qtl_id = self.option("main_id")
        if self.option("pop_type") == "F1":
            total_map = os.path.join(self.option("total_map").prop["path"], "total.sexAver.map")
            type = "map"
            parent_source = ["female", "male", "sexvar"]
        else:
            total_map = os.path.join(self.option("total_map").prop["path"], "total.map")
            type = "csv"
            parent_source = ["total"]
            data_type = "total"
        qtl_api.update_sg_qtl(qtl_id=qtl_id, parent_source=parent_source)
        for data_type in parent_source:
            trait_list = qtl_api.add_sg_qtl_result(qtl_id=qtl_id, qtl_path=self.output_dir + "/qtl_result.xls", data_type=data_type)
            qtl_api.add_sg_qtl_region(qtl_id=qtl_id, qtl_dir=self.output_dir, data_type=data_type)
            qtl_api.add_sg_box_region(qtl_id=qtl_id, qtl_dir=self.output_dir, data_type=data_type)
            qtl_api.add_sg_manhattan(qtl_id=qtl_id, trait_list=trait_list, qtl_dir=self.output_dir, data_type=data_type)
            qtl_api.add_sg_distribution(qtl_id=qtl_id, total_map=total_map, qtl_dir=self.output_dir, type=type, data_type=data_type)
            anno_api = self.api.api("dna_gmap.region_anno")
            for trit in trait_list:
                trait = trit if data_type == "total" else trit + "_" + data_type
                anno_id = anno_api.add_sg_region_anno(qtl_id, trait)
                gene_total_path = os.path.join(self.output_dir, trit + ".gene.total")
                go_path = os.path.join(self.output_dir, trit + ".gene.go.stat")
                kegg_path = os.path.join(self.output_dir, trit + ".gene.kegg.final.stat.detail")
                eggnog_path = os.path.join(self.output_dir, trit + ".gene.eggnog.stat")
                # pathway_dir = os.path.join(self.output_dir, "pathway_dir")
                pathway_dir = os.path.join(self._sheet.output, "pathway_dir/")
                anno_api.sg_region_anno_detail(anno_id, gene_total_path)
                anno_api.sg_region_anno_go_stat(anno_id, go_path)
                anno_api.sg_region_anno_kegg_stat(anno_id, kegg_path, pathway_dir)
                anno_api.sg_region_anno_eggnog_stat(anno_id, eggnog_path)

    def run(self):
        super(QtlAnalysisModule, self).run()
        self.run_qtl()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(QtlAnalysisModule, self).end()
