# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.1008


from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
import os, re
import math
import json
import xlrd
import time, datetime
from bson.objectid import ObjectId


class GwasAnalysisModule(Module):
    """
    群体进化，GWAS关联分析接口的module
    """
    def __init__(self, work_id):
        super(GwasAnalysisModule, self).__init__(work_id)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "upload_trait_path", "type": "infile", "format": "dna_gmap.trait"},  # 上传的性状文件
            {"name": "vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "recode_vcf_path", "type": "outfile", "format": "dna_gmap.vcf"},
            {"name": "task_id", "type": "string"},
            # {"name": "recode", "type": "bool", "default": False},     # 参数不传，vcftools不接收
            # {"name": "remove_filtered_all", "type": "bool", "default": False},   # --remove-filtered-all
            # {"name": "remove_indels", "type": "bool", "default": False},     # --remove-indels
            {"name": "min_dp", "type": "int"},           # --minDP 2
            {"name": "max_dp", "type": "int"},           # --maxDP 6
            {"name": "max_missing", "type": "float", "default": "0.3"},   # --max-missing 必传
            {"name": "min_maf", "type": "float", "default": "0.05"},       # --maf 0.05
            {"name": "max_maf", "type": "float", "default": "1"},       # --max-maf 1
            {"name": "ref", "type": "string", "default": "ref"},  # 20180827代码新加 cui,可不传，pl代码默认ref
            {"name": "update_info", "type": "string"},
            {"name": "chrs_list", "type": "string"}
            # {"name": "update", "type": "string", "default": 0},
        ]
        self.add_option(options)
        self.trait_annovar = self.add_tool("dna_evolution.trait_annovar")
        self.vcftools_filter = self.add_tool("dna_evolution.vcftools_filter")
        self.vcf2hapmap = self.add_tool("dna_evolution.vcf2hapmap")
        self.vcf2trit = self.add_tool("dna_evolution.vcf2trit")
        self.gwas_mvp = self.add_module("dna_evolution.gwas_mvp")
        self.new_file_path = os.path.join(self.work_dir, "trait.xls")

    def check_options(self):
        # if not self.option("main_id"):
        #     raise OptionError("请设置main_id")
        if not self.option("upload_trait_path").prop['path']:
            raise OptionError("请上传upload_trait_path")
        if not self.option("vcf_path").is_set:
            raise OptionError("请设置vcf_path")

    def run_trait_annovar(self):
        """
        GWAS性状表格
        """
        self.trait_annovar.set_options({
            "upload_trait_path": self.new_file_path
            # "task_id": self.option("task_id"),
            # "main_id": self.option("main_id")
        })
        self.trait_annovar.on('end', self.set_output, "trait_annovar")
        self.trait_annovar.run()
        print("\n##### trait_annovar运行结束\n")

    def run_vcftools_filter(self):
        """
        /mnt/ilustre/users/dna/.env/bin/vcftools
        --recode
        --out /mnt/ilustre/users/minghao.zhang/newmdt/Project/MJ20180704017_liqiang/
            liqiangdata_20180713/Gwas_20180727/Result2/step01.vcf-filter/pop
        --remove-filtered-all --remove-indels
        --minDP 10  --max-missing 0.01
        --maf 0.05
        --vcf /mnt/ilustre/users/minghao.zhang/newmdt/Project/MJ20180704017_liqiang/
            liqiangdata_20180713/Gwas_20180727/Trit/trit2.sort.vcf

        """
        options = {
            "vcf_path": self.option("vcf_path").prop['path'],
            "recode": True,     # 未设置传参，直接写死True
            # "remove_indels": True,
            "remove_filtered_all": True,
            "max_missing": float(self.option('max_missing'))
        }
        if self.option('min_dp'):
            options["minDP"] = int(self.option('min_dp'))
        if self.option('max_dp'):
            options["maxDP"] = int(self.option('max_dp'))
        if self.option('min_maf'):
            options["min_maf"] = float(self.option('min_maf'))
        if self.option('max_maf'):
            options["max_maf"] = float(self.option('max_maf'))
        self.vcftools_filter.set_options(options)
        self.vcftools_filter.on('end', self.set_output, "vcftools_filter")
        self.vcftools_filter.on("end", self.run_vcf2hapmap)
        self.vcftools_filter.run()

    def run_vcf2hapmap(self):
        """
        pop.hapmap生成
        单倍体图谱
        """
        recode_vcf_path = os.path.join(self.vcftools_filter.output_dir, "pop.recode.vcf")
        self.vcf2hapmap.set_options({
            "recode_vcf_path": recode_vcf_path,
            "ref": self.option('ref')
        })
        self.vcf2hapmap.on("end", self.set_output, "vcf2hapmap")
        self.vcf2hapmap.run()

    def run_vcf2trit(self):
        """
        生成ref.chrlist
        """
        self.vcf2trit.set_options({
            "upload_trait_path": self.new_file_path,
        })
        self.vcf2trit.on("end", self.set_output, "vcf2trit")
        self.vcf2trit.run()

    def run_gwas_mvp(self):
        """
        MVP.single.R
        曼哈顿图
        """
        trit_list_path = os.path.join(os.path.join(self.output_dir, "vcf2trit"), "trit.list")
        pop_hapmap_path = os.path.join(os.path.join(self.output_dir, "vcf2hapmap"), "pop.hapmap")
        self.gwas_mvp.set_options({
            "trit_list_path": trit_list_path,
            "pop_hapmap_path": pop_hapmap_path
        })
        self.gwas_mvp.on("end", self.set_output, "gwas_mvp")
        self.gwas_mvp.run()

    def trait_xls(self):
        data = xlrd.open_workbook(self.option("upload_trait_path").prop["path"])
        sheet1 = data.sheet_by_index(0)
        write_lines = ""
        datas = sheet1.row_values(0)
        datas[0] = 'samples'
        write_lines = write_lines + "\t".join(datas) + "\n"
        for i in range(1, sheet1.nrows):
            datas = sheet1.row_values(i)
            new_list = []
            for x in datas:
                if str(x).upper() == 'NA' or x == '' or x == '-' or x == '/':
                    new_list.append("NA")
                else:
                    new_list.append(str(x))
            write_lines = write_lines + "\t".join(new_list) + "\n"
        with open(self.new_file_path, "w")as fw:
            fw.write(write_lines)

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'trait_annovar':
            self.linkdir(obj.output_dir, 'trait_annovar')
        elif event['data'] == 'vcftools_filter':
            self.linkdir(obj.output_dir, 'vcftools_filter')
        elif event['data'] == 'vcf2hapmap':
            self.linkdir(obj.output_dir, 'vcf2hapmap')
        elif event['data'] == 'vcf2trit':
            self.linkdir(obj.output_dir, 'vcf2trit')
        elif event['data'] == 'gwas_mvp':
            self.linkdir(obj.output_dir, 'gwas_mvp')
        else:
            pass

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
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

    def run(self):
        super(GwasAnalysisModule, self).run()
        self.trait_xls()
        self.on_rely([self.vcf2hapmap, self.vcf2trit], self.run_gwas_mvp)  # 跑完再执行gwas_mvp
        self.on_rely([self.trait_annovar, self.gwas_mvp], self.end)  # 跑完再end
        self.run_trait_annovar()
        self.run_vcftools_filter()
        self.run_vcf2trit()

    def end(self):
        self.option("recode_vcf_path").set_path(os.path.join(self.vcftools_filter.output_dir, "pop.recode.vcf"))
        seq_api = self.api.api("dna_evolution.gwas_analysis")
        gwas_id = self.option("main_id")
        bar_dir = os.path.join(self.trait_annovar.work_dir, "bar_dir")
        trait_stat = os.path.join(self.trait_annovar.output_dir, "annovar.trait.xls")
        seq_api.add_sg_gwas_analysis_stat(gwas_id=gwas_id, trait_stat=trait_stat, task_id=self.option("task_id"),
                                          bar_dir=bar_dir)
        #  曼哈顿图导表
        csv_dir = os.path.join(self.output_dir, "gwas_mvp")
        chrs_list = self.option("chrs_list")
        input_chr_list = chrs_list.strip().split(",")
        for f in os.listdir(csv_dir):
            chr_list = []
            pos_list = [0]
            farmcpu_path = os.path.join(os.path.join(csv_dir, f), ("MVP." + f + ".FarmCPU.csv"))
            glm_path = os.path.join(os.path.join(csv_dir, f), ("MVP." + f + ".GLM.csv"))
            mlm_path = os.path.join(os.path.join(csv_dir, f), ("MVP." + f + ".MLM.csv"))
            farmcpu_value_list = []
            mlm_value_list = []
            glm_value_list = []
            x_categories_list = []
            pvalue_default = 0
            with open(farmcpu_path, "r") as f1:
                lines = f1.readlines()
                last_chr = "origin"
                last_pos = 0
                farmcpu_value = []
                x_categories = []
                snp_num = len(lines)
                pvalue_default = -1 * math.log10(0.05/snp_num)
                for line in lines[1:]:
                    item = line.strip().split(",")
                    if item[0].strip("\"") in input_chr_list:
                        if item[4].strip("\"") == "NA" or item[4].strip("\"") == "NaN":
                            item_4 = 0
                        else:
                            item_4 = -1 * math.log10(float(item[4].strip("\"")))
                        if last_chr != "origin":
                            if last_chr != str(item[1].strip("\"")):
                                chr_list.append(str(item[1].strip("\"")).strip())
                                pos_list.append(last_pos)
                                x_categories_list.append(x_categories)
                                farmcpu_value_list.append(farmcpu_value)
                                farmcpu_value = []
                                x_categories = []
                                x_categories.append(float(item[2].strip("\"")))
                                farmcpu_value.append(item_4)
                            else:
                                x_categories.append(float(item[2].strip("\"")))
                                farmcpu_value.append(item_4)
                        else:
                            chr_list.append(str(item[1].strip("\"")).strip())
                            x_categories.append(float(item[2].strip("\"")))
                            farmcpu_value.append(item_4)
                        last_chr = str(item[1].strip("\""))
                        last_pos = float(item[2].strip("\""))
                pos_list.append(last_pos)
                x_categories_list.append(x_categories)
                farmcpu_value_list.append(farmcpu_value)
            json_list = []
            add_pos = 0
            print farmcpu_value_list
            for i in range(len(x_categories_list)):
                positon_list = []
                for m in range(len(x_categories_list[i])):
                    positon_list.append([(add_pos + x_categories_list[i][m]), farmcpu_value_list[i][m]])
                add_pos = add_pos + x_categories_list[i][-1]
                json_list.append(positon_list)
            s3_man_farmcpu = os.path.join(self.parent._sheet.output, ("gwas_mvp/" + f + "/Man." + f + ".FarmCPU.json"))
            man_farmcpu_json = os.path.join(os.path.join(csv_dir, f), ("Man." + f + ".FarmCPU.json"))
            with open(man_farmcpu_json, "w") as fw:
                # fw.write(str({"datas": json_list}))
                json.dump({"datas": json_list}, fw)
            seq_api.add_sg_manhattan(data_path=csv_dir, main_id=gwas_id, task_id=self.option("task_id"),
                                     ext_name="farmcpu", chr_list=chr_list, pos_list=pos_list,
                                     value_list=s3_man_farmcpu, trait=f, pvalue_default=pvalue_default)
            with open(glm_path, "r") as g1:
                lines = g1.readlines()
                last_chr = "origin"
                glm_value = []
                for line in lines[1:]:
                    item = line.strip().split(",")
                    if item[0].strip("\"") in input_chr_list:
                        if item[4].strip("\"") == "NA" or item[4].strip("\"") == "NaN":
                            item_4 = 0
                        else:
                            item_4 = -1 * math.log10(float(item[4].strip("\"")))
                        if last_chr != "origin":
                            if last_chr != str(item[1].strip("\"")):
                                glm_value_list.append(glm_value)
                                glm_value = []
                                glm_value.append(item_4)
                            else:
                                glm_value.append(item_4)
                        else:
                            glm_value.append(item_4)
                        last_chr = str(item[1].strip("\""))
                glm_value_list.append(glm_value)
            json_list = []
            add_pos = 0
            for i in range(len(x_categories_list)):
                positon_list = []
                for m in range(len(x_categories_list[i])):
                    positon_list.append([(add_pos + x_categories_list[i][m]), glm_value_list[i][m]])
                add_pos = add_pos + x_categories_list[i][-1]
                json_list.append(positon_list)
            s3_man_glm = os.path.join(self.parent._sheet.output, ("gwas_mvp/" + f + "/Man." + f + ".glm.json"))
            man_glm_json = os.path.join(os.path.join(csv_dir, f), ("Man." + f + ".glm.json"))
            with open(man_glm_json, "w") as fw:
                # fw.write(str({"datas": json_list}))
                json.dump({"datas": json_list}, fw)
            seq_api.add_sg_manhattan(data_path=csv_dir, main_id=gwas_id, task_id=self.option("task_id"),
                                     ext_name="glm", chr_list=chr_list, pos_list=pos_list, value_list=s3_man_glm,
                                     trait=f, pvalue_default=pvalue_default)
            with open(mlm_path, "r") as m1:
                lines = m1.readlines()
                last_chr = "origin"
                mlm_value = []
                for line in lines[1:]:
                    item = line.strip().split(",")
                    if item[0].strip("\"") in input_chr_list:
                        if item[4].strip("\"") == "NA" or item[4].strip("\"") == "NaN":
                            item_4 = 0
                        else:
                            item_4 = -1 * math.log10(float(item[4].strip("\"")))
                        if last_chr != "origin":
                            if last_chr != str(item[1].strip("\"")):
                                mlm_value_list.append(mlm_value)
                                mlm_value = []
                                mlm_value.append(item_4)
                            else:
                                mlm_value.append(item_4)
                        else:
                            mlm_value.append(item_4)
                        last_chr = str(item[1].strip("\""))
                mlm_value_list.append(mlm_value)
            json_list = []
            add_pos = 0
            for i in range(len(x_categories_list)):
                positon_list = []
                for m in range(len(x_categories_list[i])):
                    positon_list.append([(add_pos + x_categories_list[i][m]), mlm_value_list[i][m]])
                add_pos = add_pos + x_categories_list[i][-1]
                json_list.append(positon_list)
            s3_man_mlm = os.path.join(self.parent._sheet.output, ("gwas_mvp/" + f + "/Man." + f + ".mlm.json"))
            man_mlm_json = os.path.join(os.path.join(csv_dir, f), ("Man." + f + ".mlm.json"))
            with open(man_mlm_json, "w") as fw:
                # fw.write(str({"datas": json_list}))
                json.dump({"datas": json_list}, fw)
            seq_api.add_sg_manhattan(data_path=csv_dir, main_id=gwas_id, task_id=self.option("task_id"),
                                     ext_name="mlm", chr_list=chr_list, pos_list=pos_list, value_list=s3_man_mlm,
                                     trait=f, pvalue_default=pvalue_default)
        #  GWAS期望分布图导表
        for f in os.listdir(csv_dir):
            farmcpu_path = os.path.join(os.path.join(csv_dir, f), ("MVP." + f + ".FarmCPU.csv"))
            glm_path = os.path.join(os.path.join(csv_dir, f), ("MVP." + f + ".GLM.csv"))
            mlm_path = os.path.join(os.path.join(csv_dir, f), ("MVP." + f + ".MLM.csv"))
            farmcpu_list = []
            glm_list = []
            mlm_list = []
            farmcpu_value = []
            mlm_value = []
            glm_value = []
            with open(farmcpu_path, "r") as f1:
                lines = f1.readlines()
                for line in lines[1:]:
                    item = line.strip().split(",")
                    if item[4].strip("\"") != "NA" and item[4].strip("\"") != "NaN" and item[4].strip("\"") != 0:
                        farmcpu_list.append(float(item[4].strip("\"")))
                farmcpu_snp_num = len(farmcpu_list) - 1
                farmcpu_list = sorted(farmcpu_list)
            with open(glm_path, "r") as f1:
                lines = f1.readlines()
                for line in lines[1:]:
                    item = line.strip().split(",")
                    if item[4].strip("\"") != "NA" and item[4].strip("\"") != "NaN" and item[4].strip("\"") != 0:
                        glm_list.append(float(item[4].strip("\"")))
                glm_snp_num = len(glm_list) - 1
                glm_list = sorted(glm_list)
            with open(mlm_path, "r") as f1:
                lines = f1.readlines()
                for line in lines[1:]:
                    item = line.strip().split(",")
                    if item[4].strip("\"") != "NA" and item[4].strip("\"") != "NaN" and item[4].strip("\"") != 0:
                        mlm_list.append(float(item[4].strip("\"")))
                mlm_snp_num = len(mlm_list) - 1
                mlm_list = sorted(mlm_list)
            snp_num = max(farmcpu_snp_num, glm_snp_num, mlm_snp_num)
            x_max = -math.log10(float(1) / snp_num)
            y_max = -math.log10(min(farmcpu_list[0], glm_list[0], mlm_list[0]))
            for i in range(farmcpu_snp_num-1):
                    farmcpu_value.append([-math.log10(float(i+1)/farmcpu_snp_num), -math.log10(farmcpu_list[i])])
            farmcpu_value.append([0, -math.log10(farmcpu_list[-1])])
            farmcpu_value.reverse()
            farmcpu_json = os.path.join(os.path.join(csv_dir, f), ("MVP." + f + ".FarmCPU.json"))
            with open(farmcpu_json, "w") as fw:
                # fw.write(str(farmcpu_value))
                json.dump(farmcpu_value, fw)
            for i in range(mlm_snp_num-1):
                    mlm_value.append([-math.log10(float(i + 1) / mlm_snp_num), -math.log10(mlm_list[i])])
            mlm_value.append([0, -math.log10(mlm_list[-1])])
            mlm_value.reverse()
            mlm_json = os.path.join(os.path.join(csv_dir, f), ("MVP." + f + ".mlm.json"))
            with open(mlm_json, "w") as fw:
                # fw.write(str(mlm_value))
                json.dump(mlm_value, fw)
            for i in range(glm_snp_num-1):
                    glm_value.append([-math.log10(float(i + 1) / glm_snp_num), -math.log10(glm_list[i])])
            glm_value.append([0, -math.log10(glm_list[-1])])
            glm_value.reverse()
            glm_json = os.path.join(os.path.join(csv_dir, f), ("MVP." + f + ".glm.json"))
            with open(glm_json, "w") as fw:
                # fw.write(str(glm_value))
                json.dump(glm_value, fw)
            s3_farmcpu = os.path.join(self.parent._sheet.output, ("gwas_mvp/" + f + "/MVP." + f + ".FarmCPU.json"))
            s3_mlm = os.path.join(self.parent._sheet.output, ("gwas_mvp/" + f + "/MVP." + f + ".mlm.json"))
            s3_glm = os.path.join(self.parent._sheet.output, ("gwas_mvp/" + f + "/MVP." + f + ".glm.json"))
            seq_api.add_sg_scatter(task_id=self.option("task_id"), origin_id=gwas_id, gwas_output_dir=csv_dir,
                                   x_max=x_max, y_max=y_max, farmcpu_value=s3_farmcpu, mlm_value=s3_mlm,
                                   glm_value=s3_glm, trait=f)
        recode_vcf_path = os.path.join(os.path.join(self.output_dir, "vcftools_filter"), "pop.recode.vcf")
        trait_list = []
        with open(self.new_file_path, "r")as fr:
            line = fr.readline()
            header = line.strip().split('\t')
            for trait in header[1:]:
                trait_list.append(trait)
        self.api.api('dna_evolution.api_base').update_db_record(collection="sg_gwas_analysis",
                                                                query_dict={"_id": ObjectId(gwas_id)},
                                                                update_dict={"csv_dir": csv_dir,
                                                                             "recode_vcf_path": recode_vcf_path,
                                                                             "trait_list": trait_list})
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(GwasAnalysisModule, self).end()
