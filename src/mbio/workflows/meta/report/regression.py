# -*- coding: utf-8 -*-
# __author__ = 'gaohao guanqing' 

from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
from biocluster.core.exceptions import OptionError
from mbio.packages.beta_diversity.filter_newick import get_level_newicktree
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class RegressionWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RegressionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "group_table", "type": "infile","format": "meta.otu.group_table"},
            #{"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.otu_table"}, #guanqing.zou 20180517 
            {"name": "env_labs", "type": "string", "default": ""},
            #{"name": "PCAlabs", "type": "string", "default": ""},
            {"name": "environmental_regression_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "otu_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            {"name": "env_id","type":"string"},
            {"name": "diversity_type", "type": "string", "default": "beta"},
            {"name": "diversity_analysis_type", "type": "string"},
            {"name": "distance_type", "type": "string","default":""}
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.regression = self.add_tool('meta.regression')  #guanqing.zou 20180521 和宏基因组分开
        self.regression_calculation = self.add_module('meta.regression_calculation')  #guanqing.zou 20180517
        self.sort_tax_samples = self.add_tool("meta.otu.sort_samples_mg")
        #self.get_tax_abund_table = self.add_tool('meta.create_abund_table')
        #self.regression_list = [self.sort_tax_samples,self.sort_func_samples]

    def run_get_tax_abund_table(self):
        self.logger.info("已经运行啦")
        self.get_tax_abund_table.set_options({
        'anno_table': self.option('tax_anno_table'),
        'geneset_table': self.option('geneset_table'),
        'level_type': self.option('tax_level'),
        })
        self.get_tax_abund_table.on("end", self.run_tax_sort_samples)
        self.get_tax_abund_table.run()


    def run_tax_sort_samples(self):
        abund_table = self.option("otu_table")
        self.sort_tax_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group_table'),
        })
        self.sort_tax_samples.run()


    def run_regression_calculation(self):
        """
        获得排序分析的坐标beta多样性分析和alpha多样性分析
        fix by qingchen.zhang @20200910 增加计算unifrac距离等方法
        :return:
        """
        tax_abund_table = self.sort_tax_samples.option("out_otu_table").prop['path']
        options = ({
            "tax_abund_file": tax_abund_table,
            "diversity_analysis_type": self.option('diversity_analysis_type'),
            "diversity_type": self.option('diversity_type'),
            "distance_type": self.option('distance_type')
            })
        if 'unifrac' in self.option('distance_type'):  # fix by qingchen.zhang @20200910加unifrac距离
            if self.option('level') != 9:
                newicktree = get_level_newicktree(self.option('otu_id'), level=self.option('level'),tempdir=self.work_dir, return_file=False, bind_obj=self)
                all_find = re.findall(r'\'.+?\'', newicktree)
                for n, m in enumerate(all_find):
                    all_find[n] = m.strip('\'')
                all_find = dict((i[1], i[0]) for i in enumerate(all_find))

                def match_newname(matchname):
                    if hasattr(match_newname, 'count'):
                        match_newname.count = match_newname.count + 1
                    else:
                        match_newname.count = 1
                    return 'OTU' + str(match_newname.count)
                newline = re.sub(r'\'.+?\'', match_newname, newicktree)
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newline)
                tempfile.close()
                self.logger.info('get_newick:' + temp_tree_file)
                otu_table = tax_abund_table
                temp_otu_file = tax_abund_table + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                if len(all_lines) < 3:
                    self.logger.error('分类水平：%s,OTU表数据少于2行：%s' % (self.option('level'), len(all_lines)))
                    self.set_error("OTU表数据少于2行")
                self.logger.info(len(all_lines))
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:
                    name = line.split('\t')
                    origin_name = name[0].split("; ")[-1].strip()
                    if name[0] in all_find:
                        name[0] = 'OTU' + str(all_find[name[0]] + 1)
                    new_all.append('\t'.join(name))
                    if name[0] not in self.rename:
                        self.rename[name[0]] = origin_name
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                os.rename(otu_table, os.path.join(self.work_dir, "otu_table2.xls"))
                os.rename(temp_otu_file, otu_table)
                options['tax_abund_file'] = otu_table
                options['phy_newick'] = temp_tree_file
            else:
                newicktree = get_level_newicktree(self.option('otu_id'), level=self.option('level'),
                                                  tempdir=self.work_dir, return_file=False, bind_obj=self)
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newicktree)
                tempfile.close()
                otu_table = tax_abund_table
                temp_otu_file = tax_abund_table + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:  # OTU表中有复杂的名称OTU名称，包含进化物种类型，进化树种只有OTU名称
                    name = line.split('\t')
                    origin_name = name[0].split("; ")[-1].strip()
                    name[0] = name[0].split(';')[-1].strip()
                    new_all.append('\t'.join(name))
                    if name[0] not in self.rename:
                        self.rename[name[0]] = origin_name
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                os.rename(otu_table, os.path.join(self.work_dir, "otu_table2.xls"))
                os.rename(temp_otu_file, otu_table)
                options['tax_abund_file'] = otu_table
                options['phy_newick'] = temp_tree_file
        self.regression_calculation.set_options(options)
        self.regression_calculation.on("end", self.run_regression)#
        self.regression_calculation.run()

    def replace_name(self, input):
        """
        对文件替换为原来的名称
        :return:
        """
        out_table = os.path.join(self.work_dir, "out_table.xls")
        with open(input, "r") as f, open(out_table, "w") as w:
            lines = f.readlines()
            w.write(lines[0])
            for line in lines[1:]:
                line = line.strip().split("\t")
                sp_name = line[0].split("; ")[-1].strip()
                if sp_name in self.rename:
                    line[0] = self.rename[sp_name]
                w.write("\t".join(line) + "\n")
        os.remove(input)
        os.rename(out_table, input)

    def get_new_file(self,output_old,pc,output_new):
        if pc in ['PC1','MDS1','pc1','mds1']:
            index = 1
        elif pc in ['PC2','MDS2','pc2','mds2']:
            index = 2
        elif pc in ['PC3','pc3']:
            index = 3
        else:
            index = 1
        fr = open(output_old)
        fw = open(output_new,'w')
        lines = fr.readlines()
        for line in lines:
            spline = line.split('\t')
            fw.write("\t".join([spline[0],spline[index],'\n']))
        fw.close()
        fr.close()
    
        

    def run_regression(self):

        ###guanqing 20180517 修改使画图的x轴为PC值
        #tax_abund_table = self.regression_calculation.option("output_tax").prop['path']
        tax_abund_table = self.option("envtable")    
        #func_abund_table = self.option("envtable")
        func_abund_table =  self.regression_calculation.option("output_tax").prop['path']
        #func_abund_table_new = '/'.join(func_abund_table.split('/')[:-1])+'/pick.xls'
        #self.get_new_file(func_abund_table,self.option('PCAlabs'),func_abund_table_new)  #guanqing.zou 20180518 弃用函数get_new_file

        # by houshuang 20191009 alpha多样性分析只有一个指数
        if self.option('diversity_type') == "alpha":
            pcnum = 1
        else:
            pcnum = 3
        self.regression.set_options({
            "taxon_table": tax_abund_table,
            "func_table": func_abund_table,
            "pcnum": pcnum
            })
        self.regression.on("end", self.set_db)
        self.regression.run()


    def run(self):
        #self.run_get_tax_abund_table()
        #self.on_rely(self.regression_list, self.run_regression_calculation)
        self.sort_tax_samples.on("end",self.run_regression_calculation)
        self.run_tax_sort_samples()

        #self.run_regression_calculation()  #guanqing
        super(RegressionWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        self.api_regression = self.api.api("regression")
        data_path = self.regression.output_dir+"/Regression.data.xls"
        line_path = self.regression.output_dir+"/Regression.message.xls"
        if self.option('diversity_analysis_type') in ['nmds','NMDS']:
            mds_pc = 'MDS'
        elif self.option('diversity_type') == "alpha":  # by houshuang 20191009
            mds_pc = self.option('diversity_analysis_type').capitalize()
        else:
            mds_pc = 'PC'
        self.api_regression.add_resgression_site_detail(file_path=data_path,table_id=self.option("environmental_regression_id"),MDS_PC=mds_pc)
        self.api_regression.add_resgression_message_detail(file_path=line_path,table_id=self.option("environmental_regression_id"),MDS_PC=mds_pc)
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("environmental_regression_id"), "sg_environmental_regression")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("environmental_regression_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "rank_regression",
                "interaction": 1,
                "main_table": "sg_environmental_regression",
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.regression.output_dir))
        save_params(self.regression.output_dir, self.id)
        if self.option("diversity_type") == "alpha" or self.option("distance_type") == "pca":
            if os.path.exists(self.regression.output_dir + "/运行参数.txt"):
                os.rename(self.regression.output_dir + "/运行参数.txt", self.regression.output_dir + "/运行参数1.txt")
                with open(self.regression.output_dir + "/运行参数1.txt") as f, open(self.regression.output_dir + "/运行参数.txt", "w") as t:
                    data = f.readlines()
                    for i in data:
                        if "bray_curtis" in i:
                            pass
                        else:
                            t.write(i)
                os.remove(self.regression.output_dir + "/运行参数1.txt")
        result_dir = self.add_upload_dir(self.regression.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "回归分析结果目录", 0, "110166"],
            ["./Regression.data.xls", "xls", "散点图数据", 0, "110167"],
            ["./Regression.message.xls", "xls", "回归线和R^2的数据", 0, "110168"],
            ["./排序回归分析结果图.pdf", "pdf", "环境因子排序回归分析结果图", 0, ""]
        ])
        super(RegressionWorkflow, self).end()
