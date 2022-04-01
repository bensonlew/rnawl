# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.workflow import Workflow
import os
import re
import glob
from bson import ObjectId
from mbio.packages.beta_diversity.filter_newick import get_level_newicktree
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
from mbio.packages.meta.save_params import save_params


class MantelTestWorkflow(Workflow):
    """
    报告计算mantel检验
    last_modify by zhujuan 2018.02.02 因mantel检验的环境因子矩阵计算tool换成了distance_calc.py file类型检查跟着修改成"meta.otu.otu_table"
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(MantelTestWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},  # 输入的OTU id
            {"name": "env_file", "type": "infile", 'format': "meta.otu.otu_table"},  # 输入的环境因子表 id
            {"name": "newicktree", "type": "infile", 'format': "meta.otu.group_table"},  # 输入的环境因子表 id
            {"name": "level", "type": "int"},
            {"name": "otu_id", "type": "string"},
            {"name": "env_id", "type": "string"},
            {"name": "mantel_id", "type": "string"},
            {"name": "env_labs", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "units", "type": "string"},  # partial factor
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "env_method", "type": "string"},
            {"name": "otu_method", "type": "string"},
            {"name": "partialmatrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"},
            {"name": "dis_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"},
            {"name": "fac_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"}
        ]
        self.add_option(options)
        # print(self._sheet.options())
        self.set_options(self._sheet.options())
        self.mantel = self.add_module('statistical.mantel_test')
        self.params = {}

    def run(self):
        options = self.get_options()
        self.mantel.set_options(options)
        self.mantel.on('end', self.set_db)
        self.mantel.run()
        self.output_dir = self.mantel.output_dir
        super(MantelTestWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_mantel = self.api.meta_species_env
        mantel_result = glob.glob(self.output_dir + "/Discompare/*")[0]
        if self.option('units'):
            partial_matrix = glob.glob(self.output_dir + "/partial/*")[0]
            dis_matrix = glob.glob(self.output_dir + "/Otudistance/*")[0]
            fac_matrix = glob.glob(self.output_dir + "/Facdistance/*")[0]
            # name = "partial_mantel_test" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        else:
            dis_matrix = glob.glob(self.output_dir + "/Otudistance/*")[0]
            fac_matrix = glob.glob(self.output_dir + "/Facdistance/*")[0]

        if not os.path.isfile(mantel_result):
            self.logger.error("找不到报告文件:{}".format(mantel_result))
            self.set_error("找不到报告文件", code="12702201")
        # mantel_id = api_mantel.add_mantel_table(self.option('level'), self.option('otu_id'), self.option('env_id'), name=name, params=self.params)
        mantel_id = ObjectId(self.option("mantel_id"))
        api_mantel.add_mantel_detail(mantel_result, mantel_id)
        if self.option('units'):
            api_mantel.add_mantel_matrix(partial_matrix, "partial_matrix", mantel_id)
        api_mantel.add_mantel_matrix(dis_matrix, "species_matrix", mantel_id)
        api_mantel.add_mantel_matrix(fac_matrix, "env_matrix", mantel_id)
        # self.add_return_mongo_id('sg_species_mantel_check', mantel_id)
        self.end()

    def end(self):
        if os.path.exists(self.output_dir+"/Otudistance"):
            os.rename(self.output_dir+"/Otudistance", self.output_dir+"/Distance")
        if os.path.exists(self.output_dir+"/Distance/"+ self.option("otu_method") + "_otu_file.xls"):
            os.rename(self.output_dir+"/Distance/"+ self.option("otu_method") + "_otu_file.xls", self.output_dir+"/Distance/" + self.option("otu_method") + "_file.xls")
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "MantelTest分析结果目录", 0, "110157"],  # add 8 lines by hongdongxuan 20170324
            ["./Discompare", "", "Mantel_Test分析结果目录", 0, "110160"],
            ["./Discompare/partial_mantel_results.txt", "txt", "Mantel_Test分析结果表", 0, "110161"],
            ["./Facdistance", "", "环境因子矩阵结果目录", 0, "110158"],
            ["./Facdistance/factor_out.xls", "xls", "环境因子矩阵结果表", 0, "110159"],
            ["./Distance", "", "群落矩阵结果目录", 0, "110162"],
            ["./Distance/%s_file.xls"%self.option("otu_method"), "xls", "群落矩阵结果表", 0, "110163"],
            ["./partial", "", "限制环境因子矩阵结果目录", 0, "110164"],
            ["./partial/factor_out.xls", "xls", "限制环境因子矩阵结果表", 0, "110165"]
            # [".//mantel_results.txt", "txt", "mantel检验结果"]
        ])
        # print self.get_upload_files()
        super(MantelTestWorkflow, self).end()

    def get_options(self):
        options = {
            'otutable': self.option('otu_file'),
            'factor': self.option('env_file'),
            # 'factorselected': self.option('env_labs'),
            # 'partial_factor': self.option('units'),
            'otumatrixtype': self.option('otu_method'),
            'factormatrixtype': self.option('env_method')
        }
        if self.option('units'):
            options['partial_factor'] = self.option('units')
        if 'unifrac' in self.option('otu_method'):  # sanger_bioinfo/src/mbio/workflows/meta/report/distance_calc.py中的解释
            if self.option('level') != 9:
                newicktree = get_level_newicktree(self.option('otu_id'), level=self.option('level'),
                                                  tempdir=self.work_dir, return_file=False, bind_obj=self)
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
                otu_table = self.option('otu_file').path
                temp_otu_file = self.option('otu_file').path + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                if len(all_lines) < 3:
                    self.logger.error('分类水平：%s,otu表数据少于2行：%s' % (self.option('level'), len(all_lines)))
                    self.set_error("otu表数据少于2行", code="12702202")
                self.logger.info(len(all_lines))
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:
                    name = line.split('\t')
                    if name[0] in all_find:
                        name[0] = 'OTU' + str(all_find[name[0]] + 1)
                    new_all.append('\t'.join(name))
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                options['otutable'] = temp_otu_file
                options['newicktree'] = temp_tree_file
            else:
                newicktree = get_level_newicktree(self.option('otu_id'), level=self.option('level'),
                                                  tempdir=self.work_dir, return_file=False, bind_obj=self)
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newicktree)
                tempfile.close()
                otu_table = self.option('otu_file').path
                temp_otu_file = self.option('otu_file').path + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:  # OTU表中有复杂的名称OTU名称，包含进化物种类型，进化树种只有OTU名称
                    name = line.split('\t')
                    name[0] = name[0].split(';')[-1].strip()
                    new_all.append('\t'.join(name))
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                options['otutable'] = temp_otu_file
                options['newicktree'] = temp_tree_file
        # print(options)
        return options
