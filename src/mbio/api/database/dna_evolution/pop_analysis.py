# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180911

from api_base import ApiBase
from collections import defaultdict
import os
import re
import datetime
import json


class PopAnalysis(ApiBase):
    """
    群体进化，GWAS关联分析接口
    """
    def __init__(self, bind_object):
        super(PopAnalysis, self).__init__(bind_object)
        self._project_type = "dna_evolution"

    def add_sg_pop_pca_detail(self, pca_id, file_path):
        """
        用于导入PCA分析样本坐标表--pop.pca.eigenvec
        :param pca_id:
        :param file_path:
        :return:
        """
        data_list = []
        pca_id = self.check_objectid(pca_id)
        self.check_exists(file_path)
        pcas = ['PCA{}'.format(i) for i in range(1, len(open(file_path, 'rU').readlines()[0].strip().split(' ')) - 1)]
        with open(file_path, 'r') as r:
            datas = r.readlines()
            for line in datas:
                tmp = line.strip().split(' ')
                insert_data = {
                    "pca_id": pca_id,
                    "sample_name": tmp[0]
                }
                for i in range(2, len(tmp)):
                    insert_data.update({'PCA{}'.format(i - 1): round(float(tmp[i]), 3)})
                data_list.append(insert_data)
        if data_list:
            self.col_insert_data("sg_pop_pca_detail", data_list)
        else:
            self.bind_object.logger.info("文件{}为空不进行导表".format(file_path))
        self.update_db_record("sg_pop_pca", {"_id": pca_id}, {"pca_list": pcas})

    def add_sg_pop_pca(self, task_id, project_sn, params=None):
        """
        添加pca主表
        :param task_id:
        :param project_sn:
        :param params:
        :return:
        """
        if params:
            params_ = params
        else:
            params_ = {}
        main_id = self.add_main_table("sg_pop_pca", task_id, project_sn, params_, "origin_pca_analysis",
                                      "开始进行PCA分析")
        return main_id

    def add_sg_scatter(self, task_id, origin_id, file_path):
        """
        pca的散点图的导表，包含了主表与细节表-- pop.pca.eigenvec
        :param task_id:
        :param origin_id:
        :param file_path:
        :return:
        """
        origin_id = self.check_objectid(origin_id)
        pca_data = defaultdict(list)
        sample_list = []
        scatter_id = self.sg_scatter(task_id, origin_id, "", "1", "pop_pca_scatter")
        with open(file_path, 'r') as r:
            data = r.readlines()
            for line in data:
                tmp = line.strip().split(' ')
                sample_list.append(tmp[0])
                for i in range(2, len(tmp)):
                    pca_data[i - 1].append(round(float(tmp[i]), 3))
        for key in pca_data.keys():
            self.sg_scatter_detail(scatter_id, "PCA{}".format(key), pca_data[key])
        self.update_db_record("sg_scatter", {"_id": scatter_id}, {"sample_list": sample_list})

    def add_pca_curve(self, task_id, origin_id, file_path):
        """
        PCA累积变量解释度图-- pop.pca.eigenval
        :param task_id:
        :param origin_id:
        :param file_path:
        :return:
        """
        values = []
        categories = []
        origin_id = self.check_objectid(origin_id)
        i = 1
        with open(file_path, 'r') as r:
            data = r.readlines()
            for n in data:
                temp = n.strip().split('\t')
                categories.append(i)
                values.append(round(float(temp[0]), 3))
                i += 1
        curve_id = self.sg_curve(task_id, origin_id, '', categories, '1', "pop_pca_curve")
        self.sg_curve_detail(curve_id, '', values)

    def add_pca_bar(self, task_id, origin_id, file_path):
        """
        PCA累积变量解释度图-- pop.pca.eigenval
        :param task_id:
        :param origin_id:
        :param file_path:
        :return:
        """
        values = []
        categories = []
        origin_id = self.check_objectid(origin_id)
        i = 1
        with open(file_path, 'r') as r:
            data = r.readlines()
            for n in data:
                temp = n.strip().split('\t')
                categories.append(i)
                values.append(round(float(temp[0]), 3))
                i += 1
                if i == 21:
                    break
        bar_id = self.sg_bar(task_id, origin_id, '', categories, '4', "pop_pca_bar")
        sum_ = sum(values)
        new_value = []
        for m in values:
            new_value.append(m / sum_ * 100)
        self.sg_bar_detail(bar_id, 'PCA', new_value)

    def add_sg_pop_structure(self, task_id, project_sn, params):
        """
        添加structure的分析主表
        :param task_id:
        :param project_sn:
        :param params:
        :return:
        """
        if params:
            params_ = params
        else:
            params_ = {}
        main_id = self.add_main_table("sg_pop_structure", task_id, project_sn, params_, "origin_structure_analysis",
                                      "开始进行structure分析")
        return main_id

    def add_sg_kvalue_detail(self, structure_id, file_dir, task_id):
        """
        添加群体结构数据表--pop.14.xls, 以及群体结构直方图
        :param structure_id:
        :param file_dir:
        :param task_id:
        :return:
        """
        structure_id = self.check_objectid(structure_id)
        for m in os.listdir(file_dir):
            n = re.match("pop\.(.*)\.xls$", m)
            if n:
                self.make_kvalue(os.path.join(file_dir, m), structure_id, n.group(1), task_id)

    def make_kvalue(self, file_path, structure_id, k, task_id):
        data_list = []
        k_data = defaultdict(list)
        with open(file_path, "r") as r:
            for line in r:
                tmp = line.strip().split('\t')
                insert_data = {
                    "structure_id": structure_id,
                    "sample_id": tmp[0],
                    "name": "K{}".format(k)
                }
                tmp1 = tmp[1].strip().split(' ')
                for i in range(0, len(tmp1)):
                    k_data[tmp[0]].append(round(float(tmp1[i]), 3))
                    insert_data.update({"Q{}".format(i + 1): round(float(tmp1[i]), 3)})
                data_list.append(insert_data)
        if data_list:
            self.col_insert_data("sg_kvalue_detail", data_list)
        else:
            self.bind_object.logger.info("文件{}为空不进行导表".format(file_path))
        bar_id = self.sg_bar(task_id, structure_id, "K{}".format(k), [], '5', "pop_structure_bar")  # 导入群体结构直方图
        for key in k_data.keys():
            self.sg_bar_detail(bar_id, key, k_data[key], "false")

    def add_structure_sg_curve(self, task_id, origin_id, file_path):
        """
        K 值折线图--cv.error
        :param task_id:
        :param origin_id:
        :param file_path:
        :return:
        """
        categories = []
        values = []
        origin_id = self.check_objectid(origin_id)
        with open(file_path, "r") as r:
            for line in r:
                tmp = line.strip().split('\t')
                categories.append(tmp[0])
                values.append(round(float(tmp[1]), 3))
        curve_id = self.sg_curve(task_id, origin_id, "cverror", categories, "1", "kvalue_curve")
        self.sg_curve_detail(curve_id, "cverror", values, "false")

    def add_pop_sg_tree(self, task_id, origin_id, file_path, tree_type='nj'):
        """
        添加进化树--pop.nj.tree,
        :param task_id:
        :param tree_type: nj ml bayes
        :param file_path:
        :param file_path:
        :param origin_id:
        :return:
        """
        if tree_type not in ['nj', 'ml', 'bayes']:
            self.set_error("tree_type 必须是nj, ml, bayes中的一种")
        tree_str = ''
        origin_id = self.check_objectid(origin_id)
        with open(file_path, "r") as r:
            data = r.readlines()
            if tree_type in ['nj', 'ml']:
                tree_str = data[0].strip()
            else:
                for line in data:
                    if re.search('.*con_50_majrule.*', line):
                        tmp = line.strip().split(' ')
                        tree_str = tmp[3]
                        break
        tree_id = self.add_sg_tree(task_id, origin_id, '{}_tree'.format(tree_type), 2, 'pop_tree')
        self.add_sg_tree_detail(tree_id, '{}_tree'.format(tree_type), tree_str)

    def add_sg_pop_tree(self, task_id, project_sn, params):
        """
        导入sg_pop_tree主表
        :param task_id:
        :param project_sn:
        :param params:
        :return:
        """
        if params:
            params_ = params
        else:
            params_ = {}
        main_id = self.add_main_table("sg_pop_tree", task_id, project_sn, params_, "origin_tree_analysis",
                                      "开始进行tree分析")
        return main_id


if __name__ == "__main__":
    a = PopAnalysis(None)
    task = "tsg_32990"
    proj = "evolution_test"
    # main_id_ = a.add_sg_pop_pca("evolution_test", "evolution_test", {})
    # a.add_sg_pop_pca_detail(main_id_,
    #                         "/mnt/ilustre/users/sanger-dev/workspace/20180911/PopAnalysis_workflow_pop_0911_pca"
    #                         "/output/pca/pop.pca.eigenvec")
    # a.add_sg_scatter("evolution_test", "5badce48a4e1af7ee96c2ce7",
    #                  "/mnt/ilustre/users/sanger-dev/workspace/20180911/PopAnalysis_workflow_pop_0911_pca/"
    #                  "output/pca/pop.pca.eigenvec")
    # str_id = a.add_sg_pop_structure(task, proj, {})
    # a.add_sg_kvalue_detail(str_id, "/mnt/ilustre/users/sanger-dev/workspace/20180911/"
    #                                "PopAnalysis_workflow_pop_0911_structure/output/structure/structure", task)
    # test_id = a.add_sg_pop_tree(task, proj, {})
    # a.add_pop_sg_tree(task, test_id, "/mnt/ilustre/users/sanger-dev/workspace/20180911/"
    #                                  "PopAnalysis_workflow_pop_0911_tree/output/tree/tree/pop.nj.tree", "nj")
    # a.add_pop_sg_tree(task, test_id, "/mnt/ilustre/users/sanger-dev/workspace/20180911/"
    #                                  "PopAnalysis_workflow_pop_0911_mltree/output/tree/tree/"
    #                                  "pop.phylip.raxml.bestTree", "ml")
    # a.add_pop_sg_tree(task, test_id, "/mnt/ilustre/users/sanger-dev/workspace/20180911/"
    #                                  "PopAnalysis_workflow_pop_0911_bayes/output/tree/tree/bayes.nex.con", "bayes")
    # a.add_pca_curve(task, "5c072fb3a4e1af1f7ca4b427",
    #                 "/mnt/ilustre/users/sanger-dev/workspace/20181204/Evolution_tsg_32990/output/06.pop/pop.pca.eigenval")
    # a.add_pca_bar("tsg_32120", "5bc459a6a4e1af4405399d67",
    #               "/mnt/ilustre/users/sanger-dev/workspace/20181212/Evolution_tsg_33016/GctaPca/output/pop.pca.eigenval")
    # a.make_kvalue("/mnt/ilustre/users/sanger-dev/workspace/20181210/Evolution_tsg_33005/PopStructure/output/"
    #               "structure/pop.2.xls", "5bc459a6a4e1af4405399d67", '12', "test")
    a.add_sg_kvalue_detail("5bc459a6a4e1af4405399d67", "/mnt/ilustre/users/sanger-dev/workspace/20181210/"
                                                       "Evolution_tsg_33005/PopStructure/output/structure/", "ttt")
    print "ok"




