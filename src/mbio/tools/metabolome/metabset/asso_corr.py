# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.scripts.corr_cluster import CorrClust
import subprocess
import os, glob
import re
import pandas as pd
from mbio.packages.metabolome.common import Relation


class AssoCorrAgent(Agent):
    """
    pearsonsCorrelation:用于生成环境因子和otu/taxon之间的correlation
    version: 0.1
    author: wangbixuan
    last_modified: 20160930 by qindanhua
    """
    def __init__(self, parent):
        super(AssoCorrAgent, self).__init__(parent)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.express,sequence.profile_table",
             "required": True},
            {"name": "asso_table", "type": "infile", "format": "sequence.profile_table", "required": True},
            {"name": "coefficient", "type": "string", "default": "pearson"},
            {'name': 'sct', 'type': 'string', 'default': ''},  # 关联数据聚类算法，hierarchy、无
            {'name': 'scd', 'type': 'string', 'default': ''},  # 关联数据距离计算方式
            {'name': 'scm', 'type': 'string', 'default': ''},  # 关联数据聚类方式,"complete","average","single"
            {'name': 'mct', 'type': 'string', 'default': ''},  # 代谢物聚类算法，hierarchy、kmeans、无
            {'name': 'mcd', 'type': 'string', 'default': ''},  # 代谢物距离计算方式
            {'name': 'metab_n_cluster', 'type': 'int', 'default': 0},  # 代谢物聚类数目，kmeans时使用
            {'name': 'mcm', 'type': 'string', 'default': ''},  # 代谢物聚类方式, hierarchy时使用，"complete","average","single"
            {'name': 'asso_n_cluster', 'type': 'int', 'default': 0},  # 代谢物聚类数目，kmeans时使用
            {'name': 'metab_trans', 'type': 'infile', 'format': "metabolome.express,sequence.profile_table"},  # 转化id使用
        ]
        self.add_option(options)
        self.corrfile = ""
        self.pvaluefile = ""

    def check_options(self):
        if self.option("coefficient") not in ["pearson", "spearman", "kendall"]:
            raise OptionError('不支持该相关系数方法', code="34700801")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(AssoCorrAgent, self).end()


class AssoCorrTool(Tool):
    def __init__(self, config):
        super(AssoCorrTool, self).__init__(config)
        '''
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.perl_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/"
        self.r_path = self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/Rscript'
        self.hcluster_script_path = self.config.SOFTWARE_DIR + "/bioinfo/statistical/scripts/"
        self.Rscript_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin/"
        '''
        self.cmd_path = '{}/miniconda2/bin/python {}/statistical/pearsonsCorrelation.py'\
            .format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)
        self.python_path = 'miniconda2/bin/python'
        self.cluster_script = self.config.PACKAGE_DIR + "/metabolome/scripts/corr_cluster.py"

    def run(self):
        """
        运行
        """
        super(AssoCorrTool, self).run()
        self.run_pearsonsCorrelation()
        self.logger.info(self.option("sct"))
        self.logger.info(self.option("mct"))
        if self.option("sct") or self.option("mct"):
            self.run_heatmap()
        self.set_output()
        self.end()

    def run_pearsonsCorrelation(self):
        """
        run pearsonsCorrelation.py
        """
        cmd = self.cmd_path
        if self.option("coefficient") == "kendall":
            method = "kendalltau"
        elif self.option("coefficient") == "pearson":
            method = "pearsonr"
        elif self.option("coefficient") == "spearman":
            method = "spearmanr"
        metab_table = self.option("metab_table").prop["path"]
        asso_table = self.option("asso_table").prop["path"]
        t_asso_file = os.path.join(self.work_dir, "asso_table.xls")
        t_asso_table = pd.read_table(asso_table, sep="\t", header=0, index_col=0)
        tmp_metab_table = pd.read_table(metab_table, sep="\t", header=0, index_col=0)
        print len(t_asso_table)
        print len(tmp_metab_table)
        if len(t_asso_table) < 2 and self.option("sct")!= "no":
            self.set_error('只有一个关联物数据，无法进行聚类', code="34700801")
        if len(tmp_metab_table) < 2 and self.option("mct")!= "no":
            self.set_error('只有一个代谢物数据，无法进行聚类', code="34700802")
        t_asso_table = t_asso_table.T
        t_asso_table = t_asso_table.ix[:,(t_asso_table !=0 ).any()]  # 去掉全为0的列 zouguanqing
        t_asso_table.to_csv(t_asso_file, sep="\t", index=True, index_label="asso")
        self.corrfile = self.work_dir + "/{}_correlation.xls".format(self.option("coefficient"))
        self.pvaluefile = self.work_dir + "/{}_pvalue.xls".format(self.option("coefficient"))

        cmd += " %s %s %s %s %s" % (metab_table, t_asso_file, self.corrfile, self.pvaluefile, method)
        self.logger.info('运行pearsonsCorrelation.py计算correlation')
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('Pearsons Correlation 计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('Pearsons Correlation 计算失败')
            self.set_error('pearsonsCorrelation.py 计算失败', code="34700803")
            raise Exception('pearsonsCorrelation.py 计算失败')

        self.logger.info('运行pearsonsCorrelation.py计算correlation完成')

    def run_heatmap(self):
        corrtable = pd.read_table(self.corrfile, sep="\t", header=0)
        print len(corrtable)
        if len(corrtable) < 1:
            self.set_error('选择代谢集为空，请尝试切换参数重新运行', code="34700805")
            raise Exception('选择代谢集为空，请尝试切换参数重新运行')
        cmd = '{} {} -exp {} -out {} '.format(self.python_path, self.cluster_script, self.corrfile, self.work_dir)
        if self.option("sct"):
            cmd += " -sct " + self.option("sct")
        else:
            cmd += " --nsc"
        if self.option("mct"):
            cmd += " -gct " + self.option("mct")
        else:
            cmd += " --ngc"
        if self.option("scd"):
            cmd += " -scd " + self.option("scd")
        if self.option("mcd"):
            cmd += " -gcd " + self.option("mcd")
        if self.option("mct") == "kmeans":
            cmd += " -n_clusters " + str(self.option("metab_n_cluster"))
        if self.option("sct") == "hierarchy":
            cmd += " -scm " + self.option("scm")
        if self.option("sct") == "kmeans":
            cmd += " -sn_clusters " + str(self.option("asso_n_cluster"))
        if self.option("mct") == "hierarchy":
            cmd += " -gcm " + self.option("mcm")
        command = self.add_command("cluster", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("run_heatmap succeed")
        else:
            self.set_error("run_heatmap failed!", code="34700807")
            raise Exception("run_heatmap failed!")

    def set_output(self):
        self.logger.info("set output")
        self.logger.info("进行id转化...")
        if self.option("sct") == "hierarchy":
            self.link_file("col.cluster_tree.txt", "asso.cluster_tree.xls")
        elif self.option("sct") == "kmeans":
            subclusters = glob.glob(self.work_dir + '/col.subcluster*')
            for each in subclusters:
                if os.path.exists(each):
                    each = each.split("/")[-1]
                    newname = each.replace("col", "asso")
                    self.link_file(each, newname)
            self.link_file("col.kmeans_cluster.txt", "asso.kmeans_cluster.xls")
        table_path = self.option("metab_trans").prop["path"]
        self.metab_trans = Relation()
        map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path)
        if self.option("mct") == "hierarchy":
            self.trans_tree_name(id_name_dict, "row.cluster_tree.txt", "metab.cluster_tree.xls")
            self.link_file("row.cluster_tree.txt", "metab_id.cluster_tree.xls")
        elif self.option("mct") == "kmeans":
            subclusters = glob.glob(self.work_dir + '/row.subcluster*')
            for each in subclusters:
                each = each.split("/")[-1]
                newname = each.replace("row", "metab")
                self.trans_data(each, newname, map_table, id_name_dict)
            #self.link_file("row.kmeans_cluster.txt","metab.kmeans_cluster.xls")
            self.metab_trans.trans_col_data(id_name_dict, self.work_dir + "/row.kmeans_cluster.txt",
                                            self.output_dir + "/metab.kmeans_cluster.xls")
        final_corrfile, final_pvaluefile = self.reorder(self.corrfile, self.pvaluefile)
        self.trans_data(final_corrfile, "corr.xls", map_table, id_name_dict)
        self.trans_data(final_pvaluefile, "pvalue.xls", map_table, id_name_dict)
        self.end()

    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.work_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        if os.path.exists(oldfile):
            os.link(oldfile, newfile)

    def trans_data(self, oldfile, newfile, map_table, id_name_dict):
        #oldfile = os.path.join(self.work_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        self.metab_trans.get_trans_file(map_table, id_name_dict, oldfile, newfile)

    def trans_tree_name(self, namedict, treefile, newfile):
        """
        用正则的方式替换复原树文件中的名称
        """
        treefile = os.path.join(self.work_dir, treefile)
        newfile = os.path.join(self.output_dir, newfile)
        if not isinstance(namedict, dict):
            self.set_error('复原树的枝名称需要旧名称和当前名称的字典', code="34700809")
        if os.path.exists(treefile):
            try:
                new_names = []
                with open(treefile, 'rb') as f, open(newfile, 'wb') as w:
                    tree = f.readline().strip()
                    for item in namedict.iteritems():
                        tree = re.sub(item[0] + ':', item[1] + ':', tree)
                    w.write(tree + "\n")
                    name = f.readline().strip().split(";")
                    for each in name:
                        new_name = namedict[each]
                        new_names.append(new_name)
                    w.write(";".join(new_names) + "\n")
            except IOError, e:
                self.set_error('聚类树文件无法找到或者无法打开：%s', variables=(e), code="34700810")

    def reorder(self, corr_file, p_file):
        col_tree_file = None
        row_tree_file = None
        outDir = self.work_dir
        Corr = CorrClust()
        exp_table = pd.read_table(corr_file, sep='\t', header=0, index_col=0)
        p_exp_table = pd.read_table(p_file, sep='\t', header=0, index_col=0)
        final_corr = outDir + "/corr.xls"
        final_p = outDir + "/pvalue.xls"
        corr_table = exp_table
        p_table = p_exp_table
        if self.option("sct") == "hierarchy" and len(exp_table.columns) > 1:
            col_tree_file = outDir + "/col.cluster_tree.txt"
        if self.option("mct") == "hierarchy" and len(exp_table) > 1:
            row_tree_file = outDir + "/row.cluster_tree.txt"
        if self.option("sct") == "hierarchy" or self.option("mct") == "hierarchy":
            corr_table = Corr.order_sample(corr_table, coltree=col_tree_file, rowtree=row_tree_file)
            p_table = Corr.order_sample(p_table, coltree=col_tree_file, rowtree=row_tree_file)
        if self.option("mct") == "kmeans":
            corr_table = Corr.order_kmeans(corr_table, outDir=outDir, cluster=True)
            p_table = Corr.order_kmeans(p_table, outDir=outDir, cluster=True)
        if self.option("sct") == "kmeans":
            corr_table = Corr.order_kmeans(corr_table, outDir=outDir, cluster_col=True)
            p_table = Corr.order_kmeans(p_table, outDir=outDir, cluster_col=True)
        corr_table.to_csv(final_corr, sep="\t")
        p_table.to_csv(final_p, sep="\t")
        return final_corr, final_p
