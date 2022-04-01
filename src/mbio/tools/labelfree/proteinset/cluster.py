# -*- coding: utf-8 -*-
# __author__ = "qiuping"
# last_modify:2016.10.09

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna_v2.cluster import clust
import os
# 当前未使用该脚本--gdq

class ClusterAgent(Agent):
    """
    version v1.0
    author: qiuping
    last_modify: 2016.06.20
    """
    def __init__(self, parent):
        super(ClusterAgent, self).__init__(parent)
        options = [
            {"name": "diff_fpkm", "type": "infile", "format": "rna.express_matrix"},  # 输入文件，差异基因表达量矩阵
            {"name": "samples_distance_method", "type": "string","default":"average"},  # 计算距离的算法hclust对应的是complete single average kmeans对应的是euclidean
            {"name": "genes_distance_method", "type": "string","default":"complete"},  # 计算距离的算法hclust对应的是complete single average kmeans对应的是euclidean
            {"name": "samples_distance_algorithm","type":"string","default":"pearson"}, #基因距离算法，只对hclust，默认是euclidean
            {"name": "genes_distance_algorithm", "type":"string","default":"euclidean"},  #样本聚类算法，只对hclust, 默认是pearson
            {"name": "log", "type": "int", "default": 10},  # 画热图时对原始表进行取对数处理，底数为10或2
            {"name": "method", "type": "string", "default": "hclust"},  # 聚类方法选择
            {"name": "sub_num", "type": "int", "default": 10},  # 子聚类的数目
            {"name": "is_genelist","type":"bool","default":False}, #是否设置gene_list参数
            {"name": "diff_list_dir","type":"infile", 'format':"rna.gene_list_dir"},
            {"name": "group_method", "type":"string", "default":"mean"},
            {"name": "diff_list", "type": "outfile", "format": "rna.gene_list"}

        ]
        self.add_option(options)
        self.step.add_steps("cluster")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.cluster.start()
        self.step.update()

    def stepfinish(self):
        self.step.cluster.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("diff_fpkm").is_set:
            raise OptionError("必须设置输入文件:差异基因fpkm表")
        #if self.option("distance_method") not in ("manhattan", "euclidean"):
        #    raise OptionError("所选距离算法不在提供的范围内")
        if self.option('log') not in (10, 2):
            raise OptionError("所选log底数不在提供的范围内")
        if self.option("method") not in ("hclust", "kmeans", "both"):
            raise OptionError("所选方法不在范围内")
        if not isinstance(self.option("sub_num"), int):
            raise OptionError("子聚类数目必须为整数")
        if not (self.option("sub_num") >= 1 and self.option("sub_num") <= 35):
            raise OptionError("子聚类数目范围必须在3-35之间！")

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '50G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option('method') in ('both', 'hclust'):
            result_dir.add_relpath_rules([
                [".", "", "结果输出目录"],
                ["./hclust/", "", "hclust分析结果目录"],
                ["hclust/hc_gene_order", "txt", "按基因聚类的基因排序列表"],
                ["hclust/hc_sample_order", "txt", "按样本聚类的样本排序列表"],
                ["hclust/hclust_heatmap.xls", "xls", "层级聚类热图数据"]
            ])
            result_dir.add_regexp_rules([
                [r"hclust/subcluster_", "xls", "子聚类热图数据"]
            ])
        if self.option('method') in ('both', 'kmeans'):
            result_dir.add_relpath_rules([
                [".", "", "结果输出目录"],
                ["./kmeans/", "", "kmeans分析结果目录"],
                ["kmeans/kmeans_heatmap.xls", "xls", "层级聚类热图数据"]
            ])
            result_dir.add_regexp_rules([
                [r"kmeans/subcluster_", "xls", "子聚类热图数据"]
            ])
        super(ClusterAgent, self).end()


class ClusterTool(Tool):
    """
    表达量差异检测tool
    """
    def __init__(self, config):
        super(ClusterTool, self).__init__(config)
        self._version = '1.0.1'
        self.r_path = '/program/R-3.3.1/bin/Rscript'

    def filter_gene(self, diff_list_file):
        if self.option("is_genelist"):
            gene_list = []
            with open(diff_list_file,'r+') as m1:
                for lines in m1:
                    gene_id = lines.strip()
                    if gene_id not in gene_list:
                        gene_list.append(gene_id)

            tmp = self.work_dir + "/tmp"
            sequence = []
            with open(self.option("diff_fpkm").prop['path'],'r+') as f1,open(tmp,'w+') as f2:
                ll = f1.readline()
                samples=ll.strip().split("\t")
                f2.write(ll)
                for lines in f1:
                    line = lines.strip().split("\t")
                    if line[0] in gene_list:
                        f2.write(lines)
                        sequence.append(lines)
                    else:
                        continue
            if not sequence:
                raise Exception("{}对应的genelist不在{}中出现！".format(diff_list_file,self.option('diff_fpkm').prop['path']))
            else:
                return tmp
        else:
            raise Exception("过滤基因或转录本id时请设置is_genelist为true！")

    def run_cluster(self,fpkm_path):
        #genes_distance_method=self.option("genes_distance_method")
        #samples_distance_method=self.option("samples_distance_method")
        clust(input_matrix=fpkm_path, sub_num=self.option('sub_num'), method=self.option('method'), lognorm=self.option('log'), samples_distance_method=self.option('samples_distance_method'), genes_distance_method = self.option("genes_distance_method"), cltype="both",
              samples_distance_algorithm=self.option("samples_distance_algorithm"),genes_distance_algorithm = self.option("genes_distance_algorithm"))
        clust_cmd = self.r_path + " clust.r"
        self.logger.info("开始运行clust_cmd")
        cmd = self.add_command("clust_cmd", clust_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行clust_cmd成功")
        else:
            self.set_error("运行clust_cmd出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        try:
            if self.option('method') in ('both', 'hclust'):
                os.system('cp -r %s/hclust/ %s/' % (self.work_dir, self.output_dir))
                self.logger.info("设置hclust结果目录成功")
            if self.option('method') in ('both', 'kmeans'):
                os.system('cp -r %s/kmeans/ %s/' % (self.work_dir, self.output_dir))
                self.logger.info("设置kmeans结果目录成功")
            self.logger.info("设置聚类分析结果目录成功")
        except Exception as e:
            self.logger.info("设置聚类分析结果目录失败{}".format(e))

    def run(self):
        super(ClusterTool, self).run()
        if not self.option("sub_num"):
            if len(self.option("diff_fpkm").prop['gene']) > 200:
                self.option("sub_num", 10)
            else:
                self.option("sub_num", 5)
        if self.option("is_genelist"):
            if self.option("diff_list_dir").is_set:
                self.file_path = []
                for files in os.listdir(self.option("diff_list_dir").prop['path']):
                    diff_path = os.path.join(self.option("diff_list_dir").prop['path'], files)
                    if os.path.getsize(diff_path) != 0:
                        if diff_path not in self.file_path:
                            self.file_path.append(diff_path)
                if self.file_path:
                    fpkm_path = self.filter_gene(self.file_path[0])
                    self.option("diff_list").set_path(self.file_path[0])
                    self.run_cluster(fpkm_path)
                else:
                    self.logger.info("{}没有生成差异基因集，无法根据差异基因进行聚类分析！".format(self.option("diff_list_dir").prop['path']))
        else:
            self.run_cluster(self.option("diff_fpkm").prop['path'])
        self.set_output()
        self.end()
