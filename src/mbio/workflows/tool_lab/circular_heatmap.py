# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
import time
import unittest
from biocluster.file import getsize, exists
from biocluster.file import download


class CircularHeatmapWorkflow(Workflow):
    """
    代谢样本相关性分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CircularHeatmapWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "group_detail", "type": "string"},
            {"name": "metab_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "sam_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "group_method", "type": "string", "default": "no"},
            {"name": "metab_dist", "type": "string", "default": "euclidean"},
            {"name": "metab_cluster", "type": "string", "default": "complete"},
            {"name": "n_cluster", "type": "int"},
            {"name": "sam_dist", "type": "string", "default": "euclidean"},
            {"name": "sam_cluster", "type": "string", "default": "complete"},
            {"name": "scale", "type": "string", "default": "scale"},  ## 是否标准化聚类数据，scale, unscale
            {'name': 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'relate_id', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},

        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.profile = self.add_tool("tool_lab.circular_heatmap.select_table")
        self.metab_cluster = self.add_tool("tool_lab.circular_heatmap.metab_cluster")

    def run(self):
        self.profile.on('end', self.run_cluster)
        self.metab_cluster.on('end', self.set_db)
        if self.option('source') == 'project':
            self.metab_path, self.group_path = self.check_file_path()
            self.metab_file = self.download_s3_file(self.metab_path, 'metab_table.txt')
            self.group_file = self.download_s3_file(self.group_path, 'group_table.txt')
        if self.option('source') == 'tool_lab':
            self.metab_file = self.option('metab_table').prop['path']
            self.group_file = self.option('group_table').prop['path']
        self.select_profile()
        super(CircularHeatmapWorkflow, self).run()

    def select_profile(self):
        self.logger.info("start profile!")
        # exp_profile = self.option("metab_table").prop["path"]
        if self.option("group_method") == "average":
            group_method = 2
        elif self.option("group_method") == "median":
            group_method = 3
        elif self.option("group_method") == "sum":
            group_method = 1
        elif self.option("group_method") == "no":
            group_method = 0
        else:
            self.set_error("group_method方法不对-%s", variables=(self.option("group_method")), code="14700701")
        options = {
            "origin_table": self.metab_file,
            "group_method": group_method,
            "group": self.group_file,
        }
        if self.option("scale") == "scale":
            options["scale"] = True
        else:
            options["scale"] = False
        self.profile.set_options(options)
        self.profile.run()

    def run_cluster(self):
        self.logger.info("start run metabset_cluster !")
        profile = self.profile.option("select_table")
        # metab_abu = self.option("metab_table").prop["path"]
        # metab_des = metab_abu.replace("metab_abund.txt","metab_desc.txt")
        metab_des = self.profile.option("map_table").prop["path"]
        if not os.path.exists(metab_des):
            self.set_error('metab_trans路径不存在-%s!', variables=(metab_des), code="14700702")
        options = {
            'exp': profile,
            # 'sct': self.option("sam_cluster_method"),
            # 'mct': self.option("metab_cluster_method"),
            'metab_trans': metab_des
        }
        if self.option('sam_cluster_method') == 'no':
            options['sct'] = ''
        else:
            options['sct'] = self.option("sam_cluster_method")
        if self.option('metab_cluster_method') == 'no':
            options['mct'] = ''
        else:
            options['mct'] = self.option("metab_cluster_method")
        if self.option("metab_cluster_method") == "hierarchy":
            options['mcm'] = self.option("metab_cluster")
            options['mcd'] = self.option("metab_dist")
            options['n_cluster'] = self.option("n_cluster")
        if self.option("sam_cluster_method") == "hierarchy":
            options['scm'] = self.option("sam_cluster")
            options['scd'] = self.option("sam_dist")
        if self.option("metab_cluster_method") == "kmeans":
            options['n_cluster'] = self.option("n_cluster")
            options['mcd'] = self.option("metab_dist")
        if self.option("scale") == "scale":
            options["before_scale"] = self.profile.option("select_origin_abu")
        self.logger.info(options)
        self.metab_cluster.set_options(options)
        self.metab_cluster.run()

    def check_file_path(self):
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        conn_upset = db[collection_name]
        status = 'start'
        count_time = 0
        while status == 'start':
            if count_time > 600:
                self.set_error('超过十分钟还没有结果文件生成，请检查是否生成文件时报错')
                break
            time.sleep(10)
            print 'sleep 10s'
            try:
                upset = conn_upset.find_one(
                    {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
                status = upset['status']
            except:
                pass
            count_time += 10
        upset = conn_upset.find_one(
            {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
        metab_path = upset['metab_path']
        group_path = upset['group_path']
        return metab_path, group_path

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path,), code='13700502')
        return to_path

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_circular = self.api.api("tool_lab.circular_heatmap")
        sam_tree = os.path.join(self.metab_cluster.output_dir, "col.cluster_tree.xls")
        metab_tree = os.path.join(self.metab_cluster.output_dir, "row.cluster_tree.xls")
        if os.path.exists(sam_tree):
            sam_tree = sam_tree
        else:
            sam_tree = None
        if os.path.exists(metab_tree):
            metab_tree = metab_tree
        else:
            metab_tree = None
        if self.option("scale") == "scale":
            expression_file = os.path.join(self.metab_cluster.output_dir, "cluster_scale_exp.xls")
        else:
            expression_file = os.path.join(self.metab_cluster.output_dir, "cluster_exp.xls")
        api_circular.add_circular_heatmap(sampletree=sam_tree, featuretree=metab_tree, heatmap=expression_file,
                                          group=self.group_file, main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.metab_cluster.output_dir)
        relpath_rules =[
            [".", "", "圆圈热图结果目录", 0, "150027"],
            ["col.cluster_tree.xls", "xls", "column聚类树文件", 0, "150029"],
            ["row.cluster_tree.xls", "xls", "row聚类树文件", 0, "150028"],
            ["cluster_exp.xls", "xls", "表达量文件", 0, "150030"],
            ["row.kmeans_cluster.xls", "xls", "row kmeans分类文件", 0, "150065"],
            ["cluster_scale_exp.xls", "xls", "标准化后表达量文件", 0, "150077"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        # regexps = [
        #     [r"row.subcluster_.*\.xls", "xls", "row kmeans各子类结果表", 0, "150066"],
        # ]
        # result_dir.add_regexp_rules(regexps)
        super(CircularHeatmapWorkflow, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the workflow. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitoollabtest.py '
        cmd += 'post toollabpipeline '
        cmd += '-c {} '.format("client03")
        cmd += "-b http://bcl.tsg.com "
        cmd += "-n \"params;basis\" -d \"{"
        args = dict(
            metab_table='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_062021/dataframe.txt',
            group_table='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_062021/group.txt',
            # metab_table='/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210608/MetabsetCluster_65ch_oqs584svkuvq74854gudg6_17727_483353/remote_input/metab_table/metab_abund.txt',
            # group_table='/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210608/MetabsetCluster_65ch_oqs584svkuvq74854gudg6_17727_483353/group_table_input.group.xls',
            # group_method='no',
            scale='unscale'
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.circular_heatmap",
            main_table_name="sg_circular_heatmap",
            task_id="circular_heatmap",
            project_sn="circular_heatmap",
            submit_location="circular_heatmap"
        )
        for arg in args:
            cmd += "\\\""
            cmd += arg
            cmd += "\\\":\\\""
            cmd += args[arg]
            cmd += "\\\","
        cmd = cmd.rstrip(",")
        cmd += "};{"
        for arg in config:
            cmd += "\\\""
            cmd += arg
            cmd += "\\\":\\\""
            cmd += config[arg]
            cmd += "\\\","
        cmd = cmd.rstrip(",")
        cmd += "}\""

        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()