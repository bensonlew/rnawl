# -*- coding: utf-8 -*-
# __author__ = 'xuxi_20210901'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError
import pandas as pd


class MetabVipWorkflow(Workflow):
    """
    MetabVip base R package MetabVip
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabVipWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "diff_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "group_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "scale", "type": "string", 'default': 'yes'}, # 是否标准化 yes/no
            {"name": "group_method", "type": "string", 'default': "none"}, #分组样本计算 none/median/average/sum
            {"name": "mct", "type": "string", 'default': "none"}, # 代谢物聚类算法，hierarchy/kmeans/none
            {"name": "mcm", "type": "string", 'default': "none"}, # 代谢物层级聚类方式,只当mct为hierarchy时：可选"complete","average","single"
            {"name": "n_cluster", "type": "string", 'default': "5"}, # 当mct为kmeans时，子聚类数 [1-20]
            {"name": "mcd", "type": "string", 'default': "euclidean"}, # 代谢物距离算法 ['euclidean','braycurtis','manhattan']       ,['canberra','chebyshev','cityblock','cosine','dice','euclidean','hamming','jaccard','kulsinski','mahalanobis','matching','minkowski','rogerstanimoto','russellrao','seuclidean','sokalmichener','sokalsneath','sqeuclidean','yule',]
            {"name": "main_id", "type": "string", 'default': ""},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.MetabVip = self.add_tool("tool_lab.metab_vip")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_MetabVip()
        super(MetabVipWorkflow, self).run()

    def check_options(self):
        if not self.option("diff_file").is_set:
            raise OptionError("必须设置输入差异文件")
        if not self.option("group_file").is_set:
            raise OptionError("必须设置输入分组文件")
        diff_file_df = pd.read_table(self.option('diff_file').prop['path'],sep="\t")
        all_columns = diff_file_df.columns.values.tolist()
        if 'Vip_plsda' in all_columns:
            if 'Vip_oplsda' in all_columns:
                raise OptionError("存在两列vip，只能存在一列")
        return True

    def run_MetabVip(self):
        diff_file_df = pd.read_table(self.option('diff_file').prop['path'],sep="\t")
        all_columns = diff_file_df.columns.values.tolist()
        if 'Vip_plsda' in all_columns:
            diff_file_df['Vip_oplsda'] = diff_file_df['Vip_plsda']
        if 'Vip_oplsda' in all_columns:
            diff_file_df['Vip_plsda'] = diff_file_df['Vip_oplsda']
        group_df = pd.read_table(self.option('group_file').prop['path'],sep="\t")
        #
        metab_abund_columns = ["metab_id"]+group_df.iloc[:,0].tolist()
        diff_file_df[metab_abund_columns].to_csv(os.path.join(self.work_dir, "metab_abund.txt"), sep='\t',index=False)
        #
        test_good_columns = ["metab_id","Metabolite","Mode","Formula","m/z","RT (min)","KEGG Compound ID","HMDB_ID","CAS number"]
        test_good_columns_ = [i for i in test_good_columns if i in list(diff_file_df.columns)]
        dff = diff_file_df[test_good_columns_]
        dff[test_good_columns_] = dff[test_good_columns_].astype(object)
        dff.to_csv(os.path.join(self.work_dir, "metab_desc.txt"), sep='\t',index=False)
        #
        if not os.path.exists(os.path.join(self.work_dir, "linshi_dir")):
            os.mkdir(os.path.join(self.work_dir, "linshi_dir"))
        diff_file_name = "_vs_".join(list(set(group_df.iloc[:,1].tolist()))[::-1])+".diff.exp.xls"
        diff_file_columns = ["Metabolite","Vip_plsda","Vip_oplsda","P_value","fdr","FC"]+group_df.iloc[:,0].tolist()
        diff_file_columns2 = ["Vip_plsda","Vip_oplsda","P_value","fdr","FC"]+group_df.iloc[:,0].tolist()
        df = diff_file_df[diff_file_columns]
        df[diff_file_columns2] = df[diff_file_columns2].astype(float)
        df.to_csv(os.path.join(self.work_dir, "linshi_dir",diff_file_name), sep='\t',index=False)
        #
        opts = {
            "diff_dir": os.path.join(self.work_dir, "linshi_dir"),
            "metab_trans":os.path.join(self.work_dir, "metab_desc.txt"),
            "group":self.option('group_file').prop['path'],
            "group_method":self.option('group_method'),
            "mct":self.option('mct'),
            "mcm":self.option('mcm'),
            "n_cluster":int(self.option('n_cluster')),
            "mcd":self.option('mcd'),
            "scale":{"yes":True,"no":False}[self.option('scale')],
            "metab_table":os.path.join(self.work_dir, "metab_abund.txt")
        }
        print "xx opts is:"+str(opts)
        self.MetabVip.set_options(opts)
        self.MetabVip.on('end', self.set_db)
        self.MetabVip.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        # api_name = self.api.api("metabolome.metab_vip")
        api_name = self.api.api('tool_lab.metab_vip')
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")
        vip_dir = self.MetabVip.output_dir
        self.logger.info(vip_dir)
        print "main_id is:"+str(main_id)
        api_name.add_metab_vip(main_id, main_id=main_id)
        #table_type = self.option("table_type")
        api_name.add_metab_vip_detail(main_id, vip_dir, "oplsda", os.path.join(self.work_dir, "metab_desc.txt"),
                                         scale=self.option("scale"),out_dir=vip_dir)
        self.output_dir = self.MetabVip.output_dir
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "VIP分析结果文件夹", 0, "150040"]
        ]
        regexps = [
            [r".*/metab.kmeans_cluster.xls", "", "代谢物kmeans分类文件", 0, "150065"],
            [r".*/metab.cluster_tree.xls", "xls", "代谢物树文件", 0, "150028"],
            [r".*/Vip_exp.xls", "xls", "VIP值表", 0, "150049"],
            [r".*/metab.subcluster_.*\.xls", "xls", "代谢物kmeans各子类结果表", 0, "150066"],
            [r".*/metab_id.cluster_tree.xls", "xls", "代谢物id聚类树文件", 0, "150075"],
            [r".*/Vip_scale_exp.xls", "xls", "代谢物VIP值与标准化后表达量", 0, "150078"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(MetabVipWorkflow, self).end()


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
            vcf_path='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/bsa/pop.final.vcf',
            # wp='',
            # mp='HQS1',
            mb='XS11_1',
            wb='F44_mix',
        )
        config = dict(
            type="workflow",
            task_type="submit",
            name="tool_lab.bsa",
            main_table_name="sg_bsa",
            task_id="bsa",
            project_sn="bsa",
            submit_location="bsa"
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