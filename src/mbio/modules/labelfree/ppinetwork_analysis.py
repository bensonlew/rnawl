# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify:20170315

from biocluster.module import Module
import os
import types
from biocluster.core.exceptions import OptionError
import unittest
import time

class PpinetworkAnalysisModule(Module):
    def __init__(self, work_id):
        super(PpinetworkAnalysisModule, self).__init__(work_id)
        self.step.add_steps('map', 'map_blast', 'ppinetwork_predict', 'ppinetwork_topology')
        options = [
            {"name": "diff_exp_gene", "type": "infile", "format": "labelfree.common"},
            {"name": "species", "type": "int", "default": 9606},
            {"name": "combine_score", "type": "int", "default": 300},
            {"name": "seq", "type": "string", "default": ""}
        ]
        self.add_option(options)
        self.map = self.add_tool("labelfree.proteinset.ppinetwork_map")
        self.map_blast = self.add_tool("labelfree.proteinset.ppinetwork_mapblast")
        self.ppinetwork_predict = self.add_tool("labelfree.proteinset.ppinetwork_predict")
        self.ppinetwork_topology = self.add_tool("labelfree.proteinset.ppinetwork_topology")
        self._end_info = 0

    def check_options(self):
        species_list = [30611, 9598, 61853, 9593, 9606, 9544, 9483, 30608, 9601, 9478, 10141, 10020, 10090, 9986, 10116,
                        43179, 37347, 9685, 9913, 9739, 9669, 9796, 132908, 59463, 9646, 9823, 9785, 9813, 9371, 9361,
                        28377, 9031, 13735, 9103, 59729, 8049, 31033, 8090, 8083, 69293, 99883, 8128, 7955, 13616, 9258,
                        9305, 9315, 7897, 7757, 7719, 51511, 6239, 7227, 4932, 15368, 4513, 4641, 4533, 4538, 4555,
                        4558, 4577, 59689, 3702, 3711, 3847, 3694, 4081, 4113, 29760, 88036, 3218, 3055, 45157]
        if not self.option("diff_exp_gene").is_set:
            raise OptionError("必须输入含有gene_id的差异基因列表！")
        if not isinstance(self.option('combine_score'), int) or self.option('combine_score') < 0:
            raise OptionError("combined_score值必须是大于0的整数！")
        '''
        if int(self.option('species')) not in species_list:
            raise OptionError("不能进行蛋白质互作分析，因为string数据库中不存在该物种的蛋白质互作组数据！")
        '''
        if self.option('seq') == "":
            raise OptionError("seq参数不可以为空")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def map_run(self):
        self.map.set_options({
            "diff_exp_gene": self.option("diff_exp_gene"),
            "species": self.option("species")
        })
        self.map.on('end', self.set_output, 'map')
        self.map.on('start', self.set_step, {'start': self.step.map})
        self.map.on('end', self.set_step, {'end': self.step.map})
        self.map.on('end', self.map_blast_run)
        self.map.run()

    def map_blast_run(self):
        # mapped_table = os.path.join(self.work_dir, "PpinetworkMap/output/seq_mapped.txt")
        mapped_table = self.map.option('seq_mapped').prop['path']
        # unmapped_seq = os.path.join(self.work_dir, "PpinetworkMap/output/unmapped_seq.txt")
        unmapped_seq = self.map.option('unmapped_seq').prop['path']
        # unmapped_db = os.path.join(self.work_dir, "PpinetworkMap/output/unmapped_db.txt")
        unmapped_db = self.map.option('unmapped_db').prop['path']
        self.map_blast.set_options({
            "mapped_table": mapped_table,
            "species": self.option("species"),
            "unmapped_seq": unmapped_seq,
            "unmapped_db": unmapped_db,
            "fa": self.option("seq"),
            "blast": "blastp",
            "database": "customer_mode",
            "query_type": "prot",
        })
        self.map_blast.on('end', self.set_output, 'map_blast')
        self.map_blast.on('start', self.set_step, {'start': self.step.map_blast})
        self.map_blast.on('end', self.set_step, {'end': self.step.map_blast})
        self.map_blast.on('end', self.ppinetwork_predict_run)
        self.map_blast.run()

    def ppinetwork_predict_run(self):
        diff_exp_mapped_table = os.path.join(self.map_blast.output_dir, "mapped_all.xls")
        line = open(diff_exp_mapped_table, "r").readlines()[1:]
        if not line:
            raise Exception("基因集中的基因不能匹配到string数据库中id，请您确认基因id为Ensemble或者Entrez GeneID！")
        self.ppinetwork_predict.set_options({
            "diff_exp_mapped": diff_exp_mapped_table,
            "species": self.option("species"),
            "combine_score": self.option("combine_score")
        })
        self.ppinetwork_predict.on('end', self.set_output, 'ppinetwork_predict')
        self.ppinetwork_predict.on('start', self.set_step, {'start': self.step.ppinetwork_predict})
        self.ppinetwork_predict.on('end', self.set_step, {'end': self.step.ppinetwork_predict})
        self.ppinetwork_predict.on('end', self.ppinetwork_topology_run)
        self.ppinetwork_predict.run()

    def ppinetwork_topology_run(self):
        # ppi_table = os.path.join(self.work_dir, "PpinetworkPredict/output/interaction_detail.txt")
        ppi_table = os.path.join(self.ppinetwork_predict.work_dir, "output/interaction_detail.txt")
        # node_table = os.path.join(self.work_dir, "PpinetworkPredict/output/all_nodes.txt")
        node_table = os.path.join(self.ppinetwork_predict.work_dir, "output/all_nodes.txt")
        if not os.path.exists(ppi_table) or not os.path.exists(node_table):
            super(PpinetworkAnalysisModule, self).end()
            return
        p = open(ppi_table)
        n = open(node_table)
        if len(p.readlines()) == 1 or len(n.readlines()) == 1:
            super(PpinetworkAnalysisModule, self).end()
            return
        time.sleep(2) #好像有延迟，有时无法正常结束
        self.ppinetwork_topology.set_options({
            "ppitable": ppi_table,
            "nodetable": node_table,
            "combine_score": 0
        })
        self.ppinetwork_topology.on('end', self.set_output, 'ppinetwork_topology')
        self.ppinetwork_topology.on('end', self.end)
        self.ppinetwork_topology.run()

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
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'map':
            self.linkdir(obj.output_dir, 'map')
        elif event['data'] == 'ppinetwork_predict':
            self.linkdir(obj.output_dir, 'ppinetwork_predict')
        elif event['data'] == 'ppinetwork_topology':
            self.linkdir(obj.output_dir, 'ppinetwork_topology')
        else:
            pass

    def run(self):
        super(PpinetworkAnalysisModule, self).run()
        mapped_table = os.path.join(self.work_dir, "PpinetworkMap/output/diff_exp_mapped.txt")
        if os.path.isfile(mapped_table):
            self.ppinetwork_predict_run()
        else:
            self.map_run()


    def end(self):
        repaths = [
            [".", "", "蛋白质互作网络结果输出目录"],
            ["interaction.txt", "txt", "edges结果文件信息"],
            ["all_nodes.txt", "txt", "node结果信息"],
            ["network_stats.txt", "txt", "网络统计结果信息"],
            ["diff_exp_mapped.txt", "txt", "含有STRINGid的差异基因文件"],
            ["protein_interaction_network_centrality.txt", "txt", "PPI网络中心系数表"],
            ["protein_interaction_network_clustering.txt", "txt", "PPI网络节点聚类系数表"],
            ["protein_interaction_network_transitivity.txt", "txt", "PPI网络传递性"],
            ["protein_interaction_network_degree_distribution.txt", "txt", "PPI网络度分布表"],
            ["protein_interaction_network_by_cut.txt", "txt", "combined_score值约束后的PPI网络"],
            ["protein_interaction_network_node_degree.txt", "txt", "PPI网络节点度属性表"]
        ]

        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        # sdir.add_regexp_rules(regexps)
        super(PpinetworkAnalysisModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_itraq_and_tmt/data4'
        data = {
            "id": "ppinetwork" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "itraq_and_tmt.ppinetwork_analysis",
            "instant": False,
            "options": dict(
                diff_exp_gene=test_dir + "/" + "test_set1.list",
                species=4081,
                combine_score=400
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
