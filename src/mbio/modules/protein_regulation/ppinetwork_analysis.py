# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify:20170418

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class PpinetworkAnalysisModule(Module):
    def __init__(self, work_id):
        super(PpinetworkAnalysisModule, self).__init__(work_id)
        self.step.add_steps('ppinetwork_map', 'ppinetwork_predict', 'ppinetwork_topology')
        options = [
            {"name": "diff_exp_gene", "type": "infile", "format": "rna.ppi"},
            {"name": "species", "type": "int", "default": 9606},
            {"name": "combine_score", "type": "int", "default": 300}
        ]
        self.add_option(options)
        self.ppinetwork_map = self.add_tool("protein_regulation.ppinetwork_map")
        self.ppinetwork_predict = self.add_tool("protein_regulation.ppinetwork_predict")
        self.ppinetwork_topology = self.add_tool("protein_regulation.ppinetwork_topology")
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
        if int(self.option('species')) not in species_list:
            raise OptionError("不能进行蛋白质互作分析，因为string数据库中不存在该物种的蛋白质互作组数据！")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def map_run(self):
        self.ppinetwork_map.set_options({
            "diff_exp_gene": self.option("diff_exp_gene"),
            "species": self.option("species")
        })
        self.ppinetwork_map.on('end', self.set_output, 'ppinetwork_map')
        self.ppinetwork_map.on('start', self.set_step, {'start': self.step.ppinetwork_map})
        self.ppinetwork_map.on('end', self.set_step, {'end': self.step.ppinetwork_map})
        self.ppinetwork_map.on('end', self.ppinetwork_predict_run)
        self.ppinetwork_map.run()

    def ppinetwork_predict_run(self):
        diff_exp_mapped_table = os.path.join(self.work_dir, "PpinetworkMap/output/diff_exp_mapped.txt")
        lines = len(open(diff_exp_mapped_table).readlines())
        if lines < 2:
            raise Exception("基因集中的基因不能匹配到string数据库中id，请您确认基因id为Ensemble或者Entrez GeneID！")
            # self.end()
        elif lines <= 15:
            raise Exception("能够比对到string数据库中的基因个数过少，无法进行后面的网络预测及分析，请重新选择基因集！")
        else:
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
        ppi_table = os.path.join(self.work_dir, "PpinetworkPredict/output/interaction.txt")
        self.ppinetwork_topology.set_options({
            "ppitable": ppi_table,
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
        if event['data'] == 'ppinetwork_map':
            self.linkdir(obj.output_dir, 'ppinetwork_map')
        elif event['data'] == 'ppinetwork_predict':
            self.linkdir(obj.output_dir, 'ppinetwork_predict')
        elif event['data'] == 'ppinetwork_topology':
            self.linkdir(obj.output_dir, 'ppinetwork_topology')
        else:
            pass

    def run(self):
        super(PpinetworkAnalysisModule, self).run()
        mapped_table = os.path.join(self.work_dir, "Map/output/diff_exp_mapped.txt")
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
