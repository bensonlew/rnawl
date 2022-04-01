# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify:20170315

from biocluster.module import Module
import os
import types
from biocluster.core.exceptions import OptionError
import unittest
import pandas as pd


class PpinetworkAnalysisModule(Module):
    def __init__(self, work_id):
        super(PpinetworkAnalysisModule, self).__init__(work_id)
        self.step.add_steps('map', 'map_blast', 'ppinetwork_predict', 'ppinetwork_topology')
        options = [
            {"name": "diff_exp_gene", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "species", "type": "int", "default": 9606},
            {"name": "score", "type": "float", "default": 0.4},
            {"name": "combine_score", "type": "int", "default": 300},
        ]
        self.add_option(options)
        self.map = self.add_tool("tool_lab.ppi.ppinetwork_map")

        self.ppinetwork_predict = self.add_tool("tool_lab.ppi.ppinetwork_predict")
        self.ppinetwork_topology = self.add_tool("tool_lab.ppi.ppinetwork_topology")
        self._end_info = 0
        self.map_name_bool = True
        self.map_blast = True

    def check_options(self):
        species_list = [30611, 9598, 61853, 9593, 9606, 9544, 9483, 30608, 9601, 9478, 10141, 10020, 10090, 9986, 10116,
                        43179, 37347, 9685, 9913, 9739, 9669, 9796, 132908, 59463, 9646, 9823, 9785, 9813, 9371, 9361,
                        28377, 9031, 13735, 9103, 59729, 8049, 31033, 8090, 8083, 69293, 99883, 8128, 7955, 13616, 9258,
                        9305, 9315, 7897, 7757, 7719, 51511, 6239, 7227, 4932, 15368, 4513, 4641, 4533, 4538, 4555,
                        4558, 4577, 59689, 3702, 3711, 3847, 3694, 4081, 4113, 29760, 88036, 3218, 3055, 45157]
        if not self.option("diff_exp_gene").is_set:
            raise OptionError("必须输入含有gene_id的差异基因列表！", code = "23700801")
        if not isinstance(self.option('combine_score'), int) or self.option('combine_score') < 0:
            raise OptionError("combined_score值必须是大于0的整数！", code = "23700802")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def get_file(self):
        self.seq = os.path.join(self.work_dir, 'seq.fasta')
        self.gene_id = os.path.join(self.work_dir, 'gene_id.txt')
        self.id_name = dict()
        self.name_id = dict()
        with open(self.seq, 'w') as s:
            geneset = pd.read_table(self.option('diff_exp_gene').path, header=0, index_col=False, sep='\t')
            for i in range(0, len(geneset)):
                geneset.fillna('', inplace=True)
                if geneset.iloc[i]['gene_id'] == None or geneset.iloc[i]['gene_id'] == '':
                    geneset.iloc[i]['gene_id'] = geneset.iloc[i]['gene_name']
                if geneset.iloc[i]['seq'] != None and geneset.iloc[i]['seq'] != '':
                    if geneset.iloc[i]['gene_id']:
                        s.write('>' + geneset.iloc[i]['gene_id'] + '\n' + geneset.iloc[i]['seq'] + '\n')
                    elif geneset.iloc[i]['gene_name']:
                        s.write('>' + geneset.iloc[i]['gene_name'] + '\n' + geneset.iloc[i]['seq'] + '\n')
                if geneset.iloc[i]['gene_name']:
                    self.id_name.update({geneset.iloc[i]['gene_id']: geneset.iloc[i]['gene_name']})
                    self.name_id.update({geneset.iloc[i]['gene_name']: geneset.iloc[i]['gene_id']})
                else:
                    self.id_name.update({geneset.iloc[i]['gene_id']: ''})
            geneset['gene_id'].to_csv(self.gene_id, header=True, index=False, sep='\t')
        self.map_id_run()

    def map_id_run(self):
        self.map.set_options({
            "diff_exp_gene": self.gene_id,
            "species": self.option("species")
        })
        self.map.on('end', self.set_output, 'map')
        self.map.on('start', self.set_step, {'start': self.step.map})
        self.map.on('end', self.set_step, {'end': self.step.map})
        self.map.on('end', self.map_name_run)
        self.map.run()

    def map_name_run(self):
        unmapped_seq = os.path.join(self.map.output_dir, "unmapped_seq.txt")
        gene_name = os.path.join(self.work_dir, 'unmapped_name.txt')
        with open(unmapped_seq, 'r') as u, open(gene_name, 'w') as g:
            g.write('gene_id'+ '\n')
            for line in u.readlines():
                unmapped_id = line.strip().split('\t')[0]
                unmapped_name = self.id_name[unmapped_id]
                if unmapped_name:
                    g.write(unmapped_name + '\n')
        line = open(gene_name, "r").readlines()[1:]
        if not line:
            self.map_name_bool = False
            self.map_blast_run()
        else:
            self.map_name = self.add_tool('tool_lab.ppi.ppinetwork_map')
            self.map_name.set_options({
                "diff_exp_gene": gene_name,
                "species": self.option("species")
            })
            self.map_name.on('end', self.set_output, 'map')
            self.map_name.on('start', self.set_step, {'start': self.step.map})
            self.map_name.on('end', self.set_step, {'end': self.step.map})
            self.map_name.on('end', self.merge_id_name)
            self.map_name.run()

    def merge_id_name(self):
        map_dir = os.path.join(self.output_dir, 'map')
        self.mapped_id_table = os.path.join(self.map.output_dir, "seq_mapped.txt")
        mapped_name_table = os.path.join(self.map_name.output_dir, 'seq_mapped.txt')
        unmapped_seq = os.path.join(self.map.output_dir, 'unmapped_seq.txt')
        self.unmapped_seq_new = os.path.join(map_dir, 'unmapped_seq.txt')
        with open(mapped_name_table, 'r') as mn, open(self.mapped_id_table, 'a+') as mi:
            for line in mn.readlines()[1:]:
                if line:
                    name, string_id = line.strip().split('\t')
                    mi.write(self.name_id[name] + '\t' + string_id + '\n')
        os.mkdir(map_dir)
        target_file = os.path.join(map_dir, 'seq_mapped.txt')
        os.link(self.mapped_id_table, target_file)

        unmapped_db_id = os.path.join(self.map.output_dir, "unmapped_db.txt")
        unmapped_db_name = os.path.join(self.map_name.output_dir, "unmapped_db.txt")
        self.unmapped_seq_merge = os.path.join(map_dir, 'unmapped_db.txt')
        df_id = pd.read_table(unmapped_db_id, header=None, index_col=None, sep='\t', dtype='str')
        df_name = pd.read_table(unmapped_db_name, header=None, index_col=None, sep='\t', dtype='str')
        id_name = list(set(df_name[0].tolist()).intersection(set(df_id[0].tolist())))
        df_data = pd.DataFrame(id_name)
        df_data.to_csv(self.unmapped_seq_merge, header=False, index=False, sep='\t')

        df_id_seq = pd.read_table(unmapped_seq, header=None, index_col=None, sep='\t')
        df_name_seq = pd.read_table(mapped_name_table, header=0, index_col=None, sep='\t')
        unmapped_new = list(set(df_id_seq[0].tolist())-set([self.name_id[i] for i in df_name_seq['gene_id'].tolist()]))
        df_data1 = pd.DataFrame(unmapped_new)
        df_data1.to_csv(self.unmapped_seq_new, header=False, index=False, sep='\t')
        self.map_blast_run()

    def map_blast_run(self):
        if self.map_name_bool == True:
            mapped_table = os.path.join(self.map.output_dir, 'seq_mapped.txt')
            unmapped_seq = self.unmapped_seq_new
            unmapped_db = self.unmapped_seq_merge
        else:
            mapped_table = os.path.join(self.map.output_dir, "seq_mapped.txt")
            unmapped_seq = os.path.join(self.map.output_dir, "unmapped_seq.txt")
            unmapped_db = os.path.join(self.map.output_dir, "unmapped_db.txt")
        line = open(self.seq, "r").readlines()
        if not line:
            self.map_blast = False
            self.ppinetwork_predict_run()
        else:
            self.map_blast = self.add_tool("tool_lab.ppi.ppinetwork_mapblast")
            self.map_blast.set_options({
                "mapped_table": mapped_table,
                "species": self.option("species"),
                "unmapped_seq": unmapped_seq,
                "unmapped_db": unmapped_db,
                "fa": self.seq,
                "blast": "blastx",
                "database": "customer_mode",
                "query_type": "nucl",
            })
            self.map_blast.on('end', self.set_output, 'map_blast')
            self.map_blast.on('start', self.set_step, {'start': self.step.map_blast})
            self.map_blast.on('end', self.set_step, {'end': self.step.map_blast})
            self.map_blast.on('end', self.ppinetwork_predict_run)
            self.map_blast.run()

    def ppinetwork_predict_run(self):
        if self.map_blast:
            diff_exp_mapped_table = os.path.join(self.work_dir, "PpinetworkMapblast/output/mapped_all.xls")
        else:
            output_xls = os.path.join(self.work_dir, "mapped_all.xls")
            gene_list = []
            with open(output_xls, 'w') as f:
                # 写入根据序列名比对的结果
                with open(os.path.join(self.map.output_dir, "seq_mapped.txt"), 'r') as f_id:
                    head = f_id.readline().strip()
                    f.write(head + "\tmethod\n")
                    for line in f_id.readlines():
                        if line.split("\t")[1] in gene_list:
                            pass
                        else:
                            gene_list.append(line.split("\t")[1])
                            f.write(line.strip() + "\tby_name\n")
            diff_exp_mapped_table = output_xls
        line = open(diff_exp_mapped_table, "r").readlines()[1:]
        if not line:
            self.set_error("基因集中的基因不能匹配到string数据库中id，请您确认基因id为Ensemble或者Entrez GeneID！", code="23700802")
        self.ppinetwork_predict.set_options({
            "diff_exp_mapped": diff_exp_mapped_table,
            "species": self.option("species"),
            'score': self.option('score'),
            "combine_score": self.option("combine_score")
        })
        self.ppinetwork_predict.on('end', self.set_output, 'ppinetwork_predict')
        self.ppinetwork_predict.on('start', self.set_step, {'start': self.step.ppinetwork_predict})
        self.ppinetwork_predict.on('end', self.set_step, {'end': self.step.ppinetwork_predict})
        self.ppinetwork_predict.on('end', self.ppinetwork_topology_run)
        self.ppinetwork_predict.run()

    def ppinetwork_topology_run(self):

        ppi_table = os.path.join(self.ppinetwork_predict.output_dir, "interaction_detail.txt")
        node_table = os.path.join(self.ppinetwork_predict.output_dir, "all_nodes.txt")
        # ppi_table = os.path.join(self.work_dir, "PpinetworkPredict/output/interaction_detail.txt")
        # node_table = os.path.join(self.work_dir, "PpinetworkPredict/output/all_nodes.txt")
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
        if event['data'] == 'ppinetwork_predict':
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
            self.get_file()


    def end(self):
        repaths = [
            [".", "", "蛋白质互作网络结果输出目录"],
            ["interaction.txt", "txt", "edges结果文件信息"],
            ["all_nodes.txt", "txt", "node结果信息"],
            ["network_stats.txt", "txt", "网络统计结果信息"],
            ["diff_exp_mapped.txt", "txt", "含有STRINGid的差异基因文件"],
            ["gene_interaction_network_centrality.txt", "txt", "PPI网络中心系数表"],
            ["gene_interaction_network_clustering.txt", "txt", "PPI网络节点聚类系数表"],
            ["gene_interaction_network_transitivity.txt", "txt", "PPI网络传递性"],
            ["gene_interaction_network_degree_distribution.txt", "txt", "PPI网络度分布表"],
            ["gene_interaction_network_by_cut.txt", "txt", "combined_score值约束后的PPI网络"],
            ["gene_interaction_network_node_degree.txt", "txt", "PPI网络节点度属性表"]
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
        data = {
            "id": "ppinetwork" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "tool_lab.ppinetwork_analysis",
            "instant": False,
            "options": dict(
                diff_exp_gene='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/PPI/PPI.txt',
                species=9606,
                combine_score=300
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
