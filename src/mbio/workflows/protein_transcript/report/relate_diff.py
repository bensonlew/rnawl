# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2018-12-11"

from biocluster.workflow import Workflow
import re
from collections import defaultdict
import unittest
from biocluster.wpm.client import *
import datetime
import copy
import os
from mbio.packages.itraq_and_tmt.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import glob


class RelateDiffWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelateDiffWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='diff_protein', type='string'),
            dict(name='gene_diff_id', type='string'),
            dict(name='gene_compare_group', type='string'),
            dict(name='protein_compare_group', type='string'),
            dict(name='diff_rna', type='string'),
            dict(name='relation_tab', type='string'),
            dict(name='fc', type='string'),
            dict(name="rna_type", type='string'),
            dict(name="seq", type="infile", format="itraq_and_tmt.common"),
            dict(name="rna_type", type='string', default="ref_rna_v2"),
            {"name": "species", "type": "int", "default": 9606},
            {"name": "combine_score", "type": "int", "default": 1000},
            dict(name="main_id", type='string'),
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

        self.protein_list = list()
        self.protein_up = list()
        self.protein_down = list()
        self.protein_no = list()
        self.transcript_list = list()
        self.transcript_up = list()
        self.transcript_down = list()
        self.transcript_no = list()
        self.proteinfc = dict()
        self.transcriptfc = dict()
        with open(self.option('diff_protein'), 'r') as pro, \
                open(self.option('diff_rna'), 'r') as trans, \
                open(self.option('relation_tab'), 'r') as rel_t:
            _ = pro.readline()
            _ = trans.readline()
            for line in pro.readlines():
                line = line.strip().split('\t')
                self.protein_list.append(line[0])
                self.proteinfc[line[0]] = line[1]
                if line[3] == 'yes':
                    if line[2].lower() == 'up':
                        self.protein_up.append(line[0])
                    elif line[2].lower() == 'down':
                        self.protein_down.append(line[0])
                    else:
                        self.protein_no.append(line[0])
                else:
                    self.protein_no.append(line[0])
            for line in trans.readlines():
                line = line.strip().split('\t')
                self.transcript_list.append(line[0])
                self.transcriptfc[line[0]] = line[1]
                if line[2].lower() == 'up' and line[3].lower() == 'yes':
                # if line[2].lower() == 'up':
                    self.transcript_up.append(line[0])
                elif line[2].lower() == 'down' and line[3].lower() == 'yes':
                # elif line[2].lower() == 'down':
                    self.transcript_down.append(line[0])
                else:
                    self.transcript_no.append(line[0])

            for line in rel_t:
                line = line.strip().split('\t')
                if line[0] == 'related' and len(line) >1:
                    self.related_list = line[1].split(';')
        self.protein_diff = list(set(self.protein_up + self.protein_down))
        self.transcript_diff = list(set(self.transcript_up + self.transcript_down))
        self.logger.info(len(self.transcript_diff))
        self.logger.info(len(self.related_list))
        self.logger.info(len(self.protein_diff))
        self.copy_protein_diff = copy.copy(self.protein_diff)
        self.remove_trans = list()
        self.remove_trans_diff = list()
        self.diff_relation_list = list()
        for protein in self.copy_protein_diff:
            for rela in self.related_list:
                if rela.split('|')[0] == protein:
                    gene = rela.split('|')[1]
                    if gene in self.transcript_diff or gene in self.remove_trans_diff:
                        self.diff_relation_list.append(rela)
                        self.protein_diff.remove(protein)
                        try:
                            self.transcript_diff.remove(gene)
                            self.remove_trans_diff.append(gene)
                        except:
                            self.remove_trans_diff.append(self.transcript_list.pop())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/7_relate/02_relate_exp/01_relate_diff')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(RelateDiffWorkflow, self).send_log(data)

    def run(self):
        if len(self.diff_relation_list) < 20: # by zoujiaxun 20200923: 工作流的run函数既有self.start_listener又有super(RelateDiffWorkflow, self).run()
            self.start_listener()
            self.fire("start")
            self.protein_gene()
            self.nine_blocks()
            self.diff_venn()
            self.run_ppi()
        else:
            self.protein_gene()
            self.nine_blocks()
            self.diff_venn()
            self.run_ppi()
            super(RelateDiffWorkflow, self).run()

    def set_db(self):
        # 导入ppi的信息，直接使用了geneset那边的方法
        if len(self.diff_relation_list) >= 20:
            interaction_path = self.ppinetwork.output_dir + '/ppinetwork_predict/interaction_detail.txt'
            all_nodes_path = self.ppinetwork.output_dir + '/ppinetwork_predict/all_nodes.txt'
            p = open(interaction_path)
            n = open(all_nodes_path)
            if len(p.readlines()) > 1 and len(n.readlines()) > 1:
                api_ppinetwork = self.api.api('protein_transcript.p2g_diff_ppi')
                all_nodes_path = self.ppinetwork.output_dir + '/ppinetwork_predict/all_nodes.txt'  # 画图节点属性文件
                with open(all_nodes_path, 'r') as node_r:
                    node_header = node_r.readline()
                    node_info = node_r.readlines()
                with open(all_nodes_path, 'w') as node_w:
                    node_w.write(node_header.strip() + '\tprotein_fc\ttranscript_fc\tprotein_regulate\ttranscript_regulate\n')
                    for line in node_info:
                        tmp = line.strip().split('\t')
                        if len(tmp) > 2:
                            for rel in self.diff_relation_list:
                                if rel.split('|')[0] == tmp[3]:
                                    node_w.write(line.strip() + '\t' + '\t'.join(self.up_or_down(rel)) + '\n')
                                    break

                os.link(all_nodes_path, self.output_dir + '/all_nodes.txt')
                interaction_path = self.ppinetwork.output_dir + '/ppinetwork_predict/interaction_detail.txt'  # 画图的边文件
                os.link(interaction_path, self.output_dir + '/interaction_detail.txt')
                network_stats_path = self.ppinetwork.output_dir + '/ppinetwork_predict/network_stats.txt'  # 网络全局属性统计
                os.link(network_stats_path, self.output_dir + '/network_stats.txt')
                try:
                    network_centrality_path = self.ppinetwork.output_dir + '/ppinetwork_topology/protein_interaction_network_centrality.txt'
                    os.link(network_centrality_path, self.output_dir + '/protein_interaction_network_centrality.txt')
                    degree_distribution_path = self.ppinetwork.output_dir + '/ppinetwork_topology/protein_interaction_network_degree_distribution.txt'
                    os.link(degree_distribution_path, self.output_dir + '/protein_interaction_network_degree_distribution.txt')
                    network_node_degree_path = self.ppinetwork.output_dir + '/ppinetwork_topology/protein_interaction_network_node_degree.txt'
                    os.link(network_node_degree_path, self.output_dir + '/protein_interaction_network_node_degree.txt')
                except:
                    pass
                if not os.path.isfile(all_nodes_path):
                    self.set_error("找不到报告文件:%s", variables=(all_nodes_path), code="12501501")
                if not os.path.isfile(interaction_path):
                    self.set_error("找不到报告文件:%s", variables=(interaction_path), code="12501502")
                if not os.path.isfile(network_stats_path):
                    self.set_error("找不到报告文件:%s", variables=(network_stats_path), code="12501503")
                # if not os.path.isfile(network_centrality_path):
                #     self.set_error("找不到报告文件:%s", variables=(network_centrality_path), code="12501504")
                # if not os.path.isfile(degree_distribution_path):
                #     self.set_error("找不到报告文件:%s", variables=(degree_distribution_path), code="12501505")
                # if not os.path.isfile(network_node_degree_path):
                #     self.set_error("找不到报告文件:%s", variables=(network_node_degree_path), code="12501506")

                print('start insert')
                api_ppinetwork.add_node_table(file_path=all_nodes_path, table_id=self.option("main_id"))  # 节点的属性文件（画网络图用）
                api_ppinetwork.add_edge_table(file_path=interaction_path, table_id=self.option("main_id"))  # 边信息
                api_ppinetwork.add_network_attributes(file2_path=network_stats_path,
                                                      table_id=self.option("main_id"))  # 网络全局属性

                # 节点的聚类与degree，画折线图
                if os.path.isfile(network_centrality_path):
                    api_ppinetwork.add_network_centrality(file_path=network_centrality_path, file1_path=all_nodes_path,
                                                          table_id=self.option("main_id"))  # 中心信息
                if os.path.isfile(degree_distribution_path):
                    api_ppinetwork.add_degree_distribution(file_path=degree_distribution_path, table_id=self.option("main_id"))
                print('end insert')
            else:
                pass
        # 导入diff，画diff——venn
        relation = self.api.api("protein_transcript.protein_transcript")
        params = dict(
            rna_type=self.option('rna_type'),
            species=self.option('species'),
        )
        relation.add_relation(main_id = self.option('main_id'), params = params, relation_file= self.output_dir + '/diff_relationship_list', table = 'sg_p2g_diff', detail_id= 'p2g_diff_id')


    #     导入九宫格等
        rel_diff = self.api.api("protein_transcript.p2g_diff")
        rel_diff.add_diff_venn(main_id=self.option('main_id'), venn_file=self.output_dir + '/diff_venn.xls')
        rel_diff.add_diff_nine(main_id=self.option('main_id'), nine_file=self.output_dir + '/nine_blocks.xls', fc_file = self.option('fc'))
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        protein_transcript_diff_allvenn = os.path.join(self.output_dir, "diff_relationship_list")
        protein_transcript_diff_relatevenn = os.path.join(self.output_dir, "diff_venn.xls")
        protein_transcript_diff_relatenine = os.path.join(self.output_dir, "nine_blocks.xls")
        protein_transcript_diff_relatenine_relateparams = os.path.join(self.work_dir, "fc")     #转录组的上下调差异倍数；蛋白组的上下调差异倍数；
        # protein_transcript_diff_relatenine_relateparams = "1.2;0.83;2.0;2.0"     #蛋白组的上下调差异倍数；转录组的上下调差异倍数；

        flag = 0
        if os.path.exists(protein_transcript_diff_allvenn):
            chart.chart_protein_transcript_diff_allvenn(protein_transcript_diff_allvenn)
            flag = 1

        if os.path.exists(protein_transcript_diff_relatevenn):
            chart.chart_protein_transcript_diff_relatevenn(protein_transcript_diff_relatevenn)
            flag = 1

        if os.path.exists(protein_transcript_diff_relatenine) and os.path.exists(protein_transcript_diff_relatenine_relateparams):
            protein_transcript_diff_relatenine_relateparams_str = ";".join([i.strip().split('\t')[1] for i in open(protein_transcript_diff_relatenine_relateparams).read().strip().split('\n')])
            chart.chart_protein_transcript_diff_relatenine(protein_transcript_diff_relatenine,protein_transcript_diff_relatenine_relateparams_str)
            flag = 1
        if flag:
            chart.to_pdf()

            # move pdf to result dir
            for ori,des in [["protein_transcript_diff_allvenn.venn2.pdf","venn_all.pdf"],\
            ["protein_transcript_diff_relatevenn.venn.pdf","venn.pdf"],\
            ["protein_transcript_diff_relatenine.nine.pdf","nineblocks.pdf"]\
            ]:
                pdf_file = os.path.join(self.work_dir, ori)
                if os.path.exists(pdf_file):
                    os.link(pdf_file, os.path.join(self.output_dir, des))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["7_relate", "", "蛋白组与转录组关联分析",0],
            ["7_relate/02_relate_exp", "", "表达量信息", 0],
            ["7_relate/02_relate_exp/01_relate_diff", "", "差异表达分析", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "联合分析样本数据关联结果"],
            ["./diff_relationship_list", "txt", "联合分析差异基因与蛋白关联结果表"],
            ["./nine_blocks.xls", "txt", "联合分析九宫格分析结果表"],
            ["./ppinetwork_map", "", "基因id比对到string数据库文件目录"],
            ["./ppinetwork_predict", "", "蛋白质相互作用组预测文件目录"],
            ["./ppinetwork_topology", "", "蛋白质互作网络拓扑属性文件目录"],
            ["./diff_exp_mapped.txt", "txt", "含有STRINGid的差异基因列表"],
            ["./all_nodes.txt", "txt", "PPI网络节点信息列表"],
            ["./network_stats.txt", "txt", "PPI网络统计结果表"],
            ["./interaction.txt", "txt", "PPI网络边信息列表"],
            ["./interaction_detail.txt", "txt", "PPI网络边信息列表"],
            ["./gene_protein.txt", "txt", "基因id与蛋白质对应表"],
            ["./protein_interaction_network_centrality.txt", "txt", "PPI网络中心系数表"],
            ["./protein_interaction_network_clustering.txt", "txt", "PPI网络节点聚类系数表"],
            ["./protein_interaction_network_transitivity.txt", "txt", "PPI网络传递性"],
            ["./protein_interaction_network_degree_distribution.txt", "txt", "PPI网络度分布表"],
            ["./protein_interaction_network_by_cut.txt", "txt", "根据综合值筛选得到的PPI网络"],
            ["./protein_interaction_network_node_degree.txt", "txt", "PPI网络节点度属性表"],
            ["./venn_all.pdf", "", "总体venn图"],
            ["./venn.pdf", "", "关联数据venn图"],
            ["./nineblocks.pdf", "", "关联数据九宫图"],
            ["./PPI.pdf", "", "关联数据PPI网络图"],
        ])
        super(RelateDiffWorkflow, self).end()

    def protein_gene(self):
        # for protein in self.copy_protein_diff:
        #     for rela in self.related_list:
        #         if rela.split('|')[0] == protein:
        #             gene = rela.split('|')[1]
        #             if gene in self.transcript_diff or gene in self.remove_trans_diff:
        #                 self.diff_relation_list.append(rela)
        #                 self.protein_diff.remove(protein)
        #                 try:
        #                     self.transcript_diff.remove(gene)
        #                     self.remove_trans_diff.append(gene)
        #                 except:
        #                     self.remove_trans_diff.append(self.transcript_list.pop())
        with open(self.output_dir + '/diff_relationship_list', 'w') as rel_w:
            rel_w.write('related' + '\t' + ';'.join(self.diff_relation_list) + '\n')
            rel_w.write('failed_proteins' + '\t' + ';'.join(self.protein_diff) + '\n')
            rel_w.write('failed_transcripts' + '\t' + ';'.join(self.transcript_diff) + '\n')

    def up_or_down(self, rela):
        p, g = rela.split('|')
        p_fc = self.proteinfc[p]
        g_fc = self.transcriptfc[g]
        if p in self.protein_up:
            relg_p = 'up'
        if p in self.protein_down:
            relg_p = 'down'
        if p in self.protein_no:
            relg_p = 'no change'
        if g in self.transcript_up:
            relg_g = 'up'
        if g in self.transcript_down:
            relg_g = 'down'
        if g in self.transcript_no:
            relg_g = 'no change'
        return p_fc, g_fc, relg_p, relg_g

    def nine_blocks(self):

        relations = list()
        proteins = copy.copy(self.protein_list)
        for protein in proteins:
            for rela in self.related_list:
                if rela.split('|')[0] == protein:
                    gene = rela.split('|')[1]
                    if gene in self.transcript_list:
                        relations.append(rela)
        with open(self.output_dir + '/nine_blocks.xls', 'w') as nine_w:
            nine_w.write('rela_id\tprotein_fc\ttranscript_fc\tprotein_regulate\ttranscript_regulate\trela_id_protein\trela_id_gene\n')
            for rel in relations:
                nine_w.write(rel + '\t' + '\t'.join(self.up_or_down(rel)) + '\t' + rel.split('|')[0] + '\t' + rel.split('|')[1] +'\n')

    def diff_venn(self):
        protein_up = list()
        protein_down = list()
        transcript_up = list()
        transcript_down = list()
        for rel in self.diff_relation_list:
            protein = rel.split('|')[0]
            transcript = rel.split('|')[1]
            if protein in self.protein_up:
                protein_up.append(rel)
            if protein in self.protein_down:
                protein_down.append(rel)
            if transcript in self.transcript_up:
                transcript_up.append(rel)
            if transcript in self.transcript_down:
                transcript_down.append(rel)
        with open(self.output_dir + '/diff_venn.xls', 'w') as diff_v:
            diff_v.write('protein_up' + '\t' + ';'.join(protein_up) + '\n')
            diff_v.write('protein_down' + '\t' + ';'.join(protein_down) + '\n')
            diff_v.write('transcript_up' + '\t' + ';'.join(transcript_up) + '\n')
            diff_v.write('transcript_down' + '\t' + ';'.join(transcript_down) + '\n')

    def run_ppi(self):
        if len(self.diff_relation_list) < 20:
            self.logger.info('关联到的diff蛋白只有%s个，跳过ppi的运行'%str(len(self.diff_relation_list)))
            self.set_db()
            return
        protein_list = self.work_dir + '/protein.list'
        with open(protein_list, 'w') as p_w:
            p_w.write("accession_id" + "\n")
            for rel in self.diff_relation_list:
                p_w.write(rel.split('|')[0]+'\n')
        self.ppinetwork = self.add_module('itraq_and_tmt.ppinetwork_analysis')
        options = {
            'diff_exp_gene': protein_list,
            'species': self.option('species'),
            'seq': self.option('seq').prop['path'],
            'combine_score': self.option('combine_score')
        }
        self.ppinetwork.set_options(options)
        self.ppinetwork.on('end', self.set_db)
        self.ppinetwork.run()

class TestFunction(unittest.TestCase):
    def test(self):
        # worker = worker_client()
        # id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        from mbio.workflows.protein_transcript.report.relate_diff import RelateDiffWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            "id": "protein_transcript_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "workflow",
            "name": "protein_transcript.report.relate_diff",
            "options": {
                "gene_compare_group": "Control|NaF",
                "rna_type": "ref_rna_v2",
                "seq": "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/MJ20190929137-陈克平陈亮-家蚕中肠-TMT-12个/protein.fa",
                "diff_protein": "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/MJ20190929137-陈克平陈亮-家蚕中肠-TMT-12个/diff_protein",
                "update_info": "{\"5f6ad970ffec603c79a2ebd7\": \"sg_p2g_diff\"}",
                "main_id": "5f6ad970ffec603c79a2ebd7",
                "gene_diff_id": "5ecd3c74ec02cc43aa59b37e",
                "fc": "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/MJ20190929137-陈克平陈亮-家蚕中肠-TMT-12个/fc",
                "relation_tab": "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/MJ20190929137-陈克平陈亮-家蚕中肠-TMT-12个/relation_tab",
                "combine_score": "300",
                "diff_rna": "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/MJ20190929137-陈克平陈亮-家蚕中肠-TMT-12个/diff_rna",
                "protein_compare_group": "NaF|DZ",
                "species": "394"
            }
        }
        wsheet = Sheet(data=data)
        wf =RelateDiffWorkflow(wsheet)
        wf.sheet.id = 'diff'
        wf.sheet.project_sn = 'diff'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


    # def test(self):
    #     import random
    #     from mbio.workflows.single import SingleWorkflow
    #     from biocluster.wsheet import Sheet
    #     data = {
    #         "id": "protein_transcript",
    #         "type": "workflow",
    #         "name": "protein_transcript.protein_transcript",
    #         "options": dict(
    #             rna_type="ref_rna_v1",
    #             gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/biomart/Homo_sapiens.GRCh37.biomart",
    #             pep="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/cds/Homo_sapiens.GRCh37.pep.fa",
    #             protein_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript/blast_testdata/protein.list",
    #             transcript_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript/blast_testdata/gene.list",
    #             protein_faa="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript/blast_testdata/DB.fasta",
    #             )
    #     }
    #     # data['options']['method'] = 'rsem'
    #     # wsheet = Sheet(data=data)
    #     # wf = SingleWorkflow(wsheet)
    #     # wf.run()
    #     # #
    #     # data['id'] += '1'
    #     # data['options']['method'] = 'salmon'
    #     # wsheet = Sheet(data=data)
    #     # wf = SingleWorkflow(wsheet)
    #     # wf.run()
    #     #
    #     data['id'] += '_fyt'
    #     wsheet = Sheet(data=data)
    #     wf = SingleWorkflow(wsheet)
    #     wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
