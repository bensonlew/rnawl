# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2019-01-14"

from biocluster.workflow import Workflow
import unittest
from biocluster.wpm.client import *
import datetime
import re
import os
from biocluster.config import Config
from collections import OrderedDict, defaultdict
import copy
import json
from mbio.packages.dia_v3.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder

class PrKeggClassWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PrKeggClassWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='gene_kegg_table', type='string'),
            dict(name='gene_kegg_level_table', type='string'),
            dict(name='protein_kegg_table', type='string'),
            dict(name='protein_kegg_level_table', type='string'),
            dict(name='pr_list', type='string'),
            dict(name='add_info_rna', type='string'),
            dict(name='add_info_protein', type='string'),
            dict(name='protein_kegg_id', type='string'),
            dict(name='kegg_version', type='string', default='2017'),
            dict(name='task_id', type='string'),
            dict(name="rna_type", type='string', default="ref_rna_v2"),
            dict(name="main_id", type='string'),
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.kegg_class = self.add_tool("protein_transcript_dia_v3.kegg_class")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/7_relate/03_relate_anno/04_relate_kegg')
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
        super(PrKeggClassWorkflow, self).send_log(data)

    def merge_pro_gene_kegg(self):
        ko_genes_list = os.path.join(Config().SOFTWARE_DIR, 'database/KEGG/ko_genes.list')
        ko_gene = defaultdict(list)
        gene_ko = defaultdict(list)
        with open(ko_genes_list, 'r') as k_r:
            for line in k_r:
                line = line.strip().split('\t')
                ko_gene[line[0].split(':')[1]].append(line[1])
                gene_ko[line[1]].append(line[0].split(':')[1])
        # merge add_info
        map2pr = dict()
        # map2pr = defaultdict(dict)
        with open(self.option('add_info_protein'), 'r') as info_p:
            _ = info_p.readline()
            tmp = info_p.readline()
            species_protein = filter(str.isalpha, tmp.split('\t')[0])
            url_pre = tmp.split('\t')[1].split('?')[0]
            info_p.seek(0,0)
            _ = info_p.readline()
            for line in info_p:
                map, link = line.strip().split('\t')
                proteins = link.split('?')[1].split('/')[1:]
                proteins = [x.split('%09')[0] for x in proteins]
                tmp_dict = {x:'green' for x in proteins}
                if map not in map2pr:
                    map2pr[map] = tmp_dict
                else:
                    map2pr[map].update(tmp_dict)
        with open(self.option('add_info_rna'), 'r') as info_r:
            _ = info_r.readline()
            species_rna = filter(str.isalpha, info_r.readline().split('\t')[0])
            info_r.seek(0,0)
            _ = info_r.readline()
            for line in info_r:
                map, link = line.strip().split('\t')
                map = map.replace(species_rna, species_protein)
                genes_tmp = link.split('?')[1].split('/')[1:]
                genes_tmp = [x.split('%09')[0] for x in genes_tmp]
                genes = list()
                if species_protein != 'map':
                    for gene in genes_tmp:
                        if gene in ko_gene:
                            for i in ko_gene[gene]:
                                if i.startswith(species_protein+':'):
                                    genes.append(i)
                else:
                    genes = copy.copy(genes_tmp)
                tmp_dict = dict()
                # with open('/mnt/ilustre/users/sanger-dev/workspace/20190314/PrKeggClass_tsg_33529_4266_3986/tm','a') as t:
                #     t.write(';'.join(genes) +'\n')
                if map in map2pr:
                    for gene in genes:
                        if gene in map2pr[map]:
                            tmp_dict[gene] = 'orange'
                        else:
                            tmp_dict[gene] = 'pink'
                    map2pr[map].update(tmp_dict)
                else:
                    for gene in genes:
                        tmp_dict[gene] = 'pink'
                    map2pr[map] = tmp_dict

        col_dict = dict(
            pink = '#FFC0CB',
            green='#00CD00',
            orange='#FFA500'
        )
        self.add_info = os.path.join(self.work_dir, 'add_info_pr')
        # 这两个字典还会用来调节html.mark的导表
        self.map2pr = map2pr
        self.ko_gene = ko_gene
        self.species_protein = species_protein
        self.ko_txt = os.path.join(self.work_dir, 'pr_ko_txt')
        with open(self.add_info, 'w') as info_pr, open(self.ko_txt, 'w') as k_w:
            k_w.write('#KO\tbg\tfg\n')
            ko2bg = dict()
            ko_list = list()
            p_list = list()
            info_pr.write('pathway\thyperlink\n')
            for map in map2pr:
                info_str = map + '\t' + url_pre + '?' + map + '/'
                tmp = map2pr[map]
                tmp_list = ['{}%09{}'.format(pr, tmp[pr]) for pr in tmp]
                info_str += '/'.join(tmp_list)
                info_pr.write(info_str + '\n')
                for pr in tmp:
                    if pr in ko2bg and ko2bg[pr] == 'orange':
                        continue
                    else:
                        if species_protein != 'map':
                            prs = list()
                            # for i in ko_gene:
                            #     if pr in ko_gene[i]:
                            #         prs.append(i)
                            try:
                                prs = gene_ko[pr]
                            except:
                                pass
                            for p in prs:
                                if p not in p_list:
                                    if p not in ko2bg:
                                        ko2bg[p] = tmp[pr]
                                    else:
                                        if ko2bg[p] != 'orange':
                                            ko2bg[p] = tmp[pr]
                                    # k_w.write(p + '\t' + col_dict[tmp[pr]] + '\t' + 'NA' + '\n')
                                    p_list.append(p)
                        else:
                            if pr not in ko2bg:
                                ko2bg[pr] = tmp[pr]
                            else:
                                if ko2bg[pr] != 'orange':
                                    ko2bg[pr] = tmp[pr]

                        # ko_list.append(pr)
            for ko in ko2bg:
                k_w.write(ko + '\t' + col_dict[ko2bg[ko]] + '\t' + 'NA' + '\n')

        # merge kegg_table
        self.kegg_table = os.path.join(self.work_dir, 'pr_kegg_table.xls')
        with open(self.option('protein_kegg_table'), 'r') as pkr,\
            open(self.option('gene_kegg_table'), 'r') as rkr,\
            open(self.kegg_table, 'w') as prw:
            prw.write('#Query\tKO_ID\tKO_name\tHyperlink\tPaths\tKEGG_gene_id\n')
            _ = pkr.readline()
            for line in pkr:
                prw.write(line)
            _ = rkr.readline()
            for line in rkr:
                tmp = line.strip().split('\t')
                ko_id = tmp[1].split(';')
                kegg_gene_id = set()
                if species_protein != 'map':
                    for ko in ko_id:
                        if ko in ko_gene:
                            for i in ko_gene[ko]:
                                if i.startswith(species_protein):
                                    kegg_gene_id.add(i)
                if kegg_gene_id:
                    line = line.strip('\n') + ';'.join(kegg_gene_id) + '\n'
                prw.write(line)

        # merge kegg_level_table
        self.kegg_level_table = os.path.join(self.work_dir, 'pr_kegg_level_table.xls')
        with open(self.option('gene_kegg_level_table'), 'r') as rkr, \
                open(self.option('protein_kegg_level_table'), 'r') as pkr, \
                open(self.kegg_level_table, 'w') as prw:
            prw.write('Pathway_id\tgraph_id\tnumber_of_seqs\tpathway_definition\tfirst_category\tanno_type\thyperlink\tseq_list\tgraph_png_id\tsecond_category\n')
            map2level_info = dict()
            headers = pkr.readline().strip().split('\t')
            for line in pkr:
                line = line.strip('\n').split('\t')
                tmp_dict = {headers[n]:x for n,x in enumerate(line)}
                map2level_info[line[0]] = tmp_dict
            _ = rkr.readline()
            for line in rkr:
                line = line.strip('\n').split('\t')
                line[0] = line[0].replace(species_rna, species_protein)
                map = line[0]
                tmp_dict = {headers[n]: x for n, x in enumerate(line)}
                if map in map2level_info:
                    map2level_info[map]['seq_list'] = ';'.join(set(map2level_info[map]['seq_list'].split(';') + tmp_dict['seq_list'].split(';')))
                else:
                    tmp_dict['Pathway_id'] = map
                    map2level_info[map] = tmp_dict
            for map in map2level_info:
                info_str = url_pre + '?' + map + '/'
                tmp = map2pr[map]
                tmp_list = ['{}%09{}'.format(pr, tmp[pr]) for pr in tmp]
                info_str += '/'.join(tmp_list)
                map2level_info[map]['hyperlink'] = info_str
                tmp_str = '\t'.join(['{' + key + '}' for key in headers])
                prw.write(tmp_str.format(**map2level_info[map]) + '\n')

    def run_kegg_class(self):
        opts = {
            "proteinset_kegg": self.option("pr_list"),
            "kegg_table": self.kegg_table,
            "kegg_id": self.option("protein_kegg_id"),
            "background_links": self.add_info,
            "task_id": self.option("task_id"),
            "ko_txt": self.ko_txt,
            "gene_kegg_table": self.option("gene_kegg_table"),
            "kegg_version": self.option("kegg_version"),
        }
        self.kegg_class.set_options(opts)
        self.kegg_class.run()

    def run(self):
        self.kegg_class.on("end", self.set_db)
        self.merge_pro_gene_kegg()
        self.run_kegg_class()
        super(PrKeggClassWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_kegg = self.api.api('protein_transcript_dia_v3.relaset')

        kegg_stat = self.kegg_class.output_dir + '/kegg_stat.xls'
        pathway_class = self.kegg_level_table
        pathway_file = self.kegg_class.output_dir + '/pathways'

        self.logger.info("开始进行kegg_class的导表")
        api_kegg.add_kegg_regulate_pr(self.option("main_id"), self.option("pr_list"), kegg_stat, pathway_class, self.get_workflow_output_dir())
        api_kegg.add_kegg_regulate_pic(self.option("main_id"), self.kegg_level_table,
                                             pathway_file, self.map2pr, self.ko_gene, self.species_protein, self.ko_txt)

        self.end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if re.match(r'tsanger:', workflow_output):
            workflow_output = workflow_output.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'^\w+://\S+/.+$', workflow_output):
            workflow_output = workflow_output
        else:
            workflow_output = workflow_output.replace('sanger:', '/mnt/ilustre/data/')
        return workflow_output

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'
        # pt_keggclass_level = os.path.join(self.work_dir, "pr_kegg_level_table.xls")
        pt_proteinset_file = os.path.join(self.work_dir, "gene_protein.list")
        pt_kegg_stat_file = os.path.join(self.kegg_class.output_dir, "kegg_stat.xls")
        pt_gene_kegg_level_table_xls = os.path.join(self.work_dir, "pr_kegg_level_table.xls")

        if os.path.exists(pt_proteinset_file) and os.path.exists(pt_kegg_stat_file) and os.path.exists(
                pt_gene_kegg_level_table_xls):
            chart.chart_pt_keggclass(pt_proteinset_file, pt_kegg_stat_file, pt_gene_kegg_level_table_xls)
            chart.to_pdf()
            # move pdf to result dir
            pdf_file = os.path.join(self.work_dir, "pt_keggclass.barline.pdf")
            if os.path.exists(pdf_file):
                os.link(pdf_file, os.path.join(self.kegg_class.output_dir, "bar.pdf"))

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.kegg_class.output_dir)
        self.inter_dirs = [
            ["7_relate", "", "蛋白组与转录组关联分析", 0],
            ["7_relate/03_relate_anno", "", "功能注释信息", 0],
            ["7_relate/03_relate_anno/04_relate_kegg", "", "KEGG注释", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "基因集KEGG功能分类结果目录", 0, "220079"],
            ["pathways", " ", "KEGG分析结果通路表", 0, "220080"],
            ["kegg_analysis_of_anotate.xls", " ", "KEGG分类统计表", 0, "220081"],
            ["kegg_statistics.xls", " ", "KEGG分类统计结果表", 0, "220082"],
            ["bar.pdf", " ", "Pathway分类统计柱状图", 0],
        ])
        super(PrKeggClassWorkflow, self).end()

class TestFunction(unittest.TestCase):
    def test(self):
        worker = worker_client()
        id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        data = {
            "id": "protein_transcript_labelfree_" + id,
            "type": "workflow",
            "name": "protein_transcript_labelfree.protein_transcript_labelfree",
            "options": dict(
                rna_type="ref_rna_v1",
                gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/biomart/Homo_sapiens.GRCh37.biomart",
                pep="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/cds/Homo_sapiens.GRCh37.pep.fa",
                protein_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/protein.list",
                transcript_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/gene.list",
                protein_faa="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/DB.fasta",
                )
        }

        info = worker.add_task(data)
        print(info)


    # def test(self):
    #     import random
    #     from mbio.workflows.single import SingleWorkflow
    #     from biocluster.wsheet import Sheet
    #     data = {
    #         "id": "protein_transcript_labelfree",
    #         "type": "workflow",
    #         "name": "protein_transcript_labelfree.protein_transcript_labelfree",
    #         "options": dict(
    #             rna_type="ref_rna_v1",
    #             gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/biomart/Homo_sapiens.GRCh37.biomart",
    #             pep="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/cds/Homo_sapiens.GRCh37.pep.fa",
    #             protein_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/protein.list",
    #             transcript_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/gene.list",
    #             protein_faa="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/DB.fasta",
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
    unittest.main()
