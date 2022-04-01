#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/21 16:27
@file    : gsea.py
"""
import csv
import glob
import json
import os
import time
import traceback
import unittest
from collections import OrderedDict, defaultdict
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM
import pandas as pd
import numpy as np
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import cairosvg

from mbio.packages.ref_rna_v2.geneset.gsea_exp_matrix import GseaExpMatrix
from multiprocessing import Pool


class GseaAgent(Agent):

    def __init__(self, parent):
        super(GseaAgent, self).__init__(parent)
        self.step.add_steps("gsea")
        self.add_option([
            # -----------------------  必须参数  -----------------------
            dict(name="geneset_source", type="string"),
            # 表达矩阵
            dict(name="matrix", type="infile", format="ref_rna_v2.ref_common"),
            # 表型，分组信息
            dict(name="g1", type="string"),
            dict(name="g2", type="string"),
            dict(name="group", type="string", required=True),
            # 基因集
            dict(name="gmx", type="infile", format="ref_rna_v2.geneset_gmt"),
            # id转换相关设置，如果不需要转换，则设为false
            dict(name="collapse", type="string", default='false'),
            dict(name="gene_detail", type="infile", format="ref_rna_v2.ref_common"),

            # -----------------------  基础参数  -----------------------
            dict(name="metric", type="string", default='Signal2Noise'),
            dict(name="scoring_scheme", type="string", default='weighted'),
            dict(name="sort", type="string", default='real'),
            dict(name="order", type="string", default='descending'),
            dict(name="set_max", type="int", default=500),
            dict(name="set_min", type="int", default=1),

            # ---------------  高级参数[建议默认，除非了解算法]  ---------------
            dict(name="mode", type="string", default='Max_probe'),
            dict(name="norm", type="string", default='meandiv'),
            dict(name="nperm", type="int", default=1000),
            dict(name="permute", type="string", default='phenotype'),
            dict(name="rnd_type", type="string", default='no_balance'),
            dict(name="rpt_label", type="string", default='my_analysis'),
            dict(name="create_gcts", type="string", default='false'),
            dict(name="create_svgs", type="string", default='true'),
            dict(name="include_only_symbols", type="string", default='true'),
            dict(name="make_sets", type="string", default='true'),
            dict(name="median", type="string", default='false'),
            dict(name="num", type="int", default=100),
            dict(name="plot_top_x", type="int", default=20),
            dict(name="rnd_seed", type="string", default='timestamp'),
            dict(name="save_rnd_lists", type="string", default='false'),
            dict(name="zip_report", type="string", default='false'),
        ])
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.gsea.start()
        self.step.update()

    def stepfinish(self):
        self.step.gsea.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        for name in ('matrix', 'gmx'):
            if not self.option(name).is_set:
                raise OptionError("缺少 %s 输入文件" , variables=( name), code="33710602")

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 6
        self._memory = '30G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"]
        #     ])
        # result_dir.add_regexp_rules([
        #     [r"go_enrich_.*", "xls", "go富集结果文件"],
        #     [r"go_lineage.*", "png", "go富集有向无环图"],
        #     [r"go_lineage.*", "pdf", "go富集有向无环图"],
        #     [r"go_lineage.*", "svg", "go富集有向无环图"],
        #     ])
        super(GseaAgent, self).end()


def svg_convert(fs):
    '''
    修改图片转换为并行
    '''
    print 'fs type is {} {}'.format(type(fs), fs)
    (f, out_f) = fs
    try:
        drawing = svg2rlg(f)
        renderPDF.drawToFile(drawing, out_f)
    except:
        # 避免renderPDF 不能转换成功
        try:
            cairosvg.svg2pdf(url=f, write_to=out_f)
        except:
            pass
    return True


class GseaTool(Tool):
    """
    """

    def __init__(self, config):
        super(GseaTool, self).__init__(config)
        self.gsea_dir = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/GSEA')
        perl_dir = os.path.join(self.config.SOFTWARE_DIR, 'program/perl-5.24.0/bin')
        self.set_environ(PATH=perl_dir)
        self.perl = 'program/perl-5.24.0/bin/perl'
        self.tmp_dir = os.path.join(self.work_dir, 'tmp_dir')
        if os.path.exists(self.tmp_dir):
            os.system("rm -f {}".format(self.tmp_dir))
        if not os.path.isdir(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        self.data_outdir = None
        self.name2id = dict()

    def run_cmd(self, cmd_name, cmd, is_wait=True, shell=False, ignore_error=False):
        self.logger.debug(cmd_name + ': ' + cmd + '%s' % type(cmd_name))
        cmd_obj = self.add_command(str(cmd_name), cmd, shell=shell, ignore_error=ignore_error)
        if shell is True:
            cmd_obj.software_dir = ''
            cmd_obj._start_run_time = int(time.time())
        cmd_obj.run()
        if is_wait is True:
            self._check_stat(cmd_obj)
            return
        return cmd_obj

    def _check_stat(self, *cmd_objs):
        self.wait(*cmd_objs)
        for cmd_obj in cmd_objs:
            if cmd_obj.return_code == 0:
                self.logger.info("{} Finished successfully".format(cmd_obj.name))
            elif cmd_obj.return_code is None:
                self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_obj.name))
                cmd_obj.rerun()
                self.wait()
                if cmd_obj.return_code is 0:
                    self.logger.info("{} Finished successfully".format(cmd_obj.name))
                else:
                    self.set_error("%s Failed. >>>%s", variables=(cmd_obj.name, cmd_obj.cmd), code="33710603")
            else:
                self.set_error("运行失败没有基因集满足条件， 请重新选择基因集或调整其它参数", code="33710604")
                # 没有富集的结果也会出错
                pass
                # self.set_error("{} Failed. >>>{}".format(cmd_obj.name, cmd_obj.cmd))

    def create_cls(self, samples, out_path):
        print(self.option('group'))
        group_dic = json.loads(self.option('group'), object_pairs_hook=OrderedDict)
        # group_dic = {sub_v: k for k, v in group_dic.items() for sub_v in v}
        # group_dic = {sample: group for group, sample in self.option('group').csv_reader()}
        samples_markers = [self.option('g1') for i in group_dic[self.option('g1')]] + \
                          [self.option('g2') for i in group_dic[self.option('g2')]]
        with open(out_path, 'w') as out_handler:
            out_handler.write('%s 2 1\n' % len(samples))
            out_handler.write('# {} {}\n'.format(self.option('g1'), self.option('g2')))
            out_handler.write(' '.join(samples_markers) + '\n')

    def gene_detail(self):
        exp_pd =  self.get_gene_expmatrix()
        gene_list = list(exp_pd.index)
        tmp_dic = {}
        for dic in self.option('gene_detail').csv_dict_reader():
            if dic['gene_id'] in gene_list:
                tmp_dic[dic['gene_id']] = dic
                tmp_dic[dic['gene_name']] = dic
                tmp_dic[dic['gene_name'].upper()] = dic
            self.name2id.update({dic['gene_name']: dic['gene_id']})
            self.name2id.update({dic['gene_name'].upper(): dic['gene_id']})
        return tmp_dic
        '''
        tmp_dic = {}
        for dic in self.option('gene_detail').csv_dict_reader():
            tmp_dic[dic['gene_id']] = dic
            tmp_dic[dic['gene_name']] = dic
            self.name2id.update({dic['gene_name']: dic['gene_id']})
        return tmp_dic
        '''

    def gsea_analysis(self, exp_path, cls_path, gmt, chip=None):
        """
        java
        -cp /mnt/ilustre/centos7users/dna/.env/bin/gsea-3.0.jar
        -Xmx512m xtools.gsea.Gsea
        -res  /mnt/ilustre/centos7users/jiameng.li/learning/GSEA/test/output/GSEA.txt
        -cls  /mnt/ilustre/centos7users/jiameng.li/learning/GSEA/test/output/GSEA.cls#control_versus_test  # Phenotype Data Formats
        -gmx /mnt/ilustre/centos7users/jiameng.li/learning/GSEA/test/output/GSEA.generate.gmt  # Gene Set Database Formats
        -mode Max_probe
        -norm meandiv
        -nperm 1000
        -permute phenotype
        -rnd_type no_balance
        -scoring_scheme weighted
        -rpt_label my_analysis
        -metric Signal2Noise
        -sort real
        -order descending
        -collapse true  # 为true时提供chip，目的probe id to  id
        -chip /mnt/ilustre/centos7users/jiameng.li/learning/GSEA/test/output/GSEA.chip  # Microarray Chip Annotation Formats
        -create_gcts false
        -create_svgs false
        -include_only_symbols true
        -make_sets true
        -median false
        -num 100
        -plot_top_x 20
        -rnd_seed timestamp
        -save_rnd_lists false
        -set_max 500 # 设置基因集大小，超过排除
        -set_min 1 # 设置基因集大小，超过排除
        -zip_report false
        -out /mnt/ilustre/centos7users/jiameng.li/learning/GSEA/test/output
        -gui false  # 启动UI界面
        :return:
        """
        out_dir = os.path.join(self.work_dir, 'tmp_output')
        self.data_outdir = out_dir
        if os.path.exists(self.data_outdir):
            os.system("rm -r {}".format(self.data_outdir))
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        cmd_list = [
            '{java} '.format(java='program/sun_jdk1.8.0/bin/java'),
            '-cp {cp} -Xmx25600m xtools.gsea.Gsea '.format(cp=os.path.join(self.gsea_dir, 'gsea-3.0.jar')),
            '-res {exp} '.format(exp=exp_path),
            '-cls {cls}#{g1}_versus_{g2}'.format(cls=cls_path, g1=self.option('g1'), g2=self.option('g2')),
            '-gmx {gmx} '.format(gmx=gmt),
            '-mode {mode} '.format(mode=self.option('mode')),
            '-norm {norm} '.format(norm=self.option('norm')),
            '-nperm {nperm} '.format(nperm=self.option('nperm')),
            '-permute {permute} '.format(permute=self.option('permute')),
            '-rnd_type {rnd_type} '.format(rnd_type=self.option('rnd_type')),
            '-scoring_scheme {scoring_scheme} '.format(scoring_scheme=self.option('scoring_scheme')),
            '-collapse {collapse} '.format(collapse=self.option('collapse')),
            '-rpt_label {rpt_label} '.format(rpt_label=self.option('rpt_label')),
            '-metric {metric} '.format(metric=self.option('metric').replace("log2_Radio_of_Classes", "log2_Ratio_of_Classes")),
            '-sort {sort} '.format(sort=self.option('sort')),
            '-order {order} '.format(order=self.option('order')),
            '-create_gcts {create_gcts} '.format(create_gcts=self.option('create_gcts')),
            '-create_svgs {create_svgs} '.format(create_svgs=self.option('create_svgs')),
            '-include_only_symbols {include_only_symbols} '.format(
                include_only_symbols=self.option('include_only_symbols')),
            '-make_sets {make_sets} '.format(make_sets=self.option('make_sets')),
            '-median {median} '.format(median=self.option('median')),
            '-num {num} '.format(num=self.option('num')),
            '-plot_top_x {plot_top_x} '.format(plot_top_x=self.option('plot_top_x')),
            '-rnd_seed {rnd_seed} '.format(rnd_seed=self.option('rnd_seed')),
            '-save_rnd_lists {save_rnd_lists} '.format(save_rnd_lists=self.option('save_rnd_lists')),
            '-set_max {set_max} '.format(set_max=self.option('set_max')),
            '-set_min {set_min} '.format(set_min=self.option('set_min')),
            '-zip_report {zip_report} '.format(zip_report=self.option('zip_report')),
            '-out {out} '.format(out=os.path.relpath(out_dir)),
            '-gui {gui} '.format(gui='false')
        ]
        if chip is not None:
            cmd_list.append('-chip {chip} '.format(chip=chip))

        cmd = ' '.join(cmd_list)
        self.run_cmd('gsea_analysis', cmd, is_wait=True, shell=False, ignore_error=True)


    def get_gene_expmatrix(self):
        exp_df = pd.read_table(self.option("matrix").prop['path'], header=0, index_col=0)
        exp_df = np.log(exp_df + 1)
        return exp_df

    def genesets_detail(self, geneset_names, outdir, gene_detail, is_symbol):
        fields = (
            'geneset_name', 'gene_id', 'gene_name', 'description', 'rank_in_geneset', 'rank_metric_score', 'running_es',
            'core_enrichement')
        set_detail_out_file = os.path.join(self.output_dir, 'all_sets.detail')
        exp_detail_out_file = os.path.join(self.output_dir, 'all_exp.detail')
        exp_df = self.get_gene_expmatrix()

        fields_exp = ['geneset_name', 'gene_id', 'gene_name', 'description']
        fields_exp.extend(list(exp_df.columns))

        out_handler = open(set_detail_out_file, 'w')
        out_handler.write('\t'.join(fields) + '\n')
        exp_out_handler = open(exp_detail_out_file, 'w')
        exp_out_handler.write('\t'.join(fields_exp) + '\n')
        line_demo = '\t'.join('{%s}' % i for i in fields) + '\n'
        line_exp_demo = '\t'.join('{%s}' % i for i in fields_exp) + '\n'

        out_list = []
        leading_edge = defaultdict(list)
        for name in geneset_names:
            set_detail_path = os.path.join(outdir, name.upper() + '.xls')
            if not os.path.isfile(set_detail_path):
                self.logger.debug('基因集[%s] 结果文件不存在[%s]' % (name, set_detail_path))
                continue
            with open(set_detail_path) as in_handler:
                # NAME	PROBE	GENE SYMBOL	GENE_TITLE	RANK IN GENE LIST	RANK METRIC SCORE	RUNNING ES	CORE ENRICHMENT
                for dic in csv.DictReader(in_handler, delimiter='\t'):
                    if is_symbol:
                        gene_name = dic['PROBE']
                        if gene_name in gene_detail:
                            sub_dic = gene_detail[gene_name]
                            gene_id = sub_dic['gene_id']
                            description = sub_dic['description']
                        else:
                            gene_id = gene_name
                            description = "-"
                    else:
                        gene_id = dic['PROBE']
                        if gene_id in gene_detail:
                            sub_dic = gene_detail[gene_id]
                            gene_name = sub_dic['gene_name']
                            description = sub_dic['description']
                        else:
                            gene_name = "-"
                            description = "-"
                    rank_in_geneset = dic['RANK IN GENE LIST']
                    rank_metric_score = dic['RANK METRIC SCORE']
                    running_es = dic['RUNNING ES']
                    core_enrichement = dic['CORE ENRICHMENT']
                    out_detail = line_demo.format(
                        geneset_name=name,
                        gene_id=gene_id,
                        gene_name=gene_name,
                        description=description,
                        rank_in_geneset=rank_in_geneset,
                        rank_metric_score=rank_metric_score,
                        running_es=running_es,
                        core_enrichement=core_enrichement
                    )
                    '''
                    out_list.append(line_demo.format(
                        geneset_name=name,
                        gene_id=gene_id,
                        gene_name=gene_name,
                        description=description,
                        rank_in_geneset=rank_in_geneset,
                        rank_metric_score=rank_metric_score,
                        running_es=running_es,
                        core_enrichement=core_enrichement
                    ))
                    '''

                    exps = list(exp_df.loc[gene_id])
                    exp_line = "\t".join([
                        name,
                        gene_id,
                        gene_name,
                        description
                    ] + map(str, exps)) + "\n"
                    if core_enrichement == 'Yes':
                        if name in leading_edge:
                            leading_edge[name].append(gene_id)
                        else:
                            leading_edge[name] = [gene_id]
                        exp_out_handler.write(exp_line)
                        out_handler.write(out_detail)
                    '''
                    if len(out_list) >= 1000:
                        # out_handler.write(''.join(out_list))
                        out_list = []

        if out_list:
            out_handler.write(''.join(out_list))
        '''
        out_handler.close()
        exp_out_handler.close()
        return set_detail_out_file, leading_edge


    def get_geneset2des(self):
        geneset2des = dict()
        with open(self.option("gmx").prop['path'], 'r') as in_handler:
            geneset2des = {line.split("\t")[0].upper():line.split("\t")[1] for line in in_handler}
        return geneset2des


    def sets_stat(self, outdir, leading_edge):
        out_file = os.path.join(self.output_dir, 'gsea_report.xls')  # 基因集富集汇总表
        fields = ('gene_set_name', 'group', 'size', 'es', 'nes', 'p_value', 'fdr_q_value', 'fwer_p_value', 'rank_at_max',
                  'leading_edge_num', 'leading_edge_genes', 'description')
        line_demo = '\t'.join('{%s}' % i for i in fields) + '\n'

        geneset2des = self.get_geneset2des()
        self.logger.info("geneset2des{}".format(geneset2des))
        with open(out_file, 'w') as out_handler:
            geneset_names = []
            out_handler.write('\t'.join(fields) + '\n')
            for file_path in glob.glob(os.path.join(outdir, 'gsea_report_for*.xls')):
                group = os.path.basename(file_path).split("_")[3]
                # 示例 gsea_report_for_na_pos_1558946742505.xls 取出 na_pos
                name = '_'.join(os.path.basename(file_path).split('_')[3: -1])


                out_list = []
                with open(file_path, 'r') as in_handler:
                    # SIZE	ES	NES	NOM p-val	FDR q-val	FWER p-val	RANK AT MAX	LEADING EDGE
                    for dic in csv.DictReader(in_handler, delimiter='\t'):
                        gene_set_name = dic['NAME'].upper()
                        geneset_names.append(gene_set_name)
                        size = dic['SIZE']
                        es = dic['ES']
                        nes = dic['NES']
                        p_value = dic['NOM p-val']
                        fdr_q_value = dic['FDR q-val']
                        fwer_p_value = dic['FWER p-val']
                        rank_at_max = dic['RANK AT MAX']
                        des = geneset2des.get(gene_set_name, '')
                        leading_edge_num = len(leading_edge[gene_set_name])
                        leading_edge_genes = ";".join(leading_edge[gene_set_name])
                        if leading_edge_num != 0:
                            out_list.append(
                                line_demo.format(
                                    gene_set_name=gene_set_name,
                                    group = group,
                                    size=size,
                                    es=es,
                                    nes=nes,
                                    p_value=p_value,
                                    fdr_q_value=fdr_q_value,
                                    fwer_p_value=fwer_p_value,
                                    rank_at_max=rank_at_max,
                                    leading_edge_num=leading_edge_num,
                                    leading_edge_genes=leading_edge_genes,
                                    description=des
                                )
                            )
                        if len(out_list) >= 500:
                            out_handler.write(''.join(out_list))
                            out_list = []
                    if out_list:
                        out_handler.write(''.join(out_list))
            return geneset_names

    def unzip(self, files, flag=[0]):
        cmd = '/bin/gunzip {}'.format(' '.join(files))
        flag[0] += 1
        flag_num = flag[0]
        return self.run_cmd('gunzip_%s' % flag_num, cmd, is_wait=False, shell=True, ignore_error=True)

    def svg2png(self, last_files):
        f_list = list()
        # out_f_list = list()
        for f in last_files:
            # f 为***.svg.gz 此时已有解压后的文件
            svg_f = f

            f = f[: -3]
            f_png = f.rsplit('.', 1)[0] + b'.png'
            out_name, index = os.path.basename(f).rsplit('.', 1)[0].rsplit(b'_', 1)
            try:
                int(index)
            except ValueError:
                out_name = os.path.basename(f).rsplit('.', 1)[0]

            outpng_f = os.path.join(self.output_dir, out_name + b'.png')
            if os.path.exists(f_png):
                if os.path.exists(outpng_f):
                    os.remove(outpng_f)
                os.link(f_png, outpng_f)

            out_f = os.path.join(self.output_dir, out_name + b'.pdf')
            f_list.append((f, out_f))
            # out_f_list.append(out_f)
        thread = 6
        p = Pool(thread)
        print "f_list is {}".format(f_list)
        p.map(svg_convert, f_list)
        p.close()
        p.join()

    def pictures(self, outdir, geneset_names):
        # for name in geneset_names:
        file_paths = glob.glob(os.path.join(outdir, 'enplot_*.svg.gz'))
        start = 0
        end = 50
        last_files = None
        while True:
            files = file_paths[start: end]
            if not files:
                break
            run_obj = self.unzip(files)
            if last_files != None:
                self.svg2png(last_files)
            self._check_stat(run_obj)
            start = end
            end += 50
            last_files = files
        if last_files != None:
            self.svg2png(last_files)

    def solve_data(self, gene_detail):
        geneset_names = (item[0].upper() for item in self.option('gmx').csv_reader())
        outdir = os.path.join(self.data_outdir, os.listdir(self.data_outdir)[0])

        is_symbol = False
        if self.option('geneset_source') == 'msigdb':
            is_symbol = True
        _, leading_edge = self.genesets_detail(geneset_names, outdir, gene_detail, is_symbol)
        new_geneset_names = self.sets_stat(outdir, leading_edge)
        self.pictures(outdir, new_geneset_names)

    def run(self):
        super(GseaTool, self).run()
        # gmt_obj, gmt_path = self.create_gmt()
        # exp_obj, exp_path = self.create_matrix()
        # self.wait(gmt_obj, exp_obj)
        # 目的避开带有 '-' 的路径
        exp_path = self.option('matrix').path
        gene_detail = self.gene_detail()
        exp_path = os.path.join(self.tmp_dir, 'GSEA.txt')
        gsea_exp = GseaExpMatrix(self.option('matrix').path, gene_detail, os.path.join(self.tmp_dir, 'GSEA.txt'), self.option('geneset_source'))
        samples = gsea_exp.run()
        '''
        if self.option('geneset_source') == 'msigdb':
            samples = gsea_exp.run()
        else:
            samples = next(self.option('matrix').csv_reader())
            samples.remove('seq_id')
            self.logger.debug(str(samples) + '=============================')
        '''
        cls_path = os.path.join(self.tmp_dir, 'GSEA.cls')
        self.create_cls(samples, cls_path)
        exp_matrix = os.path.relpath(exp_path)
        cls_path = os.path.relpath(cls_path)
        gmt_path = os.path.relpath(self.option('gmx').path)

        self.gsea_analysis(exp_matrix, cls_path, gmt_path)
        self.solve_data(gene_detail)
        self.end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'
            """
            dict(name="matrix", type="infile", format="denovo_rna_v2.ref_common"),
            dict(name="group", type="infile", format="denovo_rna_v2.ref_common"),
            dict(name="g1", type="string"),
            dict(name="g2", type="string"),
            dict(name="gset", type="infile", format="denovo_rna_v2.ref_common"),
            dict(name="chip", type="infile", format="denovo_rna_v2.ref_common"),
            """
            data = {
                "id": "GSEA_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "ref_rna_v2.geneset.gsea",
                "instant": False,
                "options": dict(
                    geneset_source='msigdb',
                    matrix='/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/temp/out/rsem.gene.TMM.fpkm.matrix',
                    group=json.dumps({'control': ["A_1", "A_2", "A_3"], 'test': ["B_1", "B_2", "B_3"]}),
                    g1='control',
                    g2='test',
                    gmx='/mnt/ilustre/users/sanger-dev/workspace/20190522/Single_GSEA_7777_5494/Gsea/tmp_output/'
                        'my_analysis_svgs.Gsea_test_vs_ctrl.1558659706068/edb/gene_sets.gmt',
                    # list='/mnt/ilustre/centos7users/jiameng.li/learning/GSEA/test/set',
                    gene_detail='/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/temp/out/GSEA.chip',
                    create_svgs='true',
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
