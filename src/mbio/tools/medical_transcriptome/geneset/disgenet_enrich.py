#!/usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import unittest
import os
import glob
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class DisgenetEnrichAgent(Agent):
    def __init__(self, parent):
        super(DisgenetEnrichAgent, self).__init__(parent)
        options = [
            {"name": "gene_list", "type": "infile", "format": "medical_transcriptome.gene_list"}, # EntrezID list
            {'name': "g2e_path", "type": "infile", "format": "medical_transcriptome.common"}, #gene2enrtez file
            {"name": "padjust_method", "type": "string", "default": "BH"}, # one of BY, BH, Bonferroni, Holm
            {"name": "dsi", "type": "string", "default": ""},
            {"name": "score", "type": "string", "default": ""},
            {"name": "dpi", "type": "string", "default": ""},
            {"name": "el", "type": "string", "default": ""},   # comma between terms
            {"name": "ei", "type": "string", "default": ""},
            {"name": "universe", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "pvalue_cutoff", "type": "float", "default": 1.0},
            {"name": "qvalue_cutoff", "type": "float", "default": 1.0},
            {"name": "min", "type": "int", "default": 1},  # minimum genes annotated for testing
            {"name": "max", "type": "int", "default": 500},
        ]
        self.add_option(options)
        self.step.add_steps("disgenet_enrich")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.disgenet_enrich.start()
        self.step.update()

    def step_finish(self):
        self.step.disgenet_enrich.finish()
        self.step.update()

    def check_options(self):
        if self.option('padjust_method').lower() not in ['by', 'bh', 'bonferroni', 'holm']:
            raise OptionError('多重检验校正的方法不在提供的范围内', code = "33706602")
        if not self.option("gene_list").is_set:
            raise OptionError("必须设置输入基因集。")
        if not self.option("g2e_path").is_set:
            raise OptionError("必须设置输入gene2entrez文件。")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"]
        # ])
        # result_dir.add_regexp_rules([
        #     [r"disgenet_enrichment.xls$", "xls", "DisGeNET富集分析结果"]
        # ])
        super(DisgenetEnrichAgent, self).end()


class DisgenetEnrichTool(Tool):
    def __init__(self, config):
        super(DisgenetEnrichTool, self).__init__(config)
        self._version = "v1.0"
        self.r_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/miniconda3/bin/Rscript"
        self.enricher_path = self.config.PACKAGE_DIR + "/medical_transcriptome/enricher.r"
        self.db_path = self.config.SOFTWARE_DIR + "/database/DisGeNET/7.0"
        self.g2e_dict = dict()
        self.info_dict = dict()

    def run(self):
        super(DisgenetEnrichTool, self).run()
        self.generate_background()
        self.g2e_process()
        self.result_check()
        self.set_output()
        self.end()

    def generate_background(self):
        db_file = os.path.join(self.db_path, 'browser_source_summary_gda.tsv')
        if os.path.exists(db_file):
            disgenet = pd.read_table(db_file, header=0)
            disgenet = disgenet[disgenet.Type == "disease"]
            if self.option("dsi") == "":
                disgenet = disgenet
            else:
                disgenet = disgenet.drop(disgenet[disgenet.DSI_g < float(self.option("dsi"))].index).reset_index(
                    drop=True)
                disgenet = disgenet.reset_index(drop=True)
            if self.option("dpi") == "":
                disgenet = disgenet
            else:
                disgenet = disgenet.drop(disgenet[disgenet.DPI_g < float(self.option("dpi"))].index).reset_index(
                    drop=True)
            if self.option("score") == "":
                disgenet = disgenet.drop(disgenet[disgenet.Score_gda < 0.5].index).reset_index(drop=True)
            else:
                disgenet = disgenet.drop(disgenet[disgenet.Score_gda < float(self.option("score"))].index).reset_index(
                    drop=True)
            if self.option("ei") == "":
                disgenet = disgenet
            else:
                disgenet = disgenet.drop(disgenet[disgenet.EI_gda < float(self.option("ei"))].index).reset_index(
                    drop=True)
            if "all" in self.option("el").lower().split(","):
                disgenet = disgenet
            elif self.option("el") == "":
                disgenet = disgenet[disgenet.EL_gda == "definitive"]
            else:
                els = self.option("el").lower().split(",")
                # els = ['' if i == 'None' else i for i in els]
                disgenet = disgenet[disgenet.EL_gda.isin(els)]
            disgenet.to_csv(self.work_dir + "/disgenet_background.txt", sep="\t", header=True, index=False)
            term2gene = disgenet.loc[:, ["Disease_id", "Gene_id"]]
            term2gene = term2gene.rename(columns={"Gene_id": "Entrez_ID"})  # obtain term to gene file
            term2gene.to_csv(self.work_dir + "/term2gene.txt", sep="\t", header=True, index=False)
            term2name = disgenet.loc[:, ["Disease_id", "Disease"]]  # obtain term to name file
            # term2name = term2name.drop_duplicates(keep="last")
            term2name.to_csv(self.work_dir + "/term2name.txt", sep="\t", header=True, index=False)
            self.run_enricher()

    def run_enricher(self):
        term2gene = self.work_dir + "/term2gene.txt"
        # term2name = self.work_dir + "/term2name.txt"
        cmd = '{} {}'.format(self.r_path, self.enricher_path)
        cmd += ' -g {}'.format(self.option('gene_list').prop["path"])
        cmd += ' -p {}'.format(self.option('pvalue_cutoff'))
        cmd += ' -d {}'.format(self.option('padjust_method'))
        cmd += ' -t {}'.format(term2gene)
        if self.option('universe').is_set:
            cmd += ' -u {}'.format(self.option('universe').prop["path"])
        # cmd += ' -r {}'.format(term2name)
        cmd += ' -m {}'.format(self.option('min'))
        cmd += ' -n {}'.format(self.option('max'))
        cmd += ' -q {}'.format(self.option('qvalue_cutoff'))

        cmd_name = 'run_enrichment'
        command = self.add_command(cmd_name, cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704307")

    def result_check(self):
        term2name = self.work_dir + "/term2name.txt"
        result_path = self.work_dir + "/enrichment.txt"
        if not os.path.exists(term2name):
            self.set_error("注释文件%s未找到", variables=(term2name,),)
        if not os.path.exists(result_path):
            self.set_error("DisGeNET富集结果文件%s未找到", variables=(result_path,),)
        if os.path.getsize(result_path) > 1:
            t2n_dict = dict()
            with open(term2name, 'r') as t2n:
                t2n.readline()
                for line in t2n:
                    terms = line.strip().split("\t")
                    if terms[0] not in t2n_dict.keys():
                        t2n_dict[terms[0]] = terms[1]
                    else:
                        continue
            result = pd.read_table(result_path, header=0, sep="\t")
            result.rename(columns={'geneID': 'EntrezID'}, inplace=True)
            result['Description'] = result.apply(lambda x: t2n_dict.get(x['ID']), axis=1)
            result['GeneName'] = result.apply(lambda x: self.get_info(term=x['EntrezID'], t='name'), axis=1)
            result['GeneID'] = result.apply(lambda x: self.get_info(term=x['EntrezID'], t='id'), axis=1)
            result.to_csv(self.work_dir + "/enrichment_result.txt", sep="\t", header=True, index=False)
        if os.path.getsize(result_path) <= 1:
            os.remove(result_path)
            self.logger.warn("输入的基因集没有富集到任何结果")
        self.detailed_result()

    def g2e_process(self):
        with open(self.option('g2e_path').prop['path'], 'r') as g2e:
            g2e.readline()
            for line in g2e:
                each = line.strip('\n').split("\t")
                if each[1] not in self.g2e_dict.keys():
                    self.g2e_dict[each[1]] = each[0]
                if each[1] not in self.info_dict.keys():
                    self.info_dict[each[1]] = dict()
                    if each[2]:
                        self.info_dict[each[1]]['name'] = each[2]
                    self.info_dict[each[1]]['id'] = each[0]

    def get_info(self, term, t):
        terms = str(term).split('/')
        info = list()
        for each in terms:
            detail = self.info_dict.get(each)
            if detail:
                info.append(detail.get(t))
        if info:
            info_str = '/'.join(info)
        else:
            info_str = ''
        return info_str


    def detailed_result(self):
        result_path = self.work_dir + "/enrichment.txt"
        detail_result = self.work_dir + '/DisGeNET_enrich_detail.xls'
        g2id_dict = dict()
        if os.path.exists(result_path):
            with open(result_path, 'r') as r:
                r.readline()
                for line in r:
                    terms = line.strip().split("\t")
                    entrez_l = terms[-2].split('/')
                    for each in entrez_l:
                        if each not in g2id_dict.keys():
                            g2id_dict[each] = list()
                        g2id_dict[each].append(terms[0])
            with open(detail_result, 'w') as d:
                d.write('Query\tDisGeNET_ID\n')
                for key, value in g2id_dict.items():
                    gene = self.g2e_dict.get(key)
                    disgenet_id = ';'.join(value)
                    d.write('{}\t{}\n'.format(gene, disgenet_id))

    def set_output(self):
        enrich = glob.glob(self.work_dir + '/enrichment_result.txt')
        detail = glob.glob(self.work_dir + '/DisGeNET_enrich_detail.xls')
        files = enrich+detail
        for each in files:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "disgenet_enrich_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "medical_transcriptome.disgenet_enrich",
            "instant": False,
            "options": dict(
                gene_list="/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/gene_test.list",
                padjust_method="BH",
                dsi="0.1",
                dpi="0.1",
                score="0.1",
                el="Definitive,Limited",
                min=1,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)