# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import unittest
import pandas as pd
import numpy as np
from collections import defaultdict, OrderedDict
from sklearn import preprocessing
from scipy.stats import pearsonr, spearmanr, kendalltau
from statsmodels.stats.multitest import multipletests
import json
from concurrent.futures import ThreadPoolExecutor
import itertools


class ExpressPrCorrAgent(Agent):
    def __init__(self, parent):
        super(ExpressPrCorrAgent, self).__init__(parent)
        options = [
            {"name": "group_rna", "type": "string"},  # 转录分组文件
            {"name": "group_protein", "type": "string"},  # 蛋白分组文件
            {"name": "protein_matrix", "type": "string"},  # 蛋白表达量矩阵
            {"name": "rna_matrix", "type": "string"},  # rna表达量矩阵
            {"name": "geneset_list", "type": "string"},  # 基因集列表文件
            {"name": "proteinset_list", "type": "string"},  # 蛋白集列表文件
            {"name": "corr_method", "type": "string", "default": "spearman"},  # 计算相关性的方法
            dict(name="cor_cutoff", type='float', default=0.8),
            dict(name="pvalue_cutoff", type='float', default=0.05),
            dict(name="qvalue_cutoff", type='float', default=0.05),
            dict(name="padjust_way", type='string', default="fdr_bh"),
            dict(name="sig_type", type='int', default=1),
            ]
        self.add_option(options)
        self.step.add_steps('relation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.relation.start()
        self.step.update()

    def step_end(self):
        self.step.relation.finish()
        self.step.update()

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = 8
        self._memory = '20G'

    def end(self):
        super(ExpressPrCorrAgent, self).end()


class ExpressPrCorrTool(Tool):
    def __init__(self, config):
        super(ExpressPrCorrTool, self).__init__(config)
        self.cor_cutoff = self.option('cor_cutoff')
        self.qvalue_cutoff = self.option('qvalue_cutoff')
        self.pvalue_cutoff = self.option('pvalue_cutoff')
        self.padjust_way = self.option('padjust_way')
        self.corr_way = self.option('corr_method')
        if self.corr_way == "pearson":
            self.corr_way = pearsonr
        if self.corr_way == "spearman":
            self.corr_way = spearmanr
        if self.corr_way == "kendall":
            self.corr_way = kendalltau
        self.sig_type = self.option('sig_type')
        with open(self.option('proteinset_list'), 'r') as p_l:
            self.protein_list = p_l.read().split('\n')
            while '' in self.protein_list:
                self.protein_list.remove('')
        with open(self.option('geneset_list'), 'r') as r_l:
            self.rna_list = r_l.read().split('\n')
            while '' in self.rna_list:
                self.rna_list.remove('')
        self.sample_dict_protein = OrderedDict()
        with open(self.option('group_protein'), 'r') as g_p:
            g_p.readline()
            for line in g_p:
                line = line.strip().split('\t')
                if line[1] not in self.sample_dict_protein:
                    self.sample_dict_protein[line[1]] = list()
                self.sample_dict_protein[line[1]].append(line[0])
        self.sample_dict_rna = OrderedDict()
        with open(self.option('group_rna'), 'r') as g_r:
            g_r.readline()
            for line in g_r:
                line = line.strip().split('\t')
                if line[1] not in self.sample_dict_rna:
                    self.sample_dict_rna[line[1]] = list()
                self.sample_dict_rna[line[1]].append(line[0])

    def deal_pr_matrix(self):
        protein_matrix = pd.read_table(self.option('protein_matrix'), dtype={0:str}, index_col=0)
        protein_matrix = protein_matrix.loc[self.protein_list, :]
        # group_exp = list()
        # for g in self.sample_dict_protein:
        #     g_exp = protein_matrix.loc[:, self.sample_dict_protein[g]].mean(axis=1)
        #     g_exp.name = g
        #     group_exp.append(g_exp)
        # protein_matrix = pd.concat(group_exp, axis=1)
        protein_matrix = np.log2(protein_matrix + 0.0001)
        scaler = preprocessing.StandardScaler(with_std=False)
        scaler.fit(protein_matrix.T)
        protein_matrix1 = pd.DataFrame(scaler.transform(protein_matrix.T)).T
        protein_matrix1.columns = protein_matrix.columns
        protein_matrix1.insert(0, 'seq_id', protein_matrix.index.tolist())
        protein_matrix = protein_matrix1.set_index('seq_id')

        rna_matrix = pd.read_table(self.option('rna_matrix'), index_col=0)
        rna_matrix = rna_matrix.loc[self.rna_list, :]
        # group_exp = list()
        # for g in self.sample_dict_rna:
        #     g_exp = rna_matrix.loc[:, self.sample_dict_rna[g]].mean(axis=1)
        #     g_exp.name = g
        #     group_exp.append(g_exp)
        # rna_matrix = pd.concat(group_exp, axis=1)
        rna_matrix = np.log2(rna_matrix + 0.0001)
        scaler = preprocessing.StandardScaler(with_std=False)
        scaler.fit(rna_matrix.T)
        rna_matrix1 = pd.DataFrame(scaler.transform(rna_matrix.T)).T
        rna_matrix1.columns = rna_matrix.columns
        rna_matrix1.insert(0, 'seq_id', rna_matrix.index.tolist())
        rna_matrix = rna_matrix1.set_index('seq_id')
        return protein_matrix, rna_matrix

    def calculate_corr(self, protein_matrix, rna_matrix):
        def calculate(shape):
            i = shape[0]
            j = shape[1]
            corrtest = self.corr_way(rna_matrix[rna_matrix.columns[i]].tolist(), protein_matrix[protein_matrix.columns[j]].tolist())
            print(rna_matrix[rna_matrix.columns[i]])
            print(protein_matrix[protein_matrix.columns[j]])
            print(corrtest)
            coeffmat[i, j] = corrtest[0]
            pvalmat[i, j] = corrtest[1]
        protein_matrix = protein_matrix[(protein_matrix != 0).any(axis=1)]
        del protein_matrix.index.name
        protein_matrix = protein_matrix.transpose()
        rna_matrix = rna_matrix[(rna_matrix != 0).any(axis=1)]
        del rna_matrix.index.name
        rna_matrix = rna_matrix.transpose()
        coeffmat = np.zeros((rna_matrix.shape[1], protein_matrix.shape[1]))
        pvalmat = np.zeros((rna_matrix.shape[1], protein_matrix.shape[1]))
        with ThreadPoolExecutor(6) as pool:
            pool.map(calculate, list(itertools.product(range(rna_matrix.shape[1]), range(protein_matrix.shape[1]))))
        # for i in range(rna_matrix.shape[1]):
        #     for j in range(protein_matrix.shape[1]):
        #         corrtest = self.corr_way(rna_matrix[rna_matrix.columns[i]], protein_matrix[protein_matrix.columns[j]])
        #         coeffmat[i,j] = corrtest[0]
        #         pvalmat[i,j] = corrtest[1]
        dfcoeff = pd.DataFrame(coeffmat, columns=protein_matrix.columns, index=rna_matrix.columns)
        dfpvals = pd.DataFrame(pvalmat, columns=protein_matrix.columns, index=rna_matrix.columns)
        flat_pvalue = [item for sublist in dfpvals.values.tolist() for item in sublist]
        flat_adjust_pvalue = multipletests(flat_pvalue, alpha=0.05, method=self.padjust_way, is_sorted=False, returnsorted=False)
        adjust_p_list = flat_adjust_pvalue[1].tolist()
        size = protein_matrix.shape[1]
        nested_adjust_pvalue = [adjust_p_list[i:i+size] for i in range(0, len(adjust_p_list), size)]
        df_adjust = pd.DataFrame(nested_adjust_pvalue, index=rna_matrix.columns, columns=protein_matrix.columns)
        output_tmp = self.work_dir
        output = self.output_dir
        out_dir_file = os.path.join(output_tmp, "express_correlation_info_pre")
        with open(out_dir_file, "w") as f4:
            seq = ("gene_id", "protein_id", "cor", "p_value", "q_value")
            f4.write("\t".join(seq) + "\n")
            for i in range(rna_matrix.shape[1]):
                for j in range(protein_matrix.shape[1]):
                    seq1 = [dfcoeff.index[i], dfcoeff.columns[j], str(dfcoeff.loc[dfcoeff.index[i], dfcoeff.columns[j]])]
                    seq2 = str(dfpvals.loc[dfpvals.index[i], dfpvals.columns[j]])
                    seq3 = str(df_adjust.loc[df_adjust.index[i], df_adjust.columns[j]])
                    f4.write("\t".join(seq1) + "\t" + seq2 + "\t" + seq3 + "\n")
        df_info = pd.read_table(out_dir_file, header=0, sep="\t")
        out_dir_file1 = os.path.join(output, "express_correlation_info.xls")
        if not self.sig_type:
            df_selected_temp = df_info[((df_info["cor"] > self.cor_cutoff) | (df_info["cor"] < -self.cor_cutoff)) & (df_info["p_value"] < self.pvalue_cutoff)]
            df_selected_temp.to_csv(out_dir_file1, sep="\t", index=False)
            df_selected = pd.read_table(out_dir_file1, header=0, sep="\t")
        else:
            df_selected_temp = df_info[((df_info["cor"] > self.cor_cutoff) | (df_info["cor"] < -self.cor_cutoff)) & (df_info["q_value"] < self.qvalue_cutoff)]
            df_selected_temp.to_csv(out_dir_file1, sep="\t", index=False)
            df_selected = pd.read_table(out_dir_file1, header=0, sep="\t")

        df_source_pre = df_selected.loc[:, ['gene_id', 'protein_id', 'cor']]
        df_nodes = df_source_pre[['gene_id', 'protein_id']]
        list1 = df_nodes['gene_id'].tolist()
        list2 = df_nodes['protein_id'].tolist()
        whole_list_dict = OrderedDict()
        whole_list = list(OrderedDict.fromkeys(list1 + list2))
        group = ['gene' if i in list1 else 'protein' for i in whole_list]
        for i in range(len(whole_list)):
            whole_list_dict[whole_list[i]] = i
        df_source_column1 = [whole_list_dict[x] for x in list1]
        df_source_column2 = [whole_list_dict[x] for x in list2]
        df_source_column3 = df_source_pre['cor']
        df_links = pd.DataFrame(OrderedDict({'source': df_source_column1, 'target': df_source_column2, 'distance': df_source_column3}))
        new_list = list()
        # nodes_links_dict['links'] = df_links.to_dict('records')
        for i in df_links.to_dict('records'):
            new_dict = OrderedDict()
            for key in i.keys():
                if key != "distance":
                    value = int(i[key])
                    new_dict.update({key: value})
                else:
                    new_dict.update({"distance": i['distance']})
            new_list.append(new_dict)
            # new_dict.clear()
            del new_dict
        nodes_links_dict = OrderedDict()
        # df_nodes_new_pre = pd.DataFrame(OrderedDict({"id": whole_list, "group": group}))
        df_nodes_new_pre = pd.DataFrame(OrderedDict({"id": whole_list, "group": group}))
        nodes_links_dict['nodes'] = df_nodes_new_pre.to_dict('records')
        nodes_links_dict['links'] = new_list
        out_dir_file3 = os.path.join(output, "record.json")
        with open(out_dir_file3, "w") as f:
            json.dump(nodes_links_dict, f, separators=(',', ': '), indent=2)

    def run(self):
        super(ExpressPrCorrTool, self).run()
        protein_matrix, rna_matrix = self.deal_pr_matrix()
        self.calculate_corr(protein_matrix, rna_matrix)
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87'
        data = {
            "id": "Extract_relation_" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "protein_transcript_labelfree.extract_relation",
            "instant": False,
            "options": dict(
                pep = test_dir + "/" + "cds/Homo_sapiens.GRCh37.pep.fa",
                relation_file = test_dir + "/" + "biomart/Homo_sapiens.GRCh37.biomart",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
