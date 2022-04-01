# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import unittest
import re
import copy
import pandas as pd
import numpy as np
from collections import defaultdict, OrderedDict
from sklearn import preprocessing
from sklearn.neighbors import KernelDensity
from mako.template import Template
from scipy.stats import pearsonr, spearmanr, kendalltau


class ProGeneCorrAgent(Agent):
    def __init__(self, parent):
        super(ProGeneCorrAgent, self).__init__(parent)
        options = [
            {"name": "group", "type": "string"},  # 分组文件
            {"name": "use_group", "type": "string"},  # 分组样本计算方式
            {"name": "exp_deal", "type": "string"},  # 表达量处理方法
            {"name": "protein_matrix", "type": "string"},  # 蛋白表达量矩阵
            {"name": "rna_matrix", "type": "string"},  # rna表达量矩阵
            {"name": "relaset_list", "type": "string"},  # 关联集列表文件
            {"name": "corr_method", "type": "string", "default": "spearman"},  # 计算相关性的方法
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
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(ProGeneCorrAgent, self).end()


class ProGeneCorrTool(Tool):
    def __init__(self, config):
        super(ProGeneCorrTool, self).__init__(config)
        self._version = '1.0.1'
        if '/mnt/ilustre/users/isanger' in os.getcwd():
            self.r_path = 'program/R-3.3.3/bin/Rscript'
        else:
            self.r_path = 'program/R-3.3.1/bin/Rscript'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/lib')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/lib64')
        self.corr_way = self.option('corr_method')
        if self.corr_way == "pearson":
            self.corr_way = pearsonr
        if self.corr_way == "spearman":
            self.corr_way = spearmanr
        if self.corr_way == "kendall":
            self.corr_way = kendalltau

    def get_sample_dicts(self):
        self.group_dict = OrderedDict()
        self.p_r_samples = OrderedDict()
        self.p_r_samples['protein'] = list()
        self.p_r_samples['rna'] = list()
        with open(self.option('group'), 'r') as g_r:
            g_r.readline()
            for line in g_r:
                line = line.strip().split('\t')
                if line[1] not in self.group_dict:
                    self.group_dict[line[1]] = list()
                self.group_dict[line[1]].append(line[0])
                self.group_dict[line[1]].append(line[0].replace('_rna', '_gene'))
                if line[1].endswith('_rna'):
                    self.p_r_samples['rna'].append(line[0].replace('_rna', '_gene'))
                else:
                    self.p_r_samples['protein'].append(line[0])

    def merge_matrix(self):
        def get_rel_col(inlist, rel_list, n):
            rel_copy = copy.copy(rel_list)
            outlist = list()
            for i in inlist:
                for rel in rel_copy:
                    if i == rel.split('|')[n]:
                        outlist.append(rel)
                        rel_copy.remove(rel)
                        break
            return outlist
        with open(self.option('relaset_list'),'r') as rel_r:
            rel=rel_r.readlines()
            proteins = [x.strip().split('|')[0] for x in rel]
            transcripts = [x.strip().split('|')[1] for x in rel]
            rels = [x.strip() for x in rel]
        self.merge_out = os.path.join(self.output_dir,'merge_matrix')
        protein_matrix = pd.read_table(self.option('protein_matrix'), index_col=0, dtype={0: str})
        protein_matrix = protein_matrix.loc[proteins, :]
        rna_matrix = pd.read_table(self.option('rna_matrix'), index_col=0)
        rna_matrix = rna_matrix.loc[list(set(transcripts)), :]
        rename_protein = OrderedDict()
        rename_rna = OrderedDict()
        for col in protein_matrix.columns:
            rename_protein[col] = col + '_protein'
        protein_matrix.rename(columns=rename_protein, inplace=True)
        for col in rna_matrix.columns:
            rename_rna[col] = col + '_gene'
        rna_matrix.rename(columns=rename_rna, inplace=True)
        merge_matrix = pd.DataFrame(columns=protein_matrix.columns.tolist() + rna_matrix.columns.tolist())
        for rel in rels:
            rel_dict = protein_matrix.loc[rel.split('|')[0], ].to_dict()
            rel_dict.update(rna_matrix.loc[rel.split('|')[1],].to_dict())
            rel_dict['rela_id'] = rel
            merge_matrix = merge_matrix.append(rel_dict, ignore_index=True)
        merge_matrix = merge_matrix.set_index('rela_id')
        #
        #
        # protein_matrix=pd.read_table(self.option('protein_matrix'),index_col=0, dtype={0:str})
        # protein_matrix=protein_matrix.loc[proteins,:]
        # rna_matrix = pd.read_table(self.option('rna_matrix'), index_col=0)
        # rna_matrix = rna_matrix.loc[transcripts, :]
        # protein_matrix['rela_id'] = get_rel_col(protein_matrix.index.tolist(), rels, 0)
        # rna_matrix['rela_id'] = get_rel_col(rna_matrix.index.tolist(), protein_matrix['rela_id'].tolist(),1)
        # protein_matrix = protein_matrix.set_index('rela_id')
        # rna_matrix = rna_matrix.set_index('rela_id')
        # rename_protein = OrderedDict()
        # rename_rna=OrderedDict()
        # for col in protein_matrix.columns:
        #     rename_protein[col] = col + '_protein'
        # protein_matrix.rename(columns=rename_protein, inplace=True)
        # for col in rna_matrix.columns:
        #     rename_rna[col] = col + '_gene'
        # rna_matrix.rename(columns=rename_rna, inplace=True)
        # merge_matrix = protein_matrix.join(rna_matrix, how='inner')
        # tmp_rna = list()
        # tmp_protein = list()
        # tmp_rela = list()
        # for i in merge_matrix.index:
        #     rna,protein = i.split('|')
        #     if rna not in tmp_rna and protein not in tmp_protein:
        #         tmp_protein.append(protein)
        #         tmp_rna.append(rna)
        #         tmp_rela.append(i)
        # merge_matrix = merge_matrix.loc[tmp_rela,:]
        merge_matrix.to_csv(self.merge_out, sep='\t', header=True, index=True)
        return merge_matrix

    def produce_two_matrix(self, matrix, scale = False):
        if self.group_dict is not None:
            if self.option('use_group') == 'mean':
                group_exp = list()
                for g in self.group_dict:
                    g_exp = matrix.loc[:, self.group_dict[g]].mean(axis=1)
                    g_exp.name = g.replace('_rna', '_gene')
                    group_exp.append(g_exp)
                all_exp_pd = pd.concat(group_exp, axis=1)
            elif self.option('use_group') == 'median':
                group_exp = list()
                for g in self.group_dict:
                    g_exp = matrix.loc[:, self.group_dict[g]].median(axis=1)
                    g_exp.name = g.replace('_rna', '_gene')
                    group_exp.append(g_exp)
                all_exp_pd = pd.concat(group_exp, axis=1)
            else:
                all_exp_pd = matrix
            if self.option('exp_deal') == 'log2':
                all_exp_pd = np.log2(all_exp_pd + 0.0001)
            elif self.option('exp_deal') == 'log10':
                all_exp_pd = np.log10(all_exp_pd + 0.0001)
            else:
                pass
            scaler = preprocessing.StandardScaler(with_std=False)
            scaler.fit(all_exp_pd.T)
            all_exp_pd1 = pd.DataFrame(scaler.transform(all_exp_pd.T)).T
            all_exp_pd1.columns = all_exp_pd.columns
            all_exp_pd1.insert(0, 'seq_id', all_exp_pd.index.tolist())
            all_exp_pd1 = all_exp_pd1.set_index('seq_id')
            all_exp_pd1.to_csv(os.path.join(self.output_dir, 'scaled_group_matrix.xls'), sep='\t', header=True, index=True)

        if self.p_r_samples:
            p_r_exp = list()
            p_exp = matrix.loc[:, self.p_r_samples['protein']].mean(axis=1)
            p_exp.name = 'protein'
            p_r_exp.append(p_exp)
            r_exp = matrix.loc[:, self.p_r_samples['rna']].mean(axis=1)
            r_exp.name = 'rna'
            p_r_exp.append(r_exp)
            p_r_exp = pd.concat(p_r_exp, axis=1)
            if self.option('exp_deal') == 'log2':
                p_r_exp = np.log2(p_r_exp + 0.0001)
            elif self.option('exp_deal') == 'log10':
                p_r_exp = np.log10(p_r_exp + 0.0001)
            else:
                pass
            if scale:
                scaler = preprocessing.StandardScaler(with_std=False)
                scaler.fit(p_r_exp.T)
                p_r_exp1 = pd.DataFrame(scaler.transform(p_r_exp.T)).T
                p_r_exp1.columns = p_r_exp.columns
                p_r_exp1.insert(0, 'seq_id', p_r_exp.index.tolist())
                p_r_exp1 = p_r_exp1.set_index('seq_id')
            else:
                p_r_exp1 = p_r_exp
            density = KernelDensity(bandwidth=0.01)
            p_r_exp1 = p_r_exp1.sort_values(by="protein")
            density.fit(p_r_exp1[["protein"]])
            p_r_exp1["protein density"] = density.score_samples(p_r_exp1[["protein"]])
            density = KernelDensity(bandwidth=0.01)
            p_r_exp1 = p_r_exp1.sort_values(by="rna")
            density.fit(p_r_exp1[["rna"]])
            p_r_exp1["rna density"] = density.score_samples(p_r_exp1[["rna"]])
            p_r_exp1["density"] = p_r_exp1["protein density"].apply(np.exp) + p_r_exp1["rna density"].apply(np.exp)
            p_r_exp1 = p_r_exp1.loc[:, ['protein', 'rna', 'density']]
            # p_r_exp1 = p_r_exp1.sort_values(by=p_r_exp1.columns.tolist())
            # density.fit(p_r_exp1)
            # p_r_exp1['density'] = density.score_samples(p_r_exp1)
            p_r_exp1.to_csv(os.path.join(self.output_dir, 'protein_rna_density_matrix.xls'), sep='\t', header=True,
                               index=True)
            corrtest = self.corr_way(p_r_exp1['protein'],
                                     p_r_exp1['rna'])
            with open(os.path.join(self.output_dir, 'corr_info'), 'w') as cor_w:
                cor_w.write('corr' + '\t' + str(corrtest[0]) + '\n')
                cor_w.write('p_value' + '\t' + str(corrtest[1]) + '\n')

    def run_mulit_corr(self):
        r_cmd = '''library(corrgram)

panel.raters <- function (x, y, corr = NULL, col.regions, cor.method, ...) {
  if (!is.null(corr)) 
    return()
  plot.xy(xy.coords(x, y), type = "p", col = 'blue', pch = 16,...)
  # lines(spline(x,y), col = 'brown4')
  abline(lm(y ~ x), col = 'brown4')
  box(col = "darkgray")
}

panel.text <- function (x, y, corr = NULL, col.regions, cor.method, digits = 2, 
    cex.cor, ...) 
{
    auto <- missing(cex.cor)
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    if (!is.null(corr)) {
        est <- corr
        absest <- formatC(abs(est), digits = digits, format = "f")
        est <- formatC(est, digits = digits, format = "f")
        if (auto) 
            cex.cor <- 0.7/strwidth(absest)
        text(0.5, 0.6, est, cex = cex.cor)
    }
    else {
        if (sum(complete.cases(x, y)) < 4) {
            warning("Need at least 4 complete cases for cor.test()")
        }
        else {
            results <- cor.test(x, y, alternative = "two.sided", 
                method = cor.method)
            cor_value <- round(cor(x,y,method = cor.method),4)
            est <- cor_value
            absest <- formatC(abs(est), digits = digits, format = "f")
            est <- formatC(est, digits = digits, format = "f")
            if (auto) 
                cex.cor <- 0.6/strwidth(absest)
            text(0.5, 0.6, est, cex = cex.cor)
            p <- results$p.value
            p <- formatC(p, digits = 3, format = "f")
            if (auto) 
                cex.cor <- 0.4/strwidth(p)
            text(0.5, 0.3, p, cex = cex.cor)
            box(col = "darkgray")
        }
    }
}

merge<-read.table('${merge_file}', header=T, sep='\t')
pdf('${pdf_file}', w=3*ncol(merge), h=3*ncol(merge))
corrgram(merge[,-1], diag.panel=panel.density, lower.panel=panel.raters, upper.panel=panel.text, cex.labels=2,cor.method = '${cor_method}')
dev.off()
svg('${svg_file}')
corrgram(merge[,-1], diag.panel=panel.density, lower.panel=panel.raters, upper.panel=panel.text, cex.labels=0.5,cor.method = '${cor_method}')
dev.off()
png('${png_file}', type="cairo", width=100*ncol(merge), height=100*ncol(merge))
corrgram(merge[,-1], diag.panel=panel.density, lower.panel=panel.raters, upper.panel=panel.text, cex.labels=1,cor.method = '${cor_method}')
dev.off()
cor_result<-cor(merge[,-1],method = '${cor_method}')
write.table(cor_result,file='${corr_file}',sep='\t')
'''

        f = Template(r_cmd)
        corr_file = os.path.join(self.output_dir, 'corr_result.xls')
        r_cmd = f.render(merge_file=os.path.join(self.output_dir, 'scaled_group_matrix.xls'),
                         pdf_file=os.path.join(self.output_dir, 'multily_corr.pdf'),
                         svg_file=os.path.join(self.output_dir, 'multily_corr.svg'),
                         png_file=os.path.join(self.output_dir, 'multily_corr.png'),
                         cor_method=self.option('corr_method'),
                         corr_file=corr_file
                             )
        r_path = os.path.join(self.work_dir, 'multily_corr.r')
        with open(r_path, 'w') as cor_r:
            cor_r.write(r_cmd)

        cmd = self.r_path + ' ' + r_path
        cmd_name = 'multily_corr'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>> %s", variables=(cmd_name, cmd))
        else:
            self.set_error("%s Failed. >>> %s", variables=(cmd_name, cmd))
        with open(corr_file) as corr_r:
            c=corr_r.read()
        with open(corr_file, 'w') as corr_w:
            corr_w.write('\t' + c.replace('"', ''))

    def run(self):
        super(ProGeneCorrTool, self).run()
        self.get_sample_dicts()
        matrix = self.merge_matrix()
        self.produce_two_matrix(matrix)
        self.run_mulit_corr()
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
            "name": "protein_transcript.extract_relation",
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
