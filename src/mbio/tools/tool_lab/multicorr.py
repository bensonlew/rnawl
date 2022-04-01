# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re
import pandas as pd
import scipy.stats as stats
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
# from matplotlib import pyplot as plt
# %matplotlib inline
import seaborn as sns



class MulticorrAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(MulticorrAgent, self).__init__(parent)
        options = [
            {"name": "infile_a", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "infile_b", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "infile_group", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "r_threshold", "type": "string", 'default': "0.8"},
            {"name": "p_threshold", "type": "string", 'default': "0.05"},
            {"name": "corr_method", "type": "string", 'default': "pearson"} # "pearson","spearman","kendall"
        ]
        self.add_option(options)
        self.step.add_steps('multicorr')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.multicorr.start()
        self.step.update()

    def step_end(self):
        self.step.multicorr.finish()
        self.step.update()

    def check_options(self):
        if not self.option('infile_a').is_set:
            raise OptionError('必须输入表达量文件。')
        if not self.option('infile_b').is_set:
            raise OptionError('必须输入表达量文件。')
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(MulticorrAgent, self).end()


class MulticorrTool(Tool):
    def __init__(self, config):
        super(MulticorrTool, self).__init__(config)

    def scatter_plot(self, cor_file, group_file, cor_value, p_value, plot_out_dir):
        df = pd.read_csv(cor_file, sep="\t", index_col=0)
        df_group = pd.read_csv(group_file, sep="\t", index_col=0)
        df = df.join(df_group, how='outer')
        # df

        sns.set(rc={'figure.figsize':(10,7)})
        print df
        print "__xuxi__"
        print df.columns[0],df.columns[1],df.columns[2]
        ax = sns.scatterplot(x=df.columns[0], y=df.columns[1], data=df, hue=df.columns[2])
        sns.regplot(data=df, x=df.columns[0], y=df.columns[1], scatter=False, ax=ax)
        matplotlib.pyplot.legend(loc='upper right')
        x_margin = (max(df.iloc[:,0])-min(df.iloc[:,0]))*0.2
        y_margin = (max(df.iloc[:,1])-min(df.iloc[:,1]))*0.2
        matplotlib.pyplot.xlim(min(df.iloc[:,0])-x_margin, max(df.iloc[:,0])+x_margin*1.5)
        matplotlib.pyplot.ylim(min(df.iloc[:,1])-y_margin, max(df.iloc[:,1])+y_margin)
        # matplotlib.pyplot.title('Attack&Defense of Pokemons')
        matplotlib.pyplot.text(min(df.iloc[:,0])-x_margin, max(df.iloc[:,1]+y_margin*0.5), "cor={}\npval={}".format(str(cor_value),str(p_value)), horizontalalignment='left', size='medium', color='black')
        matplotlib.pyplot.savefig(os.path.join(plot_out_dir,"{}.svg".format(os.path.basename(cor_file).split('.')[0])))
        matplotlib.pyplot.savefig(os.path.join(plot_out_dir,"{}.pdf".format(os.path.basename(cor_file).split('.')[0])))
        matplotlib.pyplot.clf()

    def run_multicorr(self):
        table_file1 = self.option('infile_a').prop['path']
        table_file2 = self.option('infile_b').prop['path']
        output_dir = self.output_dir
        r_threshold = float(self.option('r_threshold'))
        p_threshold = float(self.option('p_threshold'))
        # corr_method = "pearson"
        corr_method = self.option('corr_method')
        # corr_method = "kendall"

        self.logger.info("开始进行 组间相关性分析！")
        table1 = pd.read_table(table_file1, sep='\t', index_col=0)
        table2 = pd.read_table(table_file2, sep='\t', index_col=0)
        if list(table1.columns) == list(table2.columns):
            flag = 1
            good_cor = []
            for table1_index in table1.index:
                for table2_index in table2.index:
                    if corr_method == "pearson":
                        r,p = stats.pearsonr(table1.loc[table1_index,:], table2.loc[table2_index,:])
                    if corr_method == "spearman":
                        r,p = stats.spearmanr(table1.loc[table1_index,:], table2.loc[table2_index,:])  # 相关系数和P值
                    if corr_method == "kendall":
                        r,p = stats.kendalltau(table1.loc[table1_index,:], table2.loc[table2_index,:])  # 相关系数和P值
                    r = float('%.4f' % r)
                    p = float('%.4e' % p)
                    print(r,p)
                    if abs(r) >= r_threshold and p < p_threshold:
                        tmp_df = pd.concat([table1.loc[table1_index,:], table2.loc[table2_index,:]], axis=1)
                        tmp_df.columns = [tmp_df.columns[0]+"("+table1.index.name+")", tmp_df.columns[1]+"("+table2.index.name+")"]
                        tmp_df.to_csv(os.path.join(output_dir, "cor_"+str(flag)+".txt"), sep='\t')
                        # one_correlation = {"cor_index":"cor_"+str(flag), str(table1.index.name):str(table1_index), str(table2.index.name):str(table2_index), "P_value":p, "corr":r}
                        one_correlation = OrderedDict([("cor_index","cor_"+str(flag)),(str(table1.index.name),str(table1_index)),(str(table2.index.name),str(table2_index)),("P_value",p),("corr",r)])
                        good_cor.append(one_correlation)
                        # print('相关系数r为 = %6.3f，p值为 = %6.3f'%(r,p))
                        self.scatter_plot(os.path.join(output_dir, "cor_"+str(flag)+".txt"),\
                         self.option('infile_group').prop['path'], r, p, output_dir)
                        flag += 1
            if good_cor:
                good_cor_df = pd.DataFrame(good_cor)
                good_cor_df.to_csv(os.path.join(output_dir, "result_cor.txt"), sep='\t',index=False)
            else:
                with open(os.path.join(output_dir, "result_cor.txt"),'w') as w:
                    w.write("cor_index\t{}\t{}\tP_value\tcorr\t".format(str(table1.index.name), str(table2.index.name)))
            self.logger.info("完成组间相关性分析！")
        else:
            # print "ERROR: two files have different columns,can't do correlation analyse"
            self.set_error("ERROR: two files have different columns,can't do correlation analyse !")

    def run(self):
        super(MulticorrTool, self).run()
        self.run_multicorr()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "multicorr_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.multicorr",
            "instant": False,
            "options": dict(
                infile='/mnt/lustre/users/sanger-dev/sg-users/xuxi/multicorr_test_dir/test_data.txt',
                target_lab='species',
                show_apart="heat-tree"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)