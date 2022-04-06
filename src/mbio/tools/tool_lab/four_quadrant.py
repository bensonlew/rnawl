# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import os
import pandas as pd
import glob
import fitz

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class FourQuadrantAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(FourQuadrantAgent, self).__init__(parent)
        options = [
            {'name': 'dataframe1', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'dataframe2', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'id_map', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'fc1', 'type': 'float', 'default': 2},
            {'name': 'fc2', 'type': 'float', 'default': 1.2},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(FourQuadrantAgent, self).end()


class FourQuadrantTool(Tool):
    def __init__(self, config):
        super(FourQuadrantTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        # self.gcc = software_dir + '/gcc/5.1.0/bin'
        # self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=os.path.join(software_dir, 'program/Python/bin'))
        self.program = {
            'python': 'miniconda2/bin/python',
            'rscript': '/program/R-3.3.1/bin/Rscript'
        }
        self.script = {
            'four_quadrant': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/four_quadrant.r')
        }
        self.file = {
            'log2fc_table': os.path.join(self.work_dir, 'log2fc.txt'),
            'result_table': os.path.join(self.output_dir, 'log2fc.txt'),
            'pdf': os.path.join(self.work_dir, 'four_quadrant.pdf'),
            'png': os.path.join(self.output_dir, 'four_quadrant.png')
        }

    def run(self):
        super(FourQuadrantTool, self).run()
        self.data_preparation()
        self.run_four_quadrant()
        self.set_output()
        self.end()

    def sig_test(self, i, fc='', pval=''):
        fc_threshold = float(self.option('fc{}'.format(i)))
        p_threshold = 0.05
        if fc != '' and pval != '':
            if fc > fc_threshold and pval < p_threshold:
                sig = 'DE' + i
            else:
                sig = 'NDE' + i
        elif fc != '':
            if abs(fc) > fc_threshold:
                sig = 'DE' + i
            else:
                sig = 'NDE' + i
        elif pval != '':
            if pval < p_threshold:
                sig = 'DE' + i
            else:
                sig = 'NDE' + i
        return sig

    def col_select(self, df, order_i):
        col_list = df.columns
        select_col = [col_list[0]]
        select_dict = dict()
        if len(col_list) == 1:
            col_list_l = 'log2fc'
        else:
            col_list_l = [str(i).lower() for i in col_list]
        item_dict = {'log2fc': 'logfc', 'pvalue': 'padjust'}
        for each in item_dict.keys():
            if col_list_l.count(each) > 0:
                item_i = col_list_l.index(each)
                select_col.append(col_list[item_i])
                select_dict[each] = col_list[item_i]
            elif col_list_l.count(item_dict.get(each)) > 0:
                item_i = col_list_l.index(item_dict.get(each))
                select_col.append(col_list[item_i])
                select_dict[each] = col_list[item_i]
            else:
                self.logger.info('Lack of column of {} in dataframe'.format(each))
        final_df = df[select_col]
        if len(select_dict.keys()) == 2:
            final_df['sig'] = final_df.apply(lambda x: self.sig_test(order_i, fc=x[select_dict.get('log2fc')], pval=x[select_dict.get('pvalue')]), axis=1)
            final_df = final_df[[select_col[0], select_dict.get('log2fc'), 'sig']]
        elif select_dict.keys()[0] == 'log2fc':
            final_df['sig'] = final_df.apply(lambda x: self.sig_test(order_i, fc=x[select_dict.get('log2fc')]), axis=1)
        elif select_dict.keys()[0] == 'pvalue':
            final_df['sig'] = final_df.apply(lambda x: self.sig_test(order_i, pval=x[select_dict.get('pvalue')]), axis=1)
        return final_df

    def data_preparation(self):
        df1 = pd.read_table(self.option('dataframe1').path, sep='\t', header=0)
        final_df1 = self.col_select(df1, '1')
        df2 = pd.read_table(self.option('dataframe2').path, sep='\t', header=0)
        final_df2 = self.col_select(df2, '2')
        if self.option('id_map').is_set:
            id_df = pd.read_table(self.option('id_map').prop['path'], sep='\t', header=0)
            df_for_scatter = id_df.merge(final_df1, left_on=id_df.columns[0], right_on=final_df1.columns[0])
            df_for_result = df_for_scatter.merge(final_df2, left_on=df_for_scatter.columns[1], right_on=final_df2.columns[0])
            df_for_scatter = df_for_result.iloc[:, [0, 2, 3, 4, 5]]
        else:
            df_for_result = final_df1.merge(final_df2, left_index=True, right_index=True)
            df_for_scatter = df_for_result.iloc[:, [0, 1, 2, 4, 5]]
        col_names = df_for_scatter.columns.tolist()
        df_for_scatter.rename(columns={col_names[0]: 'node', col_names[2]: 'sig1', col_names[4]: 'sig2'}, inplace=True)
        df_for_scatter.to_csv(self.file['log2fc_table'], sep='\t', header=True, index=False)
        col_names = df_for_result.columns.tolist()
        df_for_result.rename(columns={col_names[0]: 'df1_node', col_names[1]: 'df2_node', col_names[2]: 'df1_data', col_names[3]: 'sig1', col_names[4]: 'df2_data', col_names[5]: 'sig2'}, inplace=True)
        df_for_result.to_csv(self.file['result_table'], sep='\t', header=True, index=False)

    def run_four_quadrant(self):
        ylab = os.path.splitext(os.path.basename(self.option('dataframe1').path))[0] + '(LogFC)'
        xlab = os.path.splitext(os.path.basename(self.option('dataframe2').path))[0] + '(LogFC)'
        cmd = '{} {} '.format(self.program['rscript'], self.script['four_quadrant'])
        cmd += '-d {} '.format(self.file['log2fc_table'])
        cmd += '-x {} -y {}'.format(xlab, ylab)
        cmd_name = 'run_four_quadrant'
        command = self.add_command(cmd_name, cmd)
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
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")

    def pdf2png(self):
        if os.path.exists(self.file['png']):
            os.remove(self.file['png'])
        with fitz.open(self.file['pdf']) as pdf:
            mat = fitz.Matrix(2, 2).preRotate(int(0))
            pix = pdf[0].getPixmap(matrix=mat, alpha=False)
            pix.writePNG(self.file['png'])

    def set_output(self):
        self.logger.info("设置结果目录")
        old = glob.glob(os.path.join(self.work_dir, '*pdf'))[0]
        name = os.path.basename(old)
        link = os.path.join(self.output_dir, name)
        if os.path.exists(link):
            os.remove(link)
        os.link(old, link)
        if os.path.exists(self.file['pdf']):
            self.pdf2png()




