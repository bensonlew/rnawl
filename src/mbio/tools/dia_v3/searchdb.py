# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd
import time
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import unittest


class SearchdbAgent(Agent):
    def __init__(self, parent):
        super(SearchdbAgent, self).__init__(parent)
        options = [
            {"name": "protein", "type": "infile", "format": "labelfree.common"},
            {"name": "ratio_exp", "type": "infile", "format": "labelfree.ratio_exp"},
            {"name": "protein_fasta", "type": "infile", "format": "labelfree.common"},
            {"name": "report_dia", "type": "infile", "format": "labelfree.common"},
        ]
        self.add_option(options)
        self.step.add_steps("Searchdb")
        self.on("start", self.step_start)
        self.on("end", self.step_end)


    def step_start(self):
        self.step.Searchdb.start()
        self.step.update()

    def step_end(self):
        self.step.Searchdb.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数

        :return:
        """
        if not self.option("protein").is_set:
            raise OptionError("必须设置输入文件:protein基本信息文件")

        if not self.option("ratio_exp").is_set:
            raise OptionError("必须设置输入文件:蛋白的scaled文件")

        if not self.option("protein_fasta").is_set:
            raise OptionError("必须设置输入文件:蛋白序列文件")

        if not self.option("report_dia").is_set:
            raise OptionError("必须设置输入文件:DIA report文件")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"],
        #     ["search_db.xls", " ", "搜库信息统计文件"],
        # ])
        super(SearchdbAgent, self).end()


class SearchdbTool(Tool):
    def __init__(self, config):
        super(SearchdbTool, self).__init__(config)

    def get_protein_char(self):
        """
        Using ProteinAnalysis to get protein characters
        """
        seq_dict = dict()
        pep_dict = dict()
        with open(self.option("protein_fasta").prop["path"], 'r') as fa:
            for line in fa:
                if line.startswith('>'):
                    title = line.strip().split(">")[1]
                    seq_dict[title] = ''
                else:
                    line = re.sub('[XUZB]', '', line.strip('*'))    # replace special characters
                    seq_dict[title] += line.strip().upper()
        report_selected = self.extract_data()
        for row in report_selected.itertuples():
            accession = row[1].split(";")
            for i in accession:
                if i not in pep_dict.keys():
                    pep_dict[i] = set()
                    pep = row[2].strip('\n')
                    pep_dict[i].add(pep)
                else:
                    pep = row[2].strip('\n')
                    pep_dict[i].add(pep)

        result = os.path.join(self.work_dir, "protein_analysis.txt")
        with open(result, 'w') as r:
            r.write('Accession\tpeptides\tlength\tmolecular_weight\tisoelectric_point\n')
            for key, value in seq_dict.items():
                p = ProteinAnalysis(value)
                if pep_dict.get(key):
                    pep_num = len(pep_dict.get(key))
                else:
                    pep_num = 0
                r.write('{}\t{}\t{}\t{}\t{}\n'.format(key, str(pep_num), str(p.length), str(format(float(p.molecular_weight())/1000, '.1f')), str(p.isoelectric_point())))

    def extract_data(self):
        report = pd.read_table(self.option("report_dia").prop["path"], header=0, sep="\t")
        report_selected = pd.DataFrame(index=report.index)
        colnames = report.columns
        previous = ['PG.ProteinAccessions', 'PEP.GroupingKey']
        latest = ['PG.ProteinGroups', 'PEP.StrippedSequence']
        col_list = latest
        for i in range(2):  # using these two cols to work out num of peptides
            if previous[i] in colnames:
                col_list = previous
                break
            elif latest[i] in colnames:
                col_list = latest
                break
        for col in col_list:  # using these two cols to work out num of peptides
            try:
                report_selected[col] = report[col]
            except:
                self.logger.info('文件里缺少%s这一列，请检查' % col)
                report_selected[col] = '_'
        return report_selected

    def searchdb(self):
        protein = pd.read_table(self.option("protein").prop["path"], header=0, sep="\t")
        ratio_exp = pd.read_table(self.option("ratio_exp").prop["path"], header=0, sep="\t")
        ana_path = os.path.join(self.work_dir, "protein_analysis.txt")
        p_ana = pd.read_table(ana_path, header=0, sep="\t")
        protein_selected = pd.DataFrame(index=protein.index)
        for col in ['Accession', 'Description',
                    # 'Significance', 'Coverage', '#Unique', 'Avg. Mass',
                    ]:
            try:
                protein_selected[col] = protein[col]
            except:
                try:
                    protein_selected[col] = protein['# ' + col]     # compatible to self.new in main workflow
                except:
                    self.logger.info('protein文件里缺少%s这一列，请检查'%col)
                    protein_selected[col] = '_'
        # protein_selected = protein.iloc[:, [2, 3, 4, 6, -3,-2]]
        protein_selected.set_index('Accession', inplace=True)
        ratio_exp.set_index('Accession', inplace=True)
        p_ana.set_index('Accession', inplace=True)
        indexes_to_drop = [x for x in protein_selected.index if x not in ratio_exp.index]
        protein_sliced = protein_selected.drop(list(indexes_to_drop))
        p_indexes_to_drop = [x for x in p_ana.index if x not in ratio_exp.index]
        p_ana = p_ana.drop(list(p_indexes_to_drop))
        protein_sliced.reset_index(inplace=True)
        p_ana.reset_index(inplace=True)
        ratio_exp.reset_index(inplace=True)
        protein_sliced.columns = [x.replace("# ", "") for x in list(protein_sliced.columns)]
        protein_sliced_final = pd.merge(protein_sliced, p_ana, on="Accession")
        protein_sliced_final.to_csv(self.work_dir + "/" + "protein_sliced.xls", sep='\t', index=False)
        search_db = pd.merge(protein_sliced_final, ratio_exp, on="Accession")
        search_db.to_csv(self.work_dir + "/" + "search_db.xls", sep='\t', index=False)

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))

        search_db = self.work_dir + "/" + "search_db.xls"
        os.link(search_db, os.path.join(self.output_dir, "search_db.xls"))
        self.logger.info("设置搜库结果目录")

    def run(self):
        super(SearchdbTool, self).run()
        self.get_protein_char()
        self.searchdb()
        self.set_output()
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
            "id": "dia_searchdb_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "dia_v3.searchdb",
            "instant": False,
            "options": dict(
                protein="/mnt/ilustre/users/sanger-dev/workspace/20210413/Dia_ev72_4tdv4i2ms31fp8tio3ck3k/remote_input/protein/protein_.xls",
                ratio_exp="/mnt/ilustre/users/sanger-dev/workspace/20210413/Dia_ev72_4tdv4i2ms31fp8tio3ck3k/remote_input/ratio_exp/exp.txt",
                protein_fasta="/mnt/ilustre/users/sanger-dev/workspace/20210413/Dia_ev72_4tdv4i2ms31fp8tio3ck3k/remote_input/protein_fasta/exp.fasta",
                report_dia="/mnt/ilustre/users/sanger-dev/workspace/20210413/Dia_ev72_4tdv4i2ms31fp8tio3ck3k/remote_input/report_dia/20210127_121022_Dl14ZRL_Report.xls",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)

