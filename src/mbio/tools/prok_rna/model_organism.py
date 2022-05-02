# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import ConfigParser
import unittest
import glob
import pandas as pd
from mbio.packages.align.blast.xml2table import xml2table
import numpy as np


class ModelOrganismAgent(Agent):
    """
    Model organism annotation using iModulonDB/RegulonDB database
    """
    def __init__(self, parent):
        super(ModelOrganismAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", 'format': 'prok_rna.fasta'},
            {"name": "database", "type": "string", 'default': 'regulon'},
            {"name": "species", "type": "string", 'default': 'e_coli'},
            {"name": "evalue", "type": "float", 'default': 1e-5},
        ]
        self.add_option(options)
        self.step.add_steps("model_organism")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.model_organism.start()
        self.step.update()

    def stepfinish(self):
        self.step.model_organism.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('query').is_set:
            raise OptionError("必须提供query序列文件")
        if self.option('database') == 'regulon' and self.option('species') != 'e_coli':
            raise OptionError("RegulonDB暂不支持该模式物种")
        if self.option('database') == 'imodulon' and self.option('species') not in ['e_coli', 'b_subtilis', 's_aureus']:
            raise OptionError("iModulonDB暂不支持该模式物种")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(ModelOrganismAgent, self).end()


class ModelOrganismTool(Tool):
    def __init__(self, config):
        super(ModelOrganismTool, self).__init__(config)
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin'))
        self.program = {
            'blastn': 'miniconda2/bin/blastn',
            'python': 'miniconda2/bin/python',
        }
        self.file = {
            'db': os.path.join(self.config.SOFTWARE_DIR, "database/prok_rna"),
            'xml': os.path.join(self.work_dir, os.path.splitext(os.path.basename(self.option('query').prop['path']))[
                0] + '_vs_{}.xml'.format(self.option('database'))),
        }

    def run(self):
        """
        运行
        :return:
        """
        super(ModelOrganismTool, self).run()
        self.run_blast()
        self.run_stats()
        self.set_output()
        self.end()

    def run_blast(self):
        if self.option('database') == 'imodulon':
            refdb = os.path.join(self.file['db'], 'iModulonDB/{species}/{species}'.format(species=self.option('species')))
        else:
            refdb = os.path.join(self.file['db'], 'RegulonDB-10.6/regulondb')
        self.logger.info('开始进行模式物种注释')
        cmd = "{} -task blastn-short -strand plus -num_threads 5 ".format(self.program['blastn'])
        cmd += "-db {} -query {} ".format(refdb, self.option('query').prop['path'])
        cmd += "-outfmt 5 -evalue {} ".format(self.option('evalue'))
        cmd += "-num_alignments 1 -out {}".format(self.file['xml'])
        self.logger.info("使用{}进行模式物种注释".format(self.option('database')))
        command = self.add_command("model_organism_annot", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用{}进行模式物种注释完成！".format(self.option('database')))
        else:
            self.set_error("使用{}进行模式物种注释出错！".format(self.option('database')))
        self.xml2table()

    def xml2table(self):
        "xml格式转换table"
        xml_fp = self.file['xml']
        table = os.path.basename(xml_fp)[:-3] + "xls"
        self.file['table_out'] = os.path.join(self.work_dir, table)
        self.xml2table_coverage(xml_fp, self.file['table_out'], anno_head=True)

    def xml2table_coverage(self, xml_fp, table_out, head=None, anno_head=False, query=True, hit_ref=True):
        """
        :修改表头，计算覆盖度  zouguanqing 201903
        """
        if not head:
            head = ["Query-Name", "Hit-Name", "Identity-%", "HSP-Len", "Mismatch", "Gapopen_num",
                    "Q-Begin", "Q-End", "Hsp-Begin", "Hsp-End", 'E-Value', "Score", "Q-Len", "Hit-Len"]
            if hit_ref == True:
                head[1] = "Hit-Description"
        tmp_table = table_out + '_tmp'
        xml2table(xml_fp, tmp_table, header=head, anno_head=True)

        f = pd.read_table(tmp_table, sep='\t')
        if query:
            f['coverage'] = np.abs((f['Q-End'] - f['Q-Begin']) / f['Q-Len'])
        else:
            f['coverage'] = f['HSP-Len'] / f['Hit-Len']
        if anno_head:
            f.to_csv(table_out, sep='\t', index=False, header=True)
        else:
            f.to_csv(table_out, sep='\t', index=False, header=False)

    def run_stats(self):
        annot_query = os.path.join(self.output_dir, '{}_{}_detail.xls'.format(self.option('database'), self.option('species')))
        if self.option('database') == 'imodulon':
            anno_file = os.path.join(self.file['db'], 'iModulonDB/{}/merged_anno.txt'.format(self.option('species')))
            header = 'gene_id\thit_id\thit_name\toperon\tcog\tgene_product\ttf\tim_name\tim_reg\tim_func\tim_cat\tidentity\tcoverage\tevalue\tscore\n'
        else:
            anno_file = os.path.join(self.file['db'], 'RegulonDB-10.6/merged_anno.txt')
            header = 'gene_id\thit_id\thit_name\tproduct\ttf\ttf_class\ttu\tpromoter\tterminator_class\tregulon\tregulon_func\tidentity\tcoverage\tevalue\tscore\n'
        anno_dict = dict()
        blast_dict = dict()
        with open(anno_file, 'r') as d:
            d.readline()
            for line in d:
                items = line.strip().split('\t', 1)
                if items[0] not in anno_dict.keys():
                    anno_dict[items[0]] = items[1]

        with open(self.file['table_out'], 'r') as f:
            f.readline()
            for line in f:
                items = line.strip().split('\t')
                if items[0] in blast_dict.keys() and float(items[-5]) >= float(blast_dict[items[0]][-2]):
                    continue
                if self.option('database') == 'imodulon':
                    ref_id = items[1].split(' ')[0].split('gene-')[1]
                    blast_dict[items[0]] = [ref_id, items[2], items[-1], items[-5], items[-4]]
                else:
                    ref_id, ref_name = items[1].split('|')
                    blast_dict[items[0]] = [ref_id, ref_name, items[2], items[-1], items[-5], items[-4]]

        with open(annot_query, 'w') as a:
            a.write(header)
            for key, value in blast_dict.items():
                anno = anno_dict.get(value[0])
                if anno and len(value) == 5:
                    a.write(str(key) + '\t' + value[0] + '\t' + anno + '\t' + '\t'.join(value[1:]) + '\n')
                elif anno and len(value) == 6:
                    a.write(str(key) + '\t' + '\t'.join(value[0:2]) + '\t' + anno + '\t' + '\t'.join(value[2:]) + '\n')

    def set_output(self):
        blast = glob.glob(os.path.join(self.work_dir, '*_vs_*.xls'))[0]
        link = os.path.join(self.output_dir, os.path.basename(blast))
        if os.path.exists(link):
            os.remove(link)
        os.link(blast, link)


class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "model_organism_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "prok_rna_v3.model_organism",
            "instant": False,
            "options": dict(
                # query="/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210412/Prokrna_g0uk_3f5rs2hu3jkoaf37fqm5gd/RockhopperIndex/Rockhopper_Results/cds.fa",
                query='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/prok_v3.1/test.fa',
                database='regulon',
                species='e_coli',
                evalue='1e-5',
                #config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/quantifier_test/Uniq.cfg.ini"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()