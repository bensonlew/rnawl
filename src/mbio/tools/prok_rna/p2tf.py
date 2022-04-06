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
import numpy as np
from mbio.packages.align.blast.xml2table import xml2table


class P2tfAgent(Agent):
    """
    TF prediction using P2TF database
    """
    def __init__(self, parent):
        super(P2tfAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", 'format': 'prok_rna.fasta'},
            {"name": "evalue", "type": "float", 'default': 1e-5},
        ]
        self.add_option(options)
        self.step.add_steps("tf_prediction")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.tf_prediction.start()
        self.step.update()

    def stepfinish(self):
        self.step.tf_prediction.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('query').is_set:
            raise OptionError("必须提供query序列文件")
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
        super(P2tfAgent, self).end()


class P2tfTool(Tool):
    def __init__(self, config):
        super(P2tfTool, self).__init__(config)
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin'))
        self.program = {
            'blastn': 'bioinfo/ref_rna_v2/miniconda2/bin/blastn',
            'python': 'miniconda2/bin/python',
        }
        self.file = {
            'p2tf_db': os.path.join(self.config.SOFTWARE_DIR, "database/prok_rna/P2TF_202106/p2tf"),
            # 'p2tf_db': "/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/prokrna_3.1/p2tf/p2tf",
            'p2tf_anno': os.path.join(self.config.SOFTWARE_DIR, "database/prok_rna/P2TF_202106/p2tf_anno.txt"),
            'xml': os.path.join(self.work_dir, os.path.splitext(os.path.basename(self.option('query').prop['path']))[0] + '_vs_P2TF.xml'),
            'domain': os.path.join(self.config.SOFTWARE_DIR, 'database/prok_rna/P2TF_202106/family_domain.txt'),
        }

    def run(self):
        """
        运行
        :return:
        """
        super(P2tfTool, self).run()
        self.run_blast()
        self.run_stats()
        self.set_output()
        self.end()

    def run_blast(self):
        self.logger.info('开始进行转录因子预测')
        cmd = "{} -task blastn-short -strand plus -num_threads 5 ".format(self.program['blastn'])
        cmd += "-db {} -query {} ".format(self.file['p2tf_db'], self.option('query').prop['path'])
        cmd += "-outfmt 5 -evalue {} ".format(self.option('evalue'))
        cmd += "-num_alignments 1 -out {}".format(self.file['xml'])
        self.logger.info("使用P2TF进行转录因子预测")
        command = self.add_command("p2tf_predict", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用P2TF进行转录因子预测完成!")
        else:
            self.set_error("使用P2TF进行转录因子预测出错！")
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
        annot_query = os.path.join(self.output_dir, 'p2tf_detail.xls')
        domain_dict = dict()
        blast_dict = dict()
        p2tf_anno = pd.read_table(self.file['p2tf_anno'], index_col=0, header=0)
        p2tf_anno = p2tf_anno[['class', 'type', 'P2TF description']]
        p2tf_anno.rename(columns={'P2TF description': 'des'}, inplace=True)
        anno_dict = p2tf_anno.to_dict('index')
        with open(self.file['domain'], 'r') as d:
            d.readline()
            for line in d:
                items = line.strip().split('\t')
                if items[0] not in domain_dict.keys():
                    domain_dict[items[0]] = {'arch': items[1], 'id': items[2]}
        with open(self.file['table_out'], 'r') as f:
            f.readline()
            for line in f:
                items = line.strip().split('\t')
                if items[0] in blast_dict.keys() and float(items[-5]) >= float(blast_dict[items[0]][-2]):
                    continue
                blast_dict[items[0]] = [items[2], items[-1], items[-5], items[-4]]
                tf = '_'.join(items[1].split('|')[0].split('_')[0:-2])
                anno = anno_dict.get(tf)
                if anno:
                    tf_type = str(anno.get('type'))
                    if domain_dict.get(tf_type):
                        blast_dict[items[0]].insert(0, str(domain_dict[tf_type]['arch']))
                        blast_dict[items[0]].insert(0, str(domain_dict[tf_type]['id']))
                    else:
                        blast_dict[items[0]].insert(0, '-')
                        blast_dict[items[0]].insert(0, '-')
                    blast_dict[items[0]].insert(0, str(anno['des']))
                    blast_dict[items[0]].insert(0, tf_type)
                    blast_dict[items[0]].insert(0, str(anno['class']))

        with open(annot_query, 'w') as a:
            a.write('gene\tclass\ttype\tp2tf_description\tdomain_id\tdomain_arch\tidentity\tcoverage\tevalue\tscore\n')
            for key, value in blast_dict.items():
                if len(value) > 5:
                    a.write(str(key) + '\t' + '\t'.join(value) + '\n')

    def set_output(self):
        blast = glob.glob(os.path.join(self.work_dir, '*_vs_P2TF.xls'))[0]
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
            "id": "p2tf_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "prok_rna_v3.p2tf",
            "instant": False,
            "options": dict(
                # query="/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210412/Prokrna_g0uk_3f5rs2hu3jkoaf37fqm5gd/RockhopperIndex/Rockhopper_Results/cds.fa",
                query='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/prok_v3.1/test.fa',
                evalue='1e-5',
                #config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/quantifier_test/Uniq.cfg.ini"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
