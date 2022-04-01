# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
# last_modify:2019.06.24

from biocluster.module import Module
import os
import shutil
import glob
import pandas as pd
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import unittest
from Bio.Blast import NCBIXML
from collections import Counter
import re


class DiffPpiModule(Module):
    def __init__(self, work_id):
        super(DiffPpiModule, self).__init__(work_id)
        options = [
            {'name': 'diff_path', 'type': 'infile', 'format': 'labelfree.common_dir'},
            {"name": "species", "type": "int", "default": -1},
            {"name": "combine_score", "type": "int", "default": 300},  # 设定蛋白质间的相互作用可能性值前300个互作组
            {"name": "seq", "type": "infile", "format": "labelfree.common"},
            {"name": "string_xml", "type": "string", "default": ""},
            {"name": "string_identity", "type": "int", "default": 98},
        ]
        self.add_option(options)
        self.modules = list()
        self.config = Config()

    def get_most_common_specie(self, xml):
        species = list()
        with open(xml, 'r') as blast_r:
            records = NCBIXML.parse(blast_r)
            for rec in records:
                for align in rec.alignments:
                    for hsp in align.hsps:
                        hit = align.hit_id
                        des = align.hit_def
                        if u'.' in des and '|' not in des:
                            hit = des
                        ident = float(hsp.identities)
                        hit_len = float(hsp.align_length)
                        if ident / hit_len * 100 >= self.option('string_identity'):
                            # specie = re.split('.', hit, maxsplit=1)[0]
                            specie = hit.split('.')[0]
                            # self.logger.info(hit)
                            species.append(specie)
        if not species:
            return 9606
        # 有的id在string数据库里面没有，只能再加一部判断
        # most_s = Counter(species).most_common()[0][0]
        with open(self.config.SOFTWARE_DIR + '/database/Annotation/all/String/string11.5/ppi_species.v11.5.txt') as ps:
            hit_species = [line.strip().split('\t')[1] for line in ps if line.strip()]
        for spe, count in Counter(species).most_common():
            if spe in hit_species:
                self.logger.info('选择第%s个物种%s 作为优势物种进行ppi分析'%(str(count), spe))
                return int(spe)
        return 9606

    def run_ppi(self):
        if self.option('species') == -1 or not self.option('species'):
            specie = self.get_most_common_specie(self.option('string_xml'))
        else:
            specie = self.option('species')
        diff_files = glob.glob(os.path.join(self.option('diff_path').prop['path'], '*_vs_*'))
        for diff in diff_files:
            cmp = os.path.basename(diff).split('_diff.xls')[0]
            acc_de = os.path.join(self.work_dir, '%s_all_protein.list'%cmp)
            acc_up = os.path.join(self.work_dir, '%s_up_protein.list' % cmp)
            acc_down = os.path.join(self.work_dir, '%s_down_protein.list' % cmp)
            diff_df = pd.read_csv(diff, sep='\t', index_col=0)
            uplist = diff_df[((diff_df['regulate'] == 'up') & (diff_df['significant'] == 'yes'))].index.tolist()
            downlist = diff_df[((diff_df['regulate'] == 'down') & (diff_df['significant'] == 'yes'))].index.tolist()
            with open(acc_de, 'w') as dew, open(acc_up, 'w') as uw, open(acc_down, 'w') as dw:
                dew.write('accession_id\n' + '\n'.join(uplist + downlist))
                uw.write('accession_id\n' + '\n'.join(uplist))
                dw.write('accession_id\n' + '\n'.join(downlist))
            options = {
                'diff_exp_gene': acc_de,
                'species': specie,
                'seq': self.option('seq').prop['path'],
                'combine_score': self.option('combine_score')
            }
            if len(uplist + downlist) > 10:
                de_module = self.add_module('labelfree.ppinetwork_analysis')
                de_module.set_options(options)
                self.modules.append(de_module)
            if len(uplist) > 10:
                up_module = self.add_module('labelfree.ppinetwork_analysis')
                options.update(dict(diff_exp_gene=acc_up))
                up_module.set_options(options)
                self.modules.append(up_module)
            if len(downlist) > 10:
                down_module = self.add_module('labelfree.ppinetwork_analysis')
                options.update(dict(diff_exp_gene=acc_down))
                down_module.set_options(options)
                self.modules.append(down_module)

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        return True

    def set_output(self):
        shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        for ppi in self.modules:
            out_dir = os.path.join(self.output_dir, os.path.basename(ppi.option('diff_exp_gene').prop['path']).split('_protein.list')[0])
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            for file in os.listdir(ppi.output_dir):
                source = os.path.join(ppi.output_dir,file)
                target = os.path.join(out_dir,file)
                if os.path.isfile(source):
                    os.link(source,target)
                else:
                    shutil.copytree(source, target)
        self.end()

    def run(self):
        super(DiffPpiModule, self).run()
        self.run_ppi()
        if not self.modules:
            super(DiffPpiModule, self).end()
        if len(self.modules) == 1:
            self.modules[0].on("end", self.set_output)
        else:
            self.on_rely(self.modules, self.set_output)
        for module in self.modules:
            module.run()

    def end(self):
        super(DiffPpiModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the module. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "diff_ppi_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "labelfree.diff_ppi",
            "instant": False,
            "options": dict(
                diff_path="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/Diff/output",
                seq="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/remote_input/protein_fasta/exp.fasta",
                string_xml='/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/ProteinAnnotation/output/blast_xml/string.xml'
            )
        }

        data['id'] += '_fyt'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
