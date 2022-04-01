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
from mbio.packages.itraq_and_tmt.kegg_enrichment import multtest_correct


class DiffGoEnrichModule(Module):
    def __init__(self, work_id):
        super(DiffGoEnrichModule, self).__init__(work_id)
        options = [
            {'name': 'diff_path', 'type': 'infile', 'format': 'itraq_and_tmt.common_dir'},
            {'name': 'go_list', 'type': 'infile', 'format': 'itraq_and_tmt.common'},
            {"name": "go_version", "type": "string", "default": "2019"}, #pir database version
        ]
        self.add_option(options)
        self.tools = list()

    def run_go_enrich(self):
        have_go_accs = {l.strip().split('\t')[0] for l in open(self.option('go_list').prop['path'])}
        diff_files = glob.glob(os.path.join(self.option('diff_path').prop['path'], '*_vs_*'))
        # go_df = pd.read_csv(self.option('go_list').prop['path'], sep='\t', header=None, index_col=0, names=['go'])
        for diff in diff_files:
            cmp = os.path.basename(diff).split('_diff.xls')[0]
            acc_de = os.path.join(self.work_dir, '%s_all_protein.list'%cmp)
            acc_up = os.path.join(self.work_dir, '%s_up_protein.list'%cmp)
            acc_down = os.path.join(self.work_dir, '%s_down_protein.list'%cmp)
            diff_df = pd.read_csv(diff, sep='\t', index_col=0)
            uplist = diff_df[((diff_df['regulate'] == 'up') & (diff_df['significant'] == 'yes'))].index.tolist()
            downlist = diff_df[((diff_df['regulate'] == 'down') & (diff_df['significant'] == 'yes'))].index.tolist()
            with open(acc_de, 'w') as dew, open(acc_up, 'w') as uw, open(acc_down, 'w') as dw:
                dew.write('\n'.join(uplist + downlist) + "\n")
                uw.write('\n'.join(uplist) + "\n")
                dw.write('\n'.join(downlist) + "\n")
            options = dict(
                diff_list = acc_de,
                go_version = self.option('go_version'),
                go_list = self.option('go_list').prop['path']
            )
            if len(uplist + downlist) > 10 and set(uplist + downlist) & have_go_accs:
                de_tool = self.add_tool("itraq_and_tmt.proteinset.go_enrich")
                de_tool.set_options(options)
                self.tools.append(de_tool)
            if len(uplist) > 10 and set(uplist) & have_go_accs:
                up_tool = self.add_tool("itraq_and_tmt.proteinset.go_enrich")
                options.update(dict(diff_list=acc_up))
                up_tool.set_options(options)
                self.tools.append(up_tool)
            if len(downlist) > 10 and set(downlist)  & have_go_accs:
                down_tool = self.add_tool("itraq_and_tmt.proteinset.go_enrich")
                options.update(dict(diff_list=acc_down))
                down_tool.set_options(options)
                self.tools.append(down_tool)

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
        for go_enrich in self.tools:
            out_dir = os.path.join(self.output_dir, '_'.join(os.path.basename(go_enrich.option('diff_list').prop['path']).split('_')[:-1]))
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            for file in os.listdir(go_enrich.output_dir):
                source = os.path.join(go_enrich.output_dir,file)
                target = os.path.join(out_dir,file)
                if not file.endswith('.xls'):
                    os.link(source,target)
                else:
                    pf = pd.read_table(source, header=0, sep="\t")
                    pf['p_corrected'] = multtest_correct(pf['p_uncorrected'].tolist(),3)
                    pf = pf.drop(pf.columns.tolist()[9], axis=1)
                    pf.to_csv(target, sep='\t', index=False)
        self.end()

    def run(self):
        super(DiffGoEnrichModule, self).run()
        self.run_go_enrich()
        if not self.tools:
            super(DiffGoEnrichModule, self).end()
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def end(self):
        super(DiffGoEnrichModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "diff_go_enrich_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "itraq_and_tmt.diff_go_enrich",
            "instant": False,
            "options": dict(
                diff_path="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/Diff/output",
                go_list="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/ProteinAnnotation/output/go/query_gos.list",
            )
        }
        # data['options']['method'] = 'rsem'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        # #
        # data['id'] += '1'
        # data['options']['method'] = 'salmon'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        #
        data['id'] += '_fyt'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
