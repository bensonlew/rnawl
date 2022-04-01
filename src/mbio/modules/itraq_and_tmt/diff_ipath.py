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


class DiffIpathModule(Module):
    def __init__(self, work_id):
        super(DiffIpathModule, self).__init__(work_id)
        options = [
            {'name': 'diff_path', 'type': 'infile', 'format': 'itraq_and_tmt.common_dir'},
            {'name': 'kegg_table', 'type': 'infile', 'format': 'itraq_and_tmt.common'},
        ]
        self.add_option(options)
        self.tools = list()

    def get_input_files(self):
        kegg_table = os.path.join(self.work_dir, 'protein_kegg_table.xls')
        kegg_df = pd.read_csv(self.option('kegg_table').prop['path'], sep='\t', index_col=0)
        kegg_df = kegg_df.rename(index=str, columns={'KO_name (Gene name)':'KO_name(Protein name)'})
        if 'KEGG_gene_id' not in kegg_df.columns.tolist() or 'kegg_gene_id' not in kegg_df.columns.tolist():
            kegg_df['KEGG_gene_id'] = pd.Series(['']*kegg_df.shape[0], index=kegg_df.index)
        kegg_df.to_csv(kegg_table, index=True, sep='\t', header=True)

        return kegg_table

    def run_ipath(self, kegg_table):
        diff_files = glob.glob(os.path.join(self.option('diff_path').prop['path'], '*_vs_*'))
        for diff in diff_files:
            cmp = os.path.basename(diff).split('_diff.xls')[0]
            acc_de = os.path.join(self.work_dir, '%s_all_protein.list'%cmp)
            acc_up_down = os.path.join(self.work_dir, '%s_up_down_protein.list'%cmp)
            diff_df = pd.read_csv(diff, sep='\t', index_col=0)
            uplist = diff_df[((diff_df['regulate'] == 'up') & (diff_df['significant'] == 'yes'))].index.tolist()
            downlist = diff_df[((diff_df['regulate'] == 'down') & (diff_df['significant'] == 'yes'))].index.tolist()
            with open(acc_de, 'w') as dew, open(acc_up_down, 'w') as udw:
                dew.write('%s_all'%cmp + '\t' + ','.join(uplist + downlist))
                udw.write('%s_up'%cmp + '\t' + ','.join(uplist) + '\n' + '%s_down'%cmp + '\t' + ','.join(downlist))
            options = dict(
                proteinset_kegg = acc_de,
                kegg_table = kegg_table,
            )
            if len(uplist + downlist) > 50:
                de_tool = self.add_tool("itraq_and_tmt.proteinset.ipath")
                de_tool.set_options(options)
                self.tools.append(de_tool)
                up_down_tool = self.add_tool("itraq_and_tmt.proteinset.ipath")
                options.update(dict(proteinset_kegg=acc_up_down))
                up_down_tool.set_options(options)
                self.tools.append(up_down_tool)

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
        for ipath in self.tools:
            out_dir = os.path.join(self.output_dir, os.path.basename(ipath.option('proteinset_kegg')).split('_protein.list')[0])
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            for file in os.listdir(ipath.output_dir):
                source = os.path.join(ipath.output_dir,file)
                target = os.path.join(out_dir,file)
                if os.path.isfile(source):
                    os.link(source,target)
                else:
                    shutil.copytree(source, target)
        self.end()

    def run(self):
        super(DiffIpathModule, self).run()
        kegg_table = self.get_input_files()
        self.run_ipath(kegg_table)
        if not self.tools:
            super(DiffIpathModule, self).end()
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def end(self):
        super(DiffIpathModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "diff_ipath_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "itraq_and_tmt.diff_ipath",
            "instant": False,
            "options": dict(
                diff_path="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/Diff/output",
                kegg_table="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/ProteinAnnotation/output/kegg/kegg_table.xls",
            )
        }

        data['id'] += '_fyt'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()