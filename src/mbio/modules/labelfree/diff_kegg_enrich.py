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


class DiffKeggEnrichModule(Module):
    def __init__(self, work_id):
        super(DiffKeggEnrichModule, self).__init__(work_id)
        options = [
            {'name': 'diff_path', 'type': 'infile', 'format': 'labelfree.common_dir'},
            {'name': 'kegg_table', 'type': 'infile', 'format': 'labelfree.common'},
            {'name': 'pathway_table', 'type': 'infile', 'format': 'labelfree.common'},
            {'name': 'kegg_version', 'type': 'string', 'default': '2019'},
            {'name': 'png_dir', 'type': 'infile', 'format': 'itraq_and_tmt.common_dir'},
        ]
        self.add_option(options)
        self.tools = list()

    def get_input_files(self):
        add_info = os.path.join(self.work_dir, 'add_info')
        kegg_table = os.path.join(self.work_dir, 'protein_kegg_table.xls')
        pathway_table = os.path.join(self.work_dir, 'protein_kegg_level_table.xls')
        kegg_df = pd.read_csv(self.option('kegg_table').prop['path'], sep='\t', index_col=0)
        kegg_df = kegg_df.rename(index=str, columns={'KO_name (Gene name)':'KO_name(Protein name)'})
        if 'Gene_ID' in kegg_df.columns.tolist():
            kegg_df = kegg_df.rename(index=str, columns={'Gene_ID': 'KEGG_gene_id'})
        if 'KEGG_gene_id' not in kegg_df.columns.tolist() or 'kegg_gene_id' not in kegg_df.columns.tolist():
            kegg_df['KEGG_gene_id'] = pd.Series(['']*kegg_df.shape[0], index=kegg_df.index)
        kegg_df.to_csv(kegg_table, index=True, sep='\t', header=True)
        pathway_df = pd.read_csv(self.option('pathway_table').prop['path'], sep='\t', index_col=0)
        rename_dict = {
            'Pathway':'Pathway_id',
            'num_of_seqs':'number_of_seqs',
            'Pathway_definition':'pathway_definition',
            'First Category':'first_category',
            'Second Category':'second_category',
            'Hyperlink':'hyperlink',
            'seqs_kos/gene_list':'seq_list',
        }
        pathway_df = pathway_df.rename(index=str, columns=rename_dict)
        pathway_df['graph_id'] = pd.Series(['']*pathway_df.shape[0], index=pathway_df.index)
        pathway_df['graph_png_id'] = pd.Series(['']*pathway_df.shape[0], index=pathway_df.index)
        pathway_df['anno_type'] = pd.Series(['']*pathway_df.shape[0], index=pathway_df.index)
        pathway_df = pathway_df[['graph_id','number_of_seqs','pathway_definition','first_category','anno_type','hyperlink','seq_list','graph_png_id','second_category']]
        pathway_df.index.name = 'Pathway_id'
        pathway_df.to_csv(pathway_table, index=True, sep='\t', header=True)
        add_df = pathway_df['hyperlink']
        add_df.index.name = 'pathway'
        add_df.to_csv(add_info, index=True, sep='\t', header=True)
        return kegg_table, pathway_table, add_info

    def run_kegg_enrich(self, kegg_table, pathway_table, add_info):
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
                dew.write('\n'.join(uplist + downlist))
                uw.write('\n'.join(uplist))
                dw.write('\n'.join(downlist))

            acc_de_class = os.path.join(self.work_dir, '%s_all_protein_class.list'%cmp)
            acc_up_class = os.path.join(self.work_dir, '%s_up_protein_class.list' % cmp)
            acc_down_class = os.path.join(self.work_dir, '%s_down_protein_class.list' % cmp)
            with open(acc_de_class, 'w') as dew, open(acc_up_class, 'w') as uw, open(acc_down_class, 'w') as dw:
                dew.write('%s_all'%cmp + '\t' + ','.join(uplist + downlist))
                uw.write('%s_up'%cmp + '\t' + ','.join(uplist))
                dw.write('%s_down'%cmp + '\t' + ','.join(downlist))

            options = dict(
                task_id = 'itraq_and_tmt_default',
                diff_list = acc_de,
                kegg_table = kegg_table,
                kegg_version=self.option("kegg_version"),
                add_info = add_info,
                png_dir = self.option('png_dir').prop['path']
            )
            if len(uplist + downlist) > 10:
                de_tool = self.add_tool("itraq_and_tmt.proteinset.kegg_rich")
                options.update(dict(proteinset_kegg = acc_de_class))
                de_tool.set_options(options)
                self.tools.append(de_tool)
            if len(uplist) > 10:
                up_tool = self.add_tool("itraq_and_tmt.proteinset.kegg_rich")
                options.update(dict(diff_list=acc_up))
                options.update(dict(proteinset_kegg = acc_up_class))
                up_tool.set_options(options)
                self.tools.append(up_tool)
            if len(downlist) > 10:
                down_tool = self.add_tool("itraq_and_tmt.proteinset.kegg_rich")
                options.update(dict(diff_list=acc_down))
                options.update(dict(proteinset_kegg = acc_down_class))
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
        for kegg_enrich in self.tools:
            out_dir = os.path.join(self.output_dir, os.path.basename(kegg_enrich.option('diff_list').prop['path']).split('_protein.list')[0])
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            for file in os.listdir(kegg_enrich.output_dir):
                source = os.path.join(kegg_enrich.output_dir,file)
                target = os.path.join(out_dir,file)
                if os.path.isfile(source):
                    os.link(source,target)
                else:
                    shutil.copytree(source, target)

            rm_files = glob.glob(os.path.join(out_dir, 'pathways/*.png')) + \
                       glob.glob(os.path.join(out_dir, 'pathways/*.pdf'))
            for rm_file in rm_files:
                os.remove(rm_file)
        self.end()

    def run(self):
        super(DiffKeggEnrichModule, self).run()
        kegg_table, pathway_table, add_info = self.get_input_files()
        self.run_kegg_enrich(kegg_table, pathway_table, add_info)
        if not self.tools:
            super(DiffKeggEnrichModule, self).end()
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def end(self):
        super(DiffKeggEnrichModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "diff_kegg_enrich_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "labelfree.diff_kegg_enrich",
            "instant": False,
            "options": dict(
                diff_path="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/Diff/output",
                kegg_table="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/ProteinAnnotation/output/kegg/kegg_table.xls",
                pathway_table='/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/ProteinAnnotation/output/kegg/pathway_table.xls'
            )
        }

        data['id'] += '_fyt'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
