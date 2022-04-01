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


class DiffPfamStatModule(Module):
    def __init__(self, work_id):
        super(DiffPfamStatModule, self).__init__(work_id)
        options = [
            {'name': 'diff_path', 'type': 'infile', 'format': 'labelfree.common_dir'},
            {'name': 'pfam_file', 'type': 'infile', 'format': 'labelfree.common'},
        ]
        self.add_option(options)
        self.tools = list()

    def run_pfam_stat(self):
        # deal pfam_file
        pfam_df = pd.read_csv(self.option('pfam_file').prop['path'], sep='\t')
        rename_dict = {'Seq_id':'accession_id', 'Pfam_id':'pfam', 'Domain':'domain', 'DomainDescription':'description'}
        pfam_df = pfam_df.rename(columns=rename_dict)
        pfam_path = os.path.join(self.work_dir, 'pfam_file')
        pfam_df.to_csv(pfam_path, sep='\t', index=False, header=True)
        self.option('pfam_file', pfam_path)

        diff_files = glob.glob(os.path.join(self.option('diff_path').prop['path'], '*_vs_*'))
        for diff in diff_files:
            cmp = os.path.basename(diff).split('_diff.xls')[0]
            acc_de = os.path.join(self.work_dir, '%s_all_protein.list'%cmp)
            acc_up = os.path.join(self.work_dir, '%s_up_protein.list'%cmp)
            acc_down = os.path.join(self.work_dir, '%s_down_protein.list'%cmp)
            diff_df = pd.read_csv(diff, sep='\t', index_col=0)
            uplist = diff_df[((diff_df['regulate'] == 'up') & (diff_df['significant'] == 'yes'))].index.tolist()
            downlist = diff_df[((diff_df['regulate'] == 'down') & (diff_df['significant'] == 'yes'))].index.tolist()
            with open(acc_de, 'w') as dew, open(acc_up, 'w') as uw, open(acc_down, 'w') as dw:
                dew.write('\n'.join(uplist + downlist))
                uw.write('\n'.join(uplist))
                dw.write('\n'.join(downlist))
            options = dict(
                proteinsets = acc_de,
                pfam_info = self.option('pfam_file'),
            )
            if len(uplist + downlist) > 10:
                de_tool = self.add_tool("labelfree.export_pfam_stat")
                de_tool.set_options(options)
                self.tools.append(de_tool)
                up_down_tool = self.add_tool("labelfree.export_pfam_stat")
                options.update(dict(proteinsets=acc_up+','+acc_down))
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
        for pfam_stat in self.tools:
            output_file = os.path.basename(pfam_stat.option('proteinsets')).split('_protein.list')[0]
            if 'down' in output_file:
                output_file = output_file.replace('down', 'up_down')
            else:
                if 'up' in output_file:
                    output_file = output_file.replace('up', 'up_down')
            out_dir = os.path.join(self.output_dir, output_file)
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            for file in os.listdir(pfam_stat.output_dir):
                source = os.path.join(pfam_stat.output_dir,file)
                target = os.path.join(out_dir,file)
                if os.path.isfile(source):
                    os.link(source,target)
                else:
                    shutil.copytree(source, target)
        self.end()

    def run(self):
        super(DiffPfamStatModule, self).run()
        self.run_pfam_stat()
        if not self.tools:
            super(DiffPfamStatModule, self).end()
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def end(self):
        super(DiffPfamStatModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "diff_pfam_stat_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "labelfree.diff_pfam_stat",
            "instant": False,
            "options": dict(
                diff_path="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/Diff/output",
                pfam_file="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/ProteinAnnotation/output/blast_xml/pfam_domain",
            )
        }

        data['id'] += '_fyt'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()