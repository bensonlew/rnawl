# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
# last_modify:2019.06.27

from biocluster.module import Module
import os
import shutil
import glob
import pandas as pd
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import unittest


class DiffCircleModule(Module):
    def __init__(self, work_id):
        super(DiffCircleModule, self).__init__(work_id)
        options = [
            {'name': 'diff_path', 'type': 'infile', 'format': 'itraq_and_tmt.common_dir'},
            {'name': 'kegg_enrich_path', 'type': 'infile', 'format': 'itraq_and_tmt.common_dir'},
            {'name': 'go_enrich_path', 'type': 'infile', 'format': 'itraq_and_tmt.common_dir'},
        ]
        self.add_option(options)
        self.tools = list()

    def transfer_enrich_file(self, enrich_file, out_prefix, type='go'):
        out = os.path.join(self.work_dir, '%s_circle_%s_enrich_table'%(type, out_prefix))
        with open(enrich_file) as er, open(out, 'w') as ow:
            _ = er.readline()
            if type == 'go':
                header = ["go_id", "go_type", "discription", "p_corrected", "p_uncorrected", "seq_list", "depth"]
                ow.write('\t'.join(header)+'\n')
                for line in er:
                    # self.logger.info(line)
                    line = line.strip().split('\t')
                    if not line or line[2] == 'p':
                        continue
                    data = dict((
                        ('go_id', line[0]),
                        ('go_type', line[1]),
                        ('discription', line[3]),
                        ('p_uncorrected', float(line[6])),
                        ('p_corrected', float(line[-1])),
                        ('depth', int(line[7])),
                        ('seq_list', line[-2]),
                    ))
                    write_line = '{' + '}\t{'.join(header)+'}\n'
                    # self.logger.info(write_line)
                    # self.logger.info(data['go_id'])
                    write_line = write_line.format(**data)
                    ow.write(write_line)
            else:
                header = ["id", "term", "pvalue", "corrected_pvalue", "seq_list", "kegg_type"]
                ow.write('\t'.join(header) + '\n')
                for line in er:
                    line = line.strip().split('\t')
                    if not line or line[2] == 'p':
                        continue
                    data = {
                        'term': line[1],
                        'id': line[3].split("path:")[1] if "path:" in line[3] else line[3],
                        'pvalue': float(line[6]),
                        'corrected_pvalue': float(line[7]) if not line[7] == "None" else "None",
                        'seq_list': line[8],
                        'kegg_type': "".join([x[0] for x in line[11].split(' ')])
                    }
                    write_line = '{' + '}\t{'.join(header) + '}\n'
                    write_line = write_line.format(**data)
                    ow.write(write_line)
        return out

    def run_circle(self):
        diff_files = glob.glob(os.path.join(self.option('diff_path').prop['path'], '*_vs_*'))
        for diff in diff_files:
            cmp = os.path.basename(diff).split('_diff.xls')[0]
            diff_df = pd.read_csv(diff, sep='\t', index_col=0)
            fc_de = os.path.join(self.work_dir, '%s_all_fc.list'%cmp)
            fc_up = os.path.join(self.work_dir, '%s_up_fc.list' % cmp)
            fc_down = os.path.join(self.work_dir, '%s_down_fc.list' % cmp)

            up_df = diff_df[((diff_df['regulate'] == 'up') & (diff_df['significant'] == 'yes'))]['log2fc']
            down_df = diff_df[((diff_df['regulate'] == 'down') & (diff_df['significant'] == 'yes'))]['log2fc']
            de_df = diff_df[(((diff_df['regulate'] == 'up') | (diff_df['regulate'] == 'down')) & (diff_df['significant'] == 'yes'))]['log2fc']
            up_df.to_csv(fc_up, sep='\t', header=True, index=True)
            down_df.to_csv(fc_down, sep='\t', header=True, index=True)
            de_df.to_csv(fc_de, sep='\t', header=True, index=True)

            go_enrich_up_ = glob.glob(os.path.join(self.option('go_enrich_path').prop['path'], '%s_up'%cmp, '*_enrich_*'))
            self.logger.info(go_enrich_up_)
            if go_enrich_up_ and up_df.shape[0] > 100:
                go_enrich = self.transfer_enrich_file(go_enrich_up_[0], cmp+'_up', type='go')
                options = dict(
                    enrich_type = 'GO',
                    enrich_table = go_enrich,
                    fc_table = fc_up
                )
                go_up_tool = self.add_tool("itraq_and_tmt.proteinset.enrich2circ")
                go_up_tool.set_options(options)
                self.tools.append(go_up_tool)
            go_enrich_down_ = glob.glob(
                os.path.join(self.option('go_enrich_path').prop['path'], '%s_down' % cmp, '*_enrich_*'))
            if go_enrich_down_ and down_df.shape[0] > 100:
                go_enrich = self.transfer_enrich_file(go_enrich_down_[0], cmp+'_down', type='go')
                options = dict(
                    enrich_type = 'GO',
                    enrich_table = go_enrich,
                    fc_table = fc_down
                )
                go_down_tool = self.add_tool("itraq_and_tmt.proteinset.enrich2circ")
                go_down_tool.set_options(options)
                self.tools.append(go_down_tool)
            go_enrich_de_ = glob.glob(
                os.path.join(self.option('go_enrich_path').prop['path'], '%s_all' % cmp, '*_enrich_*'))
            if go_enrich_de_ and de_df.shape[0] > 100:
                go_enrich = self.transfer_enrich_file(go_enrich_de_[0], cmp + '_all', type='go')
                options = dict(
                    enrich_type='GO',
                    enrich_table=go_enrich,
                    fc_table=fc_de
                )
                go_de_tool = self.add_tool("itraq_and_tmt.proteinset.enrich2circ")
                go_de_tool.set_options(options)
                self.tools.append(go_de_tool)

            kegg_enrich_up_ = glob.glob(
                os.path.join(self.option('kegg_enrich_path').prop['path'], '%s_up' % cmp, '*_enrichment.xls'))
            self.logger.info(kegg_enrich_up_)
            if kegg_enrich_up_ and up_df.shape[0] > 100:
                kegg_enrich = self.transfer_enrich_file(kegg_enrich_up_[0], cmp + '_up', type='kegg')
                options = dict(
                    enrich_type='KEGG',
                    enrich_table=kegg_enrich,
                    fc_table=fc_up
                )
                kegg_up_tool = self.add_tool("itraq_and_tmt.proteinset.enrich2circ")
                kegg_up_tool.set_options(options)
                self.tools.append(kegg_up_tool)
            kegg_enrich_down_ = glob.glob(
                os.path.join(self.option('kegg_enrich_path').prop['path'], '%s_down' % cmp, '*_enrichment.xls'))
            if kegg_enrich_down_ and down_df.shape[0] > 100:
                kegg_enrich = self.transfer_enrich_file(kegg_enrich_down_[0], cmp + '_down', type='kegg')
                options = dict(
                    enrich_type='KEGG',
                    enrich_table=kegg_enrich,
                    fc_table=fc_down
                )
                kegg_down_tool = self.add_tool("itraq_and_tmt.proteinset.enrich2circ")
                kegg_down_tool.set_options(options)
                self.tools.append(kegg_down_tool)
            kegg_enrich_de_ = glob.glob(
                os.path.join(self.option('kegg_enrich_path').prop['path'], '%s_all' % cmp, '*_enrichment.xls'))
            if kegg_enrich_de_ and de_df.shape[0] > 100:
                kegg_enrich = self.transfer_enrich_file(kegg_enrich_de_[0], cmp + '_all', type='kegg')
                options = dict(
                    enrich_type='KEGG',
                    enrich_table=kegg_enrich,
                    fc_table=fc_de
                )
                kegg_de_tool = self.add_tool("itraq_and_tmt.proteinset.enrich2circ")
                kegg_de_tool.set_options(options)
                self.tools.append(kegg_de_tool)

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
        for circle in self.tools:
            out_dir = os.path.join(self.output_dir, os.path.basename(circle.option('enrich_table').prop['path']).split('_enrich_table')[0])
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            for file in os.listdir(circle.output_dir):
                source = os.path.join(circle.output_dir,file)
                target = os.path.join(out_dir,file)
                if os.path.isfile(source):
                    os.link(source,target)
                else:
                    shutil.copytree(source, target)
        self.end()

    def run(self):
        super(DiffCircleModule, self).run()
        self.run_circle()
        if not self.tools:
            super(DiffCircleModule, self).end()
            return
        try:
            if len(self.tools) == 1:
                self.tools[0].on("end", self.set_output)
            else:
                self.on_rely(self.tools, self.set_output)
            for tool in self.tools:
                tool.run()
        except:
            super(DiffCircleModule, self).end()
            return

    def end(self):
        super(DiffCircleModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "diff_circle_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "itraq_and_tmt.diff_circle",
            "instant": False,
            "options": dict(
                diff_path="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/Diff/output",
                kegg_enrich_path="/mnt/ilustre/users/sanger-dev/workspace/20190626/Single_diff_kegg_enrich_5882_fyt/DiffKeggEnrich/output/",
                go_enrich_path='/mnt/ilustre/users/sanger-dev/workspace/20190625/Single_diff_go_enrich_1330_fyt/DiffGoEnrich/output/'
            )
        }

        data['id'] += '_fyt'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()