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


class DiffClusterModule(Module):
    def __init__(self, work_id):
        super(DiffClusterModule, self).__init__(work_id)
        options = [
            {'name': 'diff_path', 'type': 'infile', 'format': 'itraq_and_tmt.common_dir'},
            {'name': 'exp', 'type': 'infile', 'format': 'itraq_and_tmt.ratio_exp'},
            {'name': 'group', 'type': 'infile', 'format': 'itraq_and_tmt.group_table'},
            {'name': 'sct', 'type': 'string', 'default': 'hierarchy'},
            {'name': 'gct', 'type': 'string', 'default': 'hierarchy'},
            {'name': 'scm', 'type': 'string', 'default': 'complete'},
            {'name': 'gcm', 'type': 'string', 'default': 'average'},
            {'name': 'scd', 'type': 'string', 'default': 'correlation'},
            {'name': 'gcd', 'type': 'string', 'default': 'euclidean'},
            {'name': 'output', 'type': 'string', 'default': None},
        ]
        self.add_option(options)
        self.tools = list()


    def prepare_exp(self):
        diff_list = list()
        exp_df = pd.read_csv(self.option('exp').prop["path"], sep='\t', index_col=0)
        diff_files = glob.glob(os.path.join(self.option('diff_path').prop['path'],'*_vs_*'))
        for diff_file in diff_files:
            diff_pd = pd.read_csv(diff_file, sep='\t', index_col=0)
            diff_proteins = diff_pd[(((diff_pd['regulate'] == 'up') | (diff_pd['regulate'] == 'down')) & (diff_pd['significant'] == 'yes'))].index.tolist()
            diff_list += diff_proteins

        diff_list = list(set(diff_list))
        if len(diff_list) < 5:
            self.logger.info('差异蛋白只有%s个，cluster退出运行', len(diff_list))
            return None
        else:
            new_df = exp_df.loc[diff_list,]
            exp_path = os.path.join(self.work_dir,'all_diff.exp.txt')
            new_df.to_csv(exp_path, sep='\t', header=True, index=True)
            return exp_path

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        return True

    def run_clusters(self):
        exp_path = self.option('exp').prop['path']
        exp_df = pd.read_csv(exp_path, sep='\t', index_col=0)
        # for i in [5, 10, 20, 40, 80]:
        for i in [5, 10]:
            if exp_df.shape[0] < i:
                i = exp_df.shape[0]
            cluster_yes = self.add_tool("itraq_and_tmt.exp_cluster")
            cluster_no = self.add_tool("itraq_and_tmt.exp_cluster")
            options = dict(
                exp=self.option('exp'),
                group=self.option('group'),
                n_clusters=i,
                sct=self.option('sct'),
                gct=self.option('gct'),
                scm=self.option('scm'),
                gcm=self.option('gcm'),
                scd=self.option('scd'),
                gcd=self.option('gcd'),
                use_group='yes',
            )
            cluster_yes.set_options(options)
            self.tools.append(cluster_yes)
            options.update(dict(use_group='no'))
            cluster_no.set_options(options)
            self.tools.append(cluster_no)
            if exp_df.shape[0] <= 5:
                break

    def set_output(self):
        shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        for cluster in self.tools:
            if cluster.option('use_group') == 'no':
                perfix = 'sample'
            else:
                perfix = 'group'
            out_dir = os.path.join(self.output_dir,'cluster_%s_%s'%(perfix,str(cluster.option('n_clusters'))))
            os.mkdir(out_dir)
            for file in os.listdir(cluster.output_dir):
                source = os.path.join(cluster.output_dir,file)
                target = os.path.join(out_dir,file)
                os.link(source,target)
        self.end()

    def run(self):
        super(DiffClusterModule, self).run()
        exp_path = self.prepare_exp()
        if exp_path is not None:
            self.option('exp', exp_path)
            self.run_clusters()
            if len(self.tools) == 1:
                self.tools[0].on("end", self.set_output)
            else:
                self.on_rely(self.tools, self.set_output)
            for tool in self.tools:
                tool.run()
        else:
            super(DiffClusterModule, self).end()
        
    def end(self):
        super(DiffClusterModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "diff_cluster_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "itraq_and_tmt.diff_cluster",
            "instant": False,
            "options": dict(
                diff_path="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/MJ20200807109/",
                exp = "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/MJ20200807109/3167.exp.txt",
                group="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/MJ20200807109//group.txt",
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