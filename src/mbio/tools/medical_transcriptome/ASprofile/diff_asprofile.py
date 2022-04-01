# coding=utf-8
import os
import unittest
import pandas as pd
import itertools
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
import json
__author__ = 'zoujiaxun'


class DiffAsprofileAgent(Agent):
    """

    """
    def __init__(self, parent):
        super(DiffAsprofileAgent, self).__init__(parent)
        options = [
            dict(name='AS_result_merge', type='infile', format='ref_rna_v2.common'),
            dict(name='sample', type='string'),
            dict(name='group_dict', type='string'),
            dict(name='filter', type='int'),
            dict(name='diff_as', type='outfile', format= 'ref_rna_v2.common')


        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        pass

    def set_resource(self):
        self._cpu = 1

        self._memory = '100G'

    def end(self):
        super(DiffAsprofileAgent, self).end()


class DiffAsprofileTool(Tool):
    def __init__(self, config):
        super(DiffAsprofileTool, self).__init__(config)

        self.file = {
            'diff_as': os.path.join(self.output_dir, 'AS_diff.txt')
        }

    def run(self):
        super(DiffAsprofileTool, self).run()
        self.run_diff()
        self.set_output()
        self.end()

    def run_diff(self):
        asprofile = pd.read_table(self.option('AS_result_merge').path, sep='\t')
        if self.option('sample'):
            sample = self.option('sample')
            sample_list = sample.split(',')
            asprofile_use = asprofile[['event_type', 'gene_id', 'chrom', 'event_start',
                                                                                 'event_end', 'event_pattern', 'strand']]
            event_dict = dict()
            index_num = asprofile.index
            with open(self.file['diff_as'], 'w') as diff_as:
                diff_as.write('event_id\tevent_type\tgene_id\tchrom\tevent_start\tevent_end\tevent_pattern\tstrand\t{}\n'.format(
                    '\t'.join(sample_list)
                ))

                for i in index_num:
                    a = list(asprofile_use.loc[i])
                    b = [str(c) for c in a]
                    d = ';'.join(b)
                    if d not in event_dict:
                        event_dict[d] = [asprofile.loc[i]['sample']]
                    else:
                        event_dict[d].append(asprofile.loc[i]['sample'])
                num = 1000000
                for i in event_dict:
                    list1 = i.split(';')
                    list2 = ['yes' if sample in event_dict[i] else 'no' for sample in sample_list]
                    num += 1
                    envent_id = 'event_{}'.format(num)
                    diff_as.write(envent_id+'\t'+ '\t'.join(list1) +'\t' + '\t'.join(list2) + '\n')


            #
            #
            #
            # dfs = []
            # for i in sample_list:
            #     sample_asprofile = asprofile[asprofile['sample']==i].copy()
            #     sample_asprofile['sample'] = 'yes'
            #     sample_asprofile.drop(labels=['event_id', 'group'],axis=1,inplace=True)
            #     name = '{}'.format(i)
            #     sample_asprofile.rename(columns={'sample':name},inplace=True)
            #     dfs.append(sample_asprofile)
            # merges = reduce(lambda left, right: pd.merge(left, right, how='outer', on=['event_type', 'gene_id', 'chrom', 'event_start',
            #                                                                      'event_end', 'event_pattern', 'strand']), dfs)
            # final = merges.fillna('no')
            # # final = merges
            # final['event_id']=['event_{}'.format(i) for i in range(1000000+1, 1000000+len(final)+1)]
            # final.to_csv(self.file['diff_as'], index=False, sep='\t')

        else:
            group_dict = json.loads(self.option('group_dict'))
            sample_list = []
            for group in group_dict:
                sample_list.extend(group_dict[group])
            asprofile_sample=asprofile[asprofile['sample'].isin(sample_list)]
            filter_list = []
            group_list = list(group_dict)
            for j in group_list:
                group_asprofile = asprofile_sample[asprofile_sample['group']==j]
                group_asprofile.drop(labels=['event_id'], axis=1, inplace=True)
                a = group_asprofile[group_asprofile.duplicated(keep=False,subset=['event_type', 'gene_id', 'chrom', 'event_start',
                                                                              'event_end','event_pattern', 'strand'])]
                df = a.groupby(
                    ['event_type', 'gene_id', 'chrom', 'event_start', 'event_end', 'event_pattern', 'strand']).apply(
                    lambda x: tuple(x.index)).tolist()
                row_list = []
                for m in df:
                    if len(m)>(int(self.option('filter'))/100)*len(group_dict[j]):
                        row_list.extend(m)
                print row_list
                group_asprofile_filter = group_asprofile.loc[row_list,].sort_index()
                filter_list.append(group_asprofile_filter)
            merge_filter = pd.concat(filter_list, axis=0, ignore_index=False)
            print filter_list
            print merge_filter

            asprofile_use = merge_filter[['event_type', 'gene_id', 'chrom', 'event_start',
                                                                                 'event_end', 'event_pattern', 'strand']]
            event_dict = dict()
            index_num = merge_filter.index
            with open(self.file['diff_as'], 'w') as diff_as:
                diff_as.write('event_id\tevent_type\tgene_id\tchrom\tevent_start\tevent_end\tevent_pattern\tstrand\t{}\n'.format(
                    '\t'.join(group_list)
                ))

                for i in index_num:
                    a = list(asprofile_use.loc[i])
                    b = [str(c) for c in a]
                    d = ';'.join(b)
                    if d not in event_dict:
                        event_dict[d] = [merge_filter.loc[i]['group']]
                    else:
                        event_dict[d].append(merge_filter.loc[i]['sample'])
                num = 1000000
                for i in event_dict:
                    list1 = i.split(';')
                    list2 = ['yes' if group in event_dict[i] else 'no' for group in group_list]
                    num += 1
                    envent_id = 'event_{}'.format(num)
                    diff_as.write(envent_id+'\t'+ '\t'.join(list1) +'\t' + '\t'.join(list2) + '\n')

            #
            #
            #
            #
            # dfs = []
            # for i in group_list:
            #     group_asprofile = merge_filter[merge_filter['group']==i]
            #     group_asprofile['group'] = 'yes'
            #     group_asprofile.drop(labels=['sample'],axis=1,inplace=True)
            #     name = '{}'.format(i)
            #     group_asprofile.rename(columns={'group':name},inplace=True)
            #     group_asprofile.drop_duplicates(inplace=True)
            #     dfs.append(group_asprofile)
            # merges = reduce(lambda left, right: pd.merge(left, right, how='outer', on=['event_type', 'gene_id', 'chrom', 'event_start',
            #                                                                      'event_end', 'event_pattern', 'strand']), dfs)
            # final = merges.fillna('no')
            # final['event_id']=['event_{}'.format(i) for i in range(1000000+1, 1000000+len(final)+1)]
            # final.to_csv(self.file['diff_as'], index=False, sep='\t')

    def set_output(self):
        # all ready write results to output
        self.option('diff_as').set_path(self.file['diff_as'])


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "diff_as" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.ASprofile.diff_asprofile",
            "instant": True,
            "options": {
                "AS_result_merge": '/mnt/ilustre/users/sanger-dev/workspace/20200603/Asprofile_tsg_37468_5662_6332/Asprofile/output/AS_result_merge.txt',
                'sample': 'DM_TB1,DM_TB2,NDM_TB_C2,NDM_TB_C3',
                # 'filter':50
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

