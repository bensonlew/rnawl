# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'liubinxu'


class FastqMerge2Agent(Agent):
    """
    Merge2 fastq file
    """
    def __init__(self, parent):
        super(FastqMerge2Agent, self).__init__(parent)
        options = [
            {'type': 'string', 'name': 'sample_fq_list'},
            {'type': 'infile', 'name': 'filter_sample', 'format': 'denovo_rna_v2.common'},
            {'default': '', 'type': 'string', 'name': 'filter_sample_list'},
            {'format': 'denovo_rna_v2.common', 'type': 'outfile', 'name': 'merge_left'},
            {'format': 'denovo_rna_v2.common', 'type': 'outfile', 'name': 'merge_right'},
            {"name": "group", "type": "infile", "format": "denovo_rna_v2.group_table"},
            {"name": "assemble_method", "type": "string", "default": "total"}, # 选择拼接软件sample, group, total
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('sample_fq_list'):
            raise OptionError("样本列表文件需提供依次为sample_id\tleft.fq\tright.fq", code = "32007601")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('5')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(FastqMerge2Agent, self).end()


class FastqMerge2Tool(Tool):
    """
    Merge2 fastq file
    """
    def __init__(self, config):
        super(FastqMerge2Tool, self).__init__(config)
        self.filter_samples = []
        self.right_fq = []
        self.left_fq = []
        self.merge_left = self.work_dir +  '/merge.left.fq'
        self.merge_right = self.work_dir +  '/merge.right.fq'

    def get_filter_samples(self):
        """
        获取过滤的样本信息
        """
        if self.option('filter_sample').prop.has_key('path'):
            with open(self.option('filter_sample').prop['path'], 'r') as f:
                lines = f.readlines()
                for line in lines:
                    self.filter_samples.append(line.strip())
        elif self.option('filter_sample_list'):
            self.filter_samples = self.option('filter_sample_list').split(',')
        else:
            pass

    def run_fastq_merge(self):
        with open(self.option("sample_fq_list")) as fq_list:
            lines = fq_list.readlines()
            for line in lines:
                line = line.strip().split('\t')
                if line[0] in self.filter_samples:
                    pass
                else:
                    if len(line) > 2:
                        self.right_fq.append(line[2])
                    else:
                        pass
                    self.left_fq.append(line[1])
            with open(self.merge_left, 'w') as left:
                for file in self.left_fq:
                    for fq in open(file, 'r'):
                        left.write(fq)
            if len(self.right_fq) > 0:
                with open(self.merge_right, 'w') as right:
                    for file in self.right_fq:
                        for fq in open(file, 'r'):
                            right.write(fq)

    def run_fastq_merge_sample(self):
        left_fq_dict, right_fq_dict = self.get_sample_fq()
        for sample in left_fq_dict:
            merge_left = self.work_dir +  '/{}.merge.left.fq'.format(sample)
            merge_right = self.work_dir +  '/{}.merge.right.fq'.format(sample)
            with open(merge_left, 'w') as left:
                if sample in left_fq_dict:
                    for fq in left_fq_dict[sample]:
                        for line in open(fq, 'r'):
                            left.write(line)
            with open(merge_right, 'w') as right:
                if sample in right_fq_dict:
                    for fq in right_fq_dict[sample]:
                        for line in open(fq, 'r'):
                            right.write(line)
            
        # 数据不用合并
        return True

    def get_sample_fq(self):
        left_fq_dict = {}
        right_fq_dict = {}
        with open(self.option("sample_fq_list")) as fq_list:
            lines = fq_list.readlines()
            for line in lines:
                line = line.strip().split('\t')
                if line[0] in self.filter_samples:
                    pass
                else:
                    if len(line) > 2:
                        if line[0] in right_fq_dict:
                            right_fq_dict[line[0]].append(line[2])
                        else:
                            right_fq_dict[line[0]] = [line[2]]
                    else:
                        pass
                    if line[0] in left_fq_dict:
                        left_fq_dict[line[0]].append(line[1])
                    else:
                        left_fq_dict[line[0]] = [line[1]]
        return left_fq_dict, right_fq_dict

    def run_fastq_merge_group(self):
        group_spname =  self.option("group").get_group_spname()
        left_fq_dict, right_fq_dict = self.get_sample_fq()
        for group,samples in group_spname.items():
            merge_left = self.work_dir +  '/{}.merge.left.fq'.format(group)
            merge_right = self.work_dir +  '/{}.merge.right.fq'.format(group)
            with open(merge_left, 'w') as left:
                for sample in samples:
                    if sample in left_fq_dict:
                        for fq in left_fq_dict[sample]:
                            for line in open(fq, 'r'):
                                left.write(line)
            with open(merge_right, 'w') as right:
                for sample in samples:
                    if sample in right_fq_dict:
                        for fq in right_fq_dict[sample]:
                            for line in open(fq, 'r'):
                                right.write(line)

    def set_output(self):
        left_files = glob.glob('./*.left.fq')
        right_files = glob.glob('./*.right.fq')
        all_files = left_files + right_files
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
        '''
        if os.path.exists(self.output_dir + '/merge.left.fq'):
            os.remove(self.output_dir + '/merge.left.fq')
        if os.path.exists(self.output_dir + '/merge.right.fq'):
            os.remove(self.output_dir + '/merge.right.fq')
        os.link(self.merge_left, os.path.join(self.output_dir, 'merge.left.fq'))
        if len(self.right_fq) > 0:
            os.link(self.merge_right, os.path.join(self.output_dir, 'merge.right.fq'))
        self.logger.info("设置结果目录")
        '''

        try:
            if os.path.exists(os.path.join(self.output_dir, 'merge.left.fq')):
                self.option('merge_left', os.path.join(self.output_dir, 'merge.left.fq'))
                if len(self.right_fq) > 0:
                    self.option('merge_right', os.path.join(self.output_dir, 'merge.right.fq'))
        except Exception as e:
            self.logger.info("设置合并结果目录失败{}".format(e))


    def run(self):
        super(FastqMerge2Tool, self).run()
        self.get_filter_samples()
        if self.option("assemble_method") == "total":
            self.run_fastq_merge()
        elif self.option("assemble_method") == "group":
            self.run_fastq_merge_group()
        elif self.option("assemble_method") == "sample":
            self.run_fastq_merge_sample()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'
        data = {
            "id": "FastqMerge2" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.fastq_merge",
            "instant": False,
            "options": dict(
                sample_fq_list="/mnt/ilustre/users/sanger-dev/workspace/20190617/Denovorna_tsg_34448/HiseqQc/output/sickle_dir/fq_list.txt",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
