# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import os
import unittest
__author__ = 'fengyitong'


class MergeTargetResultAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(MergeTargetResultAgent, self).__init__(parent)
        options = [
            {'name': 'anno_detail', 'type': 'string', 'default': ''},
            {'name': 'targetscan', 'type': 'string', 'default': ''},
            {'name': 'miranda', 'type': 'string', 'default': ''},
            {'name': 'rnahybrid', 'type': 'string', 'default': ''},
            {'name': 'psrobot', 'type': 'string', 'default': ''},
            {'name': 'targetfinder', 'type': 'string', 'default': ''},
            {'name': 'min_support', 'type': 'int', 'default': 0},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('anno_detail') and not os.path.exists(self.option('anno_detail')):
            raise OptionError("不能没有注释文件")
        if self.option('targetscan') and not os.path.exists(self.option('targetscan')):
            raise OptionError("有传targetscan输出文件的参数，但是文件{}不存在".format(self.option('targetscan')))
        if self.option('miranda') and not os.path.exists(self.option('miranda')):
            raise OptionError("有传miranda输出文件的参数，但是文件不存在")
        if self.option('rnahybrid') and not os.path.exists(self.option('rnahybrid')):
            raise OptionError("有传rnahybrid输出文件的参数，但是文件不存在")
        if self.option('psrobot') and not os.path.exists(self.option('psrobot')):
            raise OptionError("有传psrobot输出文件的参数，但是文件不存在")
        if self.option('targetfinder') and not os.path.exists(self.option('targetfinder')):
            raise OptionError("有传targetfinder输出文件的参数，但是文件不存在")

    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('15')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        super(MergeTargetResultAgent, self).end()

class MergeTargetResultTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(MergeTargetResultTool, self).__init__(config)
        self.result_summary = defaultdict(dict)
        self.mirtarbase = self.config.SOFTWARE_DIR + '/database/miRTarBase/miRTarBase_MTI.txt'

    def get_targetscan(self):
        couple2info = dict()
        with open(self.option('targetscan'), 'r') as f_r:
            _ = f_r.readline()
            for line in f_r.readlines():
                line = line.strip().split('\t')
                if len(line) > 4:
                    couple2info['%s|%s' % (line[0], line[1])] = '\t'.join([line[2], line[3], line[4], line[5], line[6]])
        return couple2info

    def get_miranda(self):
        couple2info = dict()
        with open(self.option('miranda'), 'r') as f_r:
            _ = f_r.readline()
            for line in f_r.readlines():
                line = line.strip().split('\t')
                if len(line) > 4:
                    couple2info['%s|%s'%(line[0], line[1])] = '\t'.join([line[6], line[7], line[2], line[3]])
            return couple2info

    def get_rnahybrid(self):
        couple2info = dict()
        with open(self.option('rnahybrid'), 'r') as f_r:
            _ = f_r.readline()
            for line in f_r.readlines():
                line = line.strip().split('\t')
                if len(line) > 4:
                    couple2info['%s|%s'%(line[1], line[0])] = '\t'.join([line[2], line[3], line[4]])
            return couple2info

    def get_psRobot(self):
        couple2info = dict()
        with open(self.option('psrobot'), 'r') as f_r:
            _ = f_r.readline()
            for line in f_r.readlines():
                line = line.strip().split('\t')
                if len(line) > 4:
                    couple2info['%s|%s'%(line[0], line[2])] = '\t'.join([line[5], line[6], line[1], "None"])
            return couple2info

    def get_targetfinder(self):
        couple2info = dict()
        with open(self.option('targetfinder'), 'r') as f_r:
            _ = f_r.readline()
            for line in f_r.readlines():
                line = line.strip().split('\t')
                if len(line) > 4:
                    couple2info['%s|%s'%(line[0], line[1])] = '\t'.join([line[2], line[3], line[4]])
            return couple2info

    def get_mirtarbase_set(self):
        mir_set = set()
        with open(self.mirtarbase, 'r') as base_r:
            _ = base_r.readline()
            for line in base_r.readlines():
                line = line.strip().split('\t')
                if len(line) > 4:
                    mir_set.add('%s|%s'%(line[1], line[3].lower()))
        return mir_set

    def get_anno_info(self):
        anno_info = dict()
        g2gn = defaultdict(set)
        with open(self.option('anno_detail'), 'r') as f_r:
            annot_type = "ref"
            header = f_r.readline()
            if header.startswith("transcript"):
                annot_type = "denovo"

            for line in f_r.readlines():
                line = line.strip("\n").split('\t')
                if len(line) > 5:
                    if annot_type == "denovo":
                        anno_info[line[0]] = '\t'.join([line[1], "-", line[12]])
                        g2gn[line[0]].add("")
                    else:
                        if line[3] == "":
                            line[3] = "-"
                        anno_info[line[1]] = '\t'.join([line[0], line[3], line[5]])
                        g2gn[line[1]].add(line[3].lower())
        return anno_info, g2gn

    def run(self):
        super(MergeTargetResultTool, self).run()
        mir_set = self.get_mirtarbase_set()
        anno_info, g2gn = self.get_anno_info()
        header = 'query\ttarget\tgene_id\tgene_name\tgene_description\t'
        if self.option('miranda'):
            self.result_summary['miranda'] = self.get_miranda()
            header += 'start_miranda\tend_miranda\tscore\tenergy\t'
        if self.option('targetscan'):
            self.result_summary['targetscan'] = self.get_targetscan()
            header += 'utr_start_targetscan\tutr_end_targetscan\tmsa_start_targetscan\tmsa_end_targetscan\tseed_match_targetscan\t'
        if self.option('psrobot'):
            self.result_summary['psrobot'] = self.get_psRobot()
            header += 'start_psrobot\tend_psrobot\tscore\tenergy\t'
        if self.option('targetfinder'):
            self.result_summary['targetfinder'] = self.get_targetfinder()
            header += 'start_targetfinder\tend_targetfinder\tscore_targetfinder\t'
        if self.option('rnahybrid'):
            self.result_summary['rnahybrid'] = self.get_rnahybrid()
            header += 'start_rnahybrid\tend_rnahybrid\tenergy_rnahybrid\t'
        header = header + 'mirtarbase\n'
        all_matched = set()
        for key in self.result_summary:
            if not all_matched:
                all_matched = set(self.result_summary[key].keys())
            else:
                all_matched = all_matched | set(self.result_summary[key].keys())
        mir_gene_set = set()
        with open(os.path.join(self.output_dir, 'all_merged.xls'), 'w') as all_w:
            all_w.write(header)
            for match in all_matched:
                line = match.split('|')[0] + '\t' + match.split('|')[1] + '\t'
                if match.split('|')[1] in anno_info:
                    if anno_info[match.split('|')[1]][1] == "":
                        anno_info[match.split('|')[1]][1] == "-"
                    line += anno_info[match.split('|')[1]] + '\t'
                else:
                    line += '_\t_\t_\t'
                support = 0
                if self.option('miranda'):
                    if match in self.result_summary['miranda']:
                        support += 1
                        line += self.result_summary['miranda'][match] + '\t'
                    else:
                        line += '\t\t\t\t'
                if self.option('targetscan'):
                    if match in self.result_summary['targetscan']:
                        support += 1
                        line += self.result_summary['targetscan'][match] + '\t'
                    else:
                        line += '\t\t\t\t\t'
                if self.option('psrobot'):
                    if match in self.result_summary['psrobot']:
                        support += 1
                        line += self.result_summary['psrobot'][match] + '\t'
                    else:
                        line += '\t\t\t\t'
                if self.option('targetfinder'):
                    if match in self.result_summary['targetfinder']:
                        support += 1
                        line += self.result_summary['targetfinder'][match] + '\t'
                    else:
                        line += '\t\t\t'
                if self.option('rnahybrid'):
                    if match in self.result_summary['rnahybrid']:
                        support += 1
                        line += self.result_summary['rnahybrid'][match] + '\t'
                    else:
                        line += '\t\t\t'

                if match.split('|')[1] in g2gn:
                    for gn in g2gn[match.split('|')[1]]:
                        judge = False
                        if '%s|%s'%(match.split('|')[0], gn) in mir_set:
                            print('%s|%s'%(match.split('|')[0], gn))
                            line += 'yes\n'
                            judge = True
                            break
                    else:
                        if not judge:
                            line += 'no\n'
                else:
                    line += 'no\n'
                if support >= self.option("min_support"):
                    if line.split("\t")[0] + "|" + line.split("\t")[2] in mir_gene_set:
                        pass
                    else:
                        mir_gene_set.add(line.split("\t")[0] + "|" + line.split("\t")[2])
                        all_w.write(line)
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/'
        data = {
            "id": "MergeTargetResult" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "small_rna.merge_target_result",
            "instant": False,
            "options": dict(
                # targetscan=test_dir + "targetscan_merge_out",
                miranda=test_dir + "miranda_merge_out",
                # rnahybrid=test_dir + "rnahybrid_merge_out",
                rnahybrid=test_dir + "RNAhybrid_merge_out_animal",
                # psrobot=test_dir + "psrobot_merge_out",
                # targetfinder=test_dir + "targetfinder_merge_out",
                anno_detail='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Rattus_norvegicus/Ensemble_release_89/Annotation_v2/annot_class/anno_stat/all_anno_detail.xls'
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
