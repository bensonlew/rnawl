# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modify:2018.11.30

from biocluster.module import Module
import os
import shutil
import re
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import unittest
import glob

class TargetPredictModule(Module):
    '''
    预测靶基因并合并结果
    novol: 预测的新small_rna fasta格式
    known: 已知small_rna fasta格式
    ref: 靶基因序列， 动物3’utr, 植物cds
    method:靶基因预测方法， 以";" 隔开
    min_support: 至少有几个软件支持
    type: animal/plant
    outtable: 输出的合并结果
    anno_detail: 注释详情表
    '''
    def __init__(self, work_id):
        super(TargetPredictModule, self).__init__(work_id)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        options = [
            {"name": "novol", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "known", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "ref", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "method", "type": "string", "default": "miranda"},
            {"name": "type", "type": "string", "default": "animal"},
            {"name": "min_support", "type": "int", "default": 1},
            {"name": "anno_detail", "type": "string", "default": None},
            {"name": "circ_detail", "type": "string", "default": None},
            {"name": "outtable", "type": "outfile", "format": "small_rna.common"},
            {"name": "species", "type": "string", "default": None},
            {"name": "miranda_score", "type": "string", "default": "160"},
            {"name": "miranda_energy", "type": "string", "default": "-20"},
            {"name": "miranda_strict", "type": "string", "default": "on"},
            {"name": "rnahybird_num", "type": "string", "default": "100"},
            {"name": "rnahybird_energy", "type": "string", "default": "-20"},
            {"name": "rnahybird_pvalue", "type": "string", "default": "0.01"},
            {"name": "rnahybrid_num", "type": "string", "default": "100"},
            {"name": "rnahybrid_energy", "type": "string", "default": "-20"},
            {"name": "rnahybrid_pvalue", "type": "string", "default": "0.01"},
            {"name": "ps_robot_score", "type": "string", "default": "2.5"},
            {"name": "targetfinder_score", "type": "string", "default": "4"},
            {"name": "version", "type": "string", "default": "v1"},
        ]
        self.add_option(options)
        self.known_predict_tools = []
        self.novol_predict_tools = []

        self.predict_known = dict()
        self.predict_novol = dict()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        if not self.option("known").is_set:
            raise OptionError("必须设置参数已知small_rnaknown")
        if len(self.option("method").split(",")) < self.option("min_support"):
            raise OptionError('选择的软件数量应该>= 支持的软件数量')
        return True

    def merge_known(self):
        options = {
            'anno_detail': self.option('anno_detail'),
            'min_support': self.option('min_support'),
        }
        if self.option("version") >= "v1.2" and self.option('circ_detail'):
            options.update({'circ_detail': self.option('circ_detail')})
        methods = self.option("method").split(",")
        for n,method in enumerate(methods):
            options.update({
                method.lower(): self.predict_known[method].output_dir + '/' + method.lower() + '_merge_out'
            })
        self.known_merge_tool.set_options(options)
        self.known_merge_tool.run()

    def merge_novol(self):
        options = {
            'anno_detail': self.option('anno_detail'),
            'min_support': self.option('min_support'),
        }
        if self.option("version") >= "v1.2" and self.option('circ_detail'):
            options.update({'circ_detail': self.option('circ_detail')})
        methods = self.option("method").split(",")
        for novol_tools in self.novol_predict_tools:
            self.logger.info("novel_predict tools out {}".format(novol_tools.output_dir))
        for n,method in enumerate(methods):
            options.update({
                method.lower(): self.predict_novol[method].output_dir + '/' + method.lower() + '_merge_out'
            })
        self.novol_merge_tool.set_options(options)
        self.novol_merge_tool.run()

    def run_target(self):
        methods = self.option("method").split(",")
        for method in methods:
            if self.option("version") >= "v1.1":
                predict_tool = self.add_tool('small_rna.mirna_target2')
            else:
                predict_tool = self.add_tool('small_rna.mirna_target')
            options = {
                'target': self.option('ref').prop['path'],
                'query': self.option('known').prop['path'],
                "spiece" : self.option('type'),
                "method": method.lower()
            }
            if method.lower() == "targetscan":
                options.update({
                    "target_species": self.option("species")
                })

            if self.option("version") >= "v1.1":
                if method.lower() == "miranda":
                    options.update({
                        "miranda_score": self.option("miranda_score"),
                        "miranda_energy": self.option("miranda_energy"),
                        "miranda_strict": self.option("miranda_strict"),
                    })
                if method.lower() == "rnahybird":
                    options.update({
                        "rnahybird_num": self.option("rnahybird_num"),
                        "rnahybird_energy": self.option("rnahybird_energy"),
                        "rnahybird_pvalue": self.option("rnahybird_pvalue"),
                    })
                if method.lower() == "rnahybrid":
                    options.update({
                        "rnahybrid_num": self.option("rnahybrid_num"),
                        "rnahybrid_energy": self.option("rnahybrid_energy"),
                        "rnahybrid_pvalue": self.option("rnahybrid_pvalue"),
                    })

                if method.lower() == "ps_robot":
                    options.update({
                        "ps_robot_score": self.option("ps_robot_score"),
                    })
                if method.lower() == "targetfinder":
                    options.update({
                        "targetfinder_score": self.option("targetfinder_score"),
                    })
            else:
                pass

            predict_tool.set_options(options)
            self.predict_known[method] = predict_tool
            self.known_predict_tools.append(predict_tool)
            if self.option("novol").is_set:
                if self.option("version") >= "v1.1":
                    predict_tool = self.add_tool('small_rna.mirna_target2')
                else:
                    predict_tool = self.add_tool('small_rna.mirna_target')
                options = {
                    'target': self.option('ref').prop['path'],
                    'query': self.option('novol').prop['path'],
                    "spiece" : self.option('type'),
                    "method": method.lower()
                }
                if method.lower() == "targetscan":
                    options.update({
                        "target_species": self.option("species")
                    })

                if self.option("version") >= "v1.1":
                    if method.lower() == "miranda":
                        options.update({
                            "miranda_score": self.option("miranda_score"),
                            "miranda_energy": self.option("miranda_energy"),
                            "miranda_strict": self.option("miranda_strict"),
                        })
                    if method.lower() == "rnahybird":
                        options.update({
                            "rnahybird_num": self.option("rnahybird_num"),
                            "rnahybird_energy": self.option("rnahybird_energy"),
                            "rnahybird_pvalue": self.option("rnahybird_pvalue"),
                        })
                    if method.lower() == "rnahybrid":
                        options.update({
                            "rnahybrid_num": self.option("rnahybrid_num"),
                            "rnahybrid_energy": self.option("rnahybrid_energy"),
                            "rnahybrid_pvalue": self.option("rnahybrid_pvalue"),
                        })

                    if method.lower() == "ps_robot":
                        options.update({
                            "ps_robot_score": self.option("ps_robot_score"),
                        })
                    if method.lower() == "targetfinder":
                        options.update({
                            "targetfinder_score": self.option("targetfinder_score"),
                        })
                else:
                    pass


                predict_tool.set_options(options)
                self.predict_novol[method] = predict_tool
                self.novol_predict_tools.append(predict_tool)
        for novol_tools in self.novol_predict_tools:
            self.logger.info("novel_predict tools out {}".format(novol_tools.output_dir))

        self.on_rely(self.known_predict_tools, self.merge_known)
        if self.option("novol").is_set:
            self.on_rely(self.novol_predict_tools, self.merge_novol)
            self.on_rely([self.known_merge_tool, self.novol_merge_tool], self.set_output)
        else:
            self.known_merge_tool.on('end', self.set_output)
        for tool in self.known_predict_tools + self.novol_predict_tools:
            tool.run()

    def import_target(self, target, target_new=None):
        '''
        导入 靶基因信息
        '''
        target_smallrna = dict()
        with open(target, 'rb') as target_f:
            target_f.readline()
            for line in target_f:
                if line.startswith("#"):
                    pass
                else:
                    mirna, target = line.strip().split("\t")[:2]
                    if target in target_smallrna:
                        target_smallrna[target].append(mirna)
                    else:
                        target_smallrna[target] = [mirna]
        if target_new:
            with open(target_new, 'rb') as target_f:
                target_f.readline()
                for line in target_f:
                    if line.startswith("#"):
                        pass
                    else:
                        mirna, target = line.strip().split("\t")[:2]
                        if target in target_smallrna:
                            target_smallrna[target].append(mirna)
                        else:
                            target_smallrna[target] = [mirna]
        return target_smallrna

    def merge_annot_smallrna(self):
        if self.option('novol').is_set:
            target_smallrna = self.import_target(self.known_merge_tool.output_dir + '/all_merged.xls', self.novol_merge_tool.output_dir + '/all_merged.xls')
        else:
            target_smallrna = self.import_target(self.known_merge_tool.output_dir + '/all_merged.xls', None)
        with open(self.option("anno_detail"), 'r') as annot_in, open(self.output_dir + "/All_annot_target.xls", 'w') as annot_out:
            header = annot_in.readline()
            cols = header.split("\t")

            annot_type = "ref"
            if header.startswith("transcript"):
                annot_type = "denovo"
            annot_out.write("target_transcript\tsmall_rnas\tgene\t" + "\t".join(cols[3:]))
            for line in annot_in:
                cols = line.split("\t")
                if annot_type == "denovo":
                    if cols[0] in target_smallrna:
                        small_rnas = ";".join(target_smallrna[cols[0]])
                    else:
                        small_rnas = ""
                else:
                    if cols[1] in target_smallrna:
                        small_rnas = ";".join(target_smallrna[cols[1]])
                    else:
                        small_rnas = ""
                annot_out.write(cols[1] + "\t" + small_rnas + "\t" + cols[0] + "\t" + "\t".join(cols[3:]))

    def set_output(self):
        self.merge_annot_smallrna()
        if os.path.exists(self.output_dir + '/known_target.xls'):
            os.remove(self.output_dir + '/known_target.xls')
        if os.path.exists(self.output_dir + '/novol_target.xls'):
            os.remove(self.output_dir + '/novol_target.xls')
        if os.path.exists(self.output_dir + '/target.fa'):
            os.remove(self.output_dir + '/target.fa')
        os.link(self.option("ref").prop['path'], self.output_dir + '/target.fa')
        for tool in self.known_predict_tools:
            detail_file = glob.glob(tool.output_dir + "/*detail.txt.gz")
            if detail_file:
                detail_out = self.output_dir + '/known_' + os.path.basename(detail_file[0])
                if os.path.exists(detail_out):
                    os.remove(detail_out)
                os.link(detail_file[0], detail_out)
        if self.option('novol').is_set:
            os.link(self.novol_merge_tool.output_dir + '/all_merged.xls', self.output_dir + '/novol_target.xls')
            for tool in self.novol_predict_tools:
                detail_file = glob.glob(tool.output_dir + "/*detail.txt.gz")
                if detail_file:
                    detail_out = self.output_dir + '/novol_' + os.path.basename(detail_file[0])
                    if os.path.exists(detail_out):
                        os.remove(detail_out)
                    os.link(detail_file[0], detail_out)
        os.link(self.known_merge_tool.output_dir + '/all_merged.xls', self.output_dir + '/known_target.xls')

        self.end()

    def run(self):
        super(TargetPredictModule, self).run()
        if self.option("version") >= "v1.2" and self.option("circ_detail"):
            self.known_merge_tool = self.add_tool('small_rna.merge_target_result3')
            self.novol_merge_tool = self.add_tool('small_rna.merge_target_result3')
        elif self.option("version") >= "v1.1":
            self.known_merge_tool = self.add_tool('small_rna.merge_target_result2')
            self.novol_merge_tool = self.add_tool('small_rna.merge_target_result2')
        else:
            self.known_merge_tool = self.add_tool('small_rna.merge_target_result')
            self.novol_merge_tool = self.add_tool('small_rna.merge_target_result')
        self.run_target()

    def end(self):
        repaths = [
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(TargetPredictModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the module Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/'
        data = {
            "id": "target_module" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "module",
            "name": "small_rna.target_predict",
            "instant": False,
            "options": dict(
                novol = test_dir + 'mature_rno_novol.dna.fa',
                known = test_dir + 'mature_rno.dna.fa',
                ref = test_dir + 'animal.exon.dna.fa',
                method = 'miranda;targetscan;rnahybrid',
                anno_detail='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Rattus_norvegicus/Ensemble_release_89/Annotation_v2/annot_class/anno_stat/all_anno_detail.xls',
                type = 'animal',
                species = 'Rattus_norvegicus'
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
