# -*- coding: utf-8 -*-
# __author__ = 'konghualei,qinjincheng'

from biocluster.agent import Agent
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import re
import json
import time
import unittest

class GeneFaAgent(Agent):
    '''
    last_modify: 2019.07.04
    '''
    def __init__(self, parent):
        super(GeneFaAgent, self).__init__(parent)
        options = [
            {'name': 'ref_new_gtf', 'type': 'infile', 'format': 'ref_rna_v2.gtf'},
            {'name': 'ref_genome_custom', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'transcript_bed', 'type': 'outfile', 'format': 'gene_structure.bed'},
            {'name': 'gene_bed', 'type': 'outfile', 'format': 'gene_structure.bed'},
            {'name': 'gene_fa','type':'outfile', 'format': 'sequence.fasta'}
        ]
        self.add_option(options)
        self.step.add_steps('gene_fa')
        self.on('start',self.step_start)
        self.on('end',self.step_end)

    def step_start(self):
        self.step.gene_fa.start()
        self.step.update()

    def step_end(self):
        self.step.gene_fa.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(os.path.getsize(self.option('ref_genome_custom').path) / 1024 ** 3 + 16)

    @toolfuncdeco
    def end(self):
        super(GeneFaAgent, self).end()

class GeneFaTool(Tool):
    def __init__(self, config):
        super(GeneFaTool, self).__init__(config)
        self.program = {
            'gtftogenepred': 'bioinfo/align/ucsc_tools/gtfToGenePred',
            'genepredtobed': 'bioinfo/align/ucsc_tools/genePredToBed',
            'python': 'miniconda2/bin/python',
            'bedtools': 'miniconda2/bin/bedtools'
        }
        self.script = {
            'get_gene_bed': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/get_gene_bed.py')
        }
        self.file = {
            'txpt_genepred': os.path.join(self.work_dir, 'transcript.genepred'),
            'txpt_bed': os.path.join(self.output_dir, 'transcript.bed'),
            'gene_bed': os.path.join(self.output_dir, 'gene.bed'),
            'gene_fasta': os.path.join(self.output_dir, 'gene.fasta'),
        }

    @toolfuncdeco
    def run(self):
        super(GeneFaTool, self).run()
        self.run_gtftogenepred()
        self.run_genepredtobed()
        self.run_get_gene_bed()
        self.run_bedtools_getfasta()
        self.set_output()
        self.end()
        # ref_new_bed = self.get_gene_bed(self.option("ref_new_gtf"), self.work_dir + "/tmp",
        #                                 assembly_method=self.option("assembly_method"), query_type="gene")
        # ref_trans_bed = self.get_gene_bed(self.option("ref_new_gtf"), self.work_dir + "/tmp_trans",
        #                                 assembly_method=self.option("assembly_method"), query_type="transcript")
        # tmp_trans_bed = os.path.join(self.work_dir, 'ref_new_trans_bed')
        # tmp_new_bed = os.path.join(self.work_dir, 'ref_new_bed')
        # cmd = """sort -r -k1n {} > {}""".format(ref_new_bed, tmp_new_bed)
        # cmd1 = """sort -r -k1n {} > {}""".format(ref_trans_bed, tmp_trans_bed)
        # os.system(cmd)
        # os.system(cmd1)
        # # if os.path.exists(ref_new_bed):
        # #     os.remove(ref_new_bed)
        # #     self.logger.info("删除{}文件成功!".format(ref_new_bed))
        # # if os.path.exists(ref_trans_bed):
        # #     os.remove(ref_trans_bed)
        # #     self.logger.info("删除{}文件成功!".format(ref_trans_bed))
        # self.option("gene_bed").set_path(tmp_new_bed)
        # self.option("transcript_bed").set_path(tmp_trans_bed)
        # self.get_gene_fasta(ref_new_bed, self.option("ref_genome_custom"), self.output_dir,'gene.fa')
        # self.option("gene_fa").set_path(self.output_dir + "/gene.fa")
        # self.end()

    @toolfuncdeco
    def run_gtftogenepred(self):
        cmd = '{} {} {}'.format(
            self.program['gtftogenepred'], self.option('ref_new_gtf').path, self.file['txpt_genepred']
        )
        cmd_name = 'run_gtftogenepred'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def run_genepredtobed(self):
        cmd = '{} {} {}'.format(self.program['genepredtobed'], self.file['txpt_genepred'], self.file['txpt_bed'])
        cmd_name = 'run_genepredtobed'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def run_get_gene_bed(self):
        cmd = '{} {}'.format(self.program['python'], self.script['get_gene_bed'])
        cmd += ' --gtf {}'.format(self.option('ref_new_gtf').path)
        cmd += ' --bed {}'.format(self.file['txpt_bed'])
        cmd += ' --output {}'.format(self.file['gene_bed'])
        cmd_name = 'run_get_gene_bed'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def run_bedtools_getfasta(self):
        cmd = '{} getfasta -name -s'.format(self.program['bedtools'])
        cmd += ' -fi {}'.format(self.option('ref_genome_custom').path)
        cmd += ' -fo {}'.format(self.file['gene_fasta'])
        cmd += ' -bed {}'.format(self.file['gene_bed'])
        cmd_name = 'run_bedtools_getfasta'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def set_output(self):
        self.option('transcript_bed').set_path(self.file['txpt_bed'])
        self.option('gene_bed').set_path(self.file['gene_bed'])
        self.option('gene_fa').set_path(self.file['gene_fasta'])

    def get_gene_bed(self, transcript_gtf = None, gene_bed = None, assembly_method = 'stringtie', query_type = None):

        """
        :param transcript_gtf: stringtie或cufflinks生成的新转录本的gtf文件
        :param gene_bed: 基因的bed文件
        :param assembly_method: 组装方法
        :param query_type: 基因或转录本
        :return:
        """
        # 此函数只适应于stringtie拼接的结果文件,cufflinks
        # cufflinks 组装生成的新基因 'XLOC_xx' 新转录本 'TCONS_xx' 样式
        start = time.time()
        print query_type
        if not os.path.exists(transcript_gtf):
            self.set_error("%s文件不存在", variables = (transcript_gtf), code = "33705601")
        else:
            gene_location_info = {}
            with open(transcript_gtf, 'r+') as f1, open(gene_bed, 'w+') as f2:
                for lines in f1:
                    line = lines.strip().split("\t")
                    if assembly_method == 'stringtie':
                        gene_id_m = re.search(r'gene_id\s+\"(.*?)\";', lines)
                        trans_id_m = re.search(r'transcript_id\s+\"(.*?)\";', lines)
                        if not gene_id_m:
                            gene_id_m = re.search(r'gene_id\s+\"(.*?)\";', lines)
                        if not trans_id_m:
                            trans_id_m = re.search(r'transcript_id\s+\"(.*?)\";', lines)
                    else:
                        gene_id_m = re.search(r'gene_id\s+\"(.*?)\"', lines)
                        trans_id_m = re.search(r'transcript_id\s+\"(.*?)\"', lines)
                    if gene_id_m:
                        gene_id = gene_id_m.group(1)
                    else:
                        continue
                    if trans_id_m:
                        trans_id = trans_id_m.group(1)
                    else:
                        continue
                    if query_type == 'gene':
                        seq_id = gene_id
                    if query_type == 'transcript':
                        seq_id = trans_id
                    start = int(line[3]) - 1  # gtf转为bed格式需要start坐标减去1
                    if seq_id not in gene_location_info.keys():
                        gene_location_info[seq_id] = {}
                        gene_location_info[seq_id]['start'] = start
                        gene_location_info[seq_id]['end'] = int(line[4])
                        gene_location_info[seq_id]['chr'] = line[0]
                        gene_location_info[seq_id]['str'] = line[6]
                    else:
                        #if query_type == 'gene':
                        if gene_location_info[seq_id]['start'] > start:
                            gene_location_info[seq_id]['start'] = start
                        if gene_location_info[seq_id]['end'] <= int(line[4]):
                            gene_location_info[seq_id]['end'] = int(line[4])
                        #else:
                            #continue
                if gene_location_info:
                    for keys, values in gene_location_info.items():
                        info = [values['chr'], values['start'], values['end'], keys, '0', values['str']]
                        info_tmp = [str(i) for i in info]
                        f2.write("\t".join(info_tmp) + "\n")
            end = time.time()
            duration = end - start
            print '提取新基因的bed序列耗时{}s'.format(str(duration))
            return gene_bed

    def get_gene_fasta(self,gene_bed, ref_fa, fa_path, filename):
        """通过getfasta得到基因的fa文件"""
        out_fa = os.path.join(fa_path, filename)
        cmd = "%s getfasta -fi %s -bed %s -fo %s -name -s" % (self.getfasta_path, ref_fa, gene_bed, out_fa)
        self.logger.info("开始打印cmd命令!")
        self.logger.info(cmd)
        fa_cmd = self.add_command("gene_fa",cmd).run()
        self.wait(fa_cmd)
        if fa_cmd.return_code == 0:
            self.logger.info("%s运行完成" % fa_cmd)
        else:
            self.set_error("%s运行出错", variables = (fa_cmd), code = "33705602")
        print "生成基因的fa序列"
        return out_fa

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'gene_fa_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.gene_fa',
            'instant': False,
            'options': {
                'ref_new_gtf': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/gene_fa/ref_and_new.gtf',
                'ref_genome_custom': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
