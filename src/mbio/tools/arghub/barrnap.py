# -*- coding: utf-8 -*-

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class BarrnapAgent(Agent):
    """
    Barrnap 进行rRNA预测
    """

    def __init__(self, parent):
        super(BarrnapAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "kingdom", "type": "string", "default": "bac"},  # genome type Kingdom: arc mito bac euk
            {"name": "lencutoff", "type": "string", "default": "0.8"},
            # Proportional length threshold to label as partial，小于该阈值预测结果会加上(partial)
            {"name": "reject", "type": "string", "default": "0.25"},
            # Proportional length threshold to reject prediction， 只有当大于该阈值的结果才会保留
            {"name": "evalue", "type": "float", "default": 1e-06},  # Similarity e-value cut-off
            {"name": "rrna", "type": "outfile", "format": "sequence.fasta"},  # 预测rRNA序列
            {"name": "r_gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式,并增加两列位置信息，用与绘制圈图
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome", code="33300201")

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(BarrnapAgent, self).end()


class BarrnapTool(Tool):
    def __init__(self, config):
        super(BarrnapTool, self).__init__(config)
        self.genome = self.option("input_genome").prop['path']
        self.barrnap_path = self.config.SOFTWARE_DIR +"/bioinfo/meta/arghub/ariba/bin/"
        self.set_environ(PATH='{0}:{0}'.format(self.barrnap_path))

    def run(self):
        super(BarrnapTool, self).run()
        self.run_barrnap()
        self.end()

    def run_barrnap(self):
        cmd = "{} --threads 2 --kingdom {} --quiet {} --lencutoff {}" +\
                " --reject {} --evalue {} --outseq output/{}.rRNA.fna > output/{}.rRNA.gff"
        cmd = cmd.format(
            self.barrnap_path + 'barrnap',
            self.option("kingdom"), self.genome,
            self.option("lencutoff"), self.option("reject"),
            self.option("evalue"), os.path.basename(self.genome),
            os.path.basename(self.genome))
        command = self.add_command('barrnap', cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('rRNA预测完成')
            self.logger.info('开始设置输出文件')
            self.set_output()
            self.logger.info('成功设置输出文件')

            
    def set_output(self):
        prefix = os.path.basename(self.genome)
        if os.path.getsize(os.path.join(self.output_dir, prefix + '.rRNA.fna')) == 0:
            return 0
        self.option('rrna', os.path.join(self.output_dir, prefix + '.rRNA.fna'))
        self.option('r_gff', os.path.join(self.output_dir, prefix + '.rRNA.gff'))
        self.format_out()

    def format_out(self):
        partner = re.compile(r'Name\=([^\;]+)')
        rna_tmp = self.option('rrna').path + '.tmp'
        gff_tmp = self.option('r_gff').path + '.tmp'
        ids_dict = {}
        loc_dict = {}
        with open(gff_tmp, 'w') as w, open(self.option('r_gff').path, 'r') as r:
            i = 0
            for line in r:
                print line
                match = partner.search(line)
                print 'match: {}'.format(match)
                if match:
                    name = match.groups()[0]
                    print match
                    li = line.split('\t')
                    s_id = '{}::{}:{}-{}({})'.format(name, li[0], int(li[3]) - 1, li[4], li[6])
                    if s_id not in ids_dict:
                        i += 1
                        ids_dict[s_id] = name + '_rrna' + str(i)
                    print line
                    line = partner.sub('ID={};Name={}'.format(ids_dict[s_id], name), line)
                    print line
                    print '######' * 5
                w.write(line)

        self.logger.info('{}'.format(ids_dict.keys()[1]))
        partner2 = re.compile(r'>(\S+)')
        with open(rna_tmp, 'w') as w, open(self.option('rrna').path, 'r') as r:
            for line in r:
                match = partner2.search(line)
                if match:
                    line = partner2.sub('>' + ids_dict[match.groups()[0]], line)
                w.write(line)
        os.remove(self.option('r_gff').path)
        os.rename(gff_tmp, self.option('r_gff').path)
        os.remove(self.option('rrna').path)
        os.rename(rna_tmp, self.option('rrna').path)

