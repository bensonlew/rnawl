# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class CatFileAgent(Agent):
    def __init__(self, parent):
        super(CatFileAgent, self).__init__(parent)
        options = [
            {'name': 'file_list', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'sample_name', 'type': 'string'},
            {'name': 'gff', 'type': 'outfile', 'format': 'gene_structure.gff3'},
            {'name': 'prot', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'cds', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'rnt', 'type': 'outfile', 'format': 'sequence.profile_table'},
            {'name': 'trna', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'rrna', 'type': 'outfile', 'format': 'sequence.fasta'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '2G'


class CatFileTool(Tool):
    def __init__(self, config):
        super(CatFileTool, self).__init__(config)
        self.sh = '/bin/sh'

    def run(self):
        super(CatFileTool, self).run()
        self.run_cat()
        self.end()

    def run_cat(self):
        outfile = {
                'gff': os.path.join(self.output_dir, self.option('sample_name') + '.gff'),
                'cds': os.path.join(self.output_dir, self.option('sample_name') + '.fnn'),
                'prot': os.path.join(self.output_dir, self.option('sample_name') + '.faa'),
                'rnt': os.path.join(self.output_dir, self.option('sample_name') + '.rnt'),
                'trna': os.path.join(self.output_dir, self.option('sample_name') + '.tRNA.fnn'),
                'rrna': os.path.join(self.output_dir, self.option('sample_name') + '.rRNA.fnn')
                }
        out_opts = []
        if os.listdir(self.output_dir):
            os.system('rm {}/*'.format(self.output_dir))
        with open(self.option('file_list').path, 'r') as r, open('cat.sh', 'w') as w:
            for l in r:
                l = l.strip().split()
                if l[0] in outfile:
                    w.write('cat {} >> {}\n'.format(l[1], outfile[l[0]]))
                    out_opts.append(l[0])
        cmd = self.sh + ' cat.sh'
        command = self.add_command('cat_sh', cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('cat done')
            [self.option(k, outfile[k]) for k in out_opts]
        else:
            self.set_error('wrong in prodigal prediction')

