# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class ProdigalAgent(Agent):
    def __init__(self, parent):
        super(ProdigalAgent, self).__init__(parent)
        options = [
            {'name': 'input', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'mode', 'type': 'string', 'default': 'single'},
            {'name': 'gff', 'type': 'outfile', 'format': 'gene_structure.gff3'},
            {'name': 'prot', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'cds', 'type': 'outfile', 'format': 'sequence.fasta'},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('input').is_set:
            raise OptionError('缺少输入序列文件')

    def set_resource(self):
        self._cpu = 1
        self._memory = '2G'


class ProdigalTool(Tool):
    def __init__(self, config):
        super(ProdigalTool, self).__init__(config)
        self.prodigal = 'bioinfo/metaGenomic/Prodigal-2.6.3/prodigal'
        self.genome = self.option('input').prop['path']

    def run(self):
        super(ProdigalTool, self).run()
        self.run_predict()
        self.set_output()
        self.end()

    def run_predict(self):
        cmd = '{0} -i {1} -o output/{2}.gff -f gff -d output/{2}.cds.fna -a output/{2}.prot.faa -p {3}'
        cmd = cmd.format(self.prodigal, self.genome,
                         os.path.basename(self.genome),
                         self.option('mode'))
        command = self.add_command('prodigal_predict', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('prodigal predcit done')
        else:
            self.set_error('wrong in prodigal prediction')

    def set_output(self):
        prefix = os.path.basename(self.genome)
        gff = prefix + '.gff'
        cds = prefix + '.cds.fna'
        prot = prefix + '.prot.faa'
        self.option('gff', os.path.join(self.output_dir, gff))
        self.option('cds', os.path.join(self.output_dir, cds))
        self.option('prot', os.path.join(self.output_dir, prot))

