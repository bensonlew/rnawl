# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna.trans_step import step_count


class MetagenemarkAgent(Agent):
    def __init__(self, parent):
        super(MetagenemarkAgent, self).__init__(parent)
        options = [
            {'name': 'input_fasta', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'sample_name', 'type': 'string', 'default': ''},
            {'name': 'mode', 'type': 'string', 'default': 'meta'},
            {'name': 'min_gene', 'type': 'int', 'default': 100},
            {'name': 'sample_info', 'type': 'string', 'default': ''},
            {'name': 'gff', 'type': 'outfile', 'format': 'gene_structure.gff3'},
            {'name': 'prot', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'cds', 'type': 'outfile', 'format': 'sequence.fasta'},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('input_fasta').is_set:
            raise OptionError('缺少输入序列文件')

    def set_resource(self):
        self._cpu = 4
        self._memory = '4G'


class MetagenemarkTool(Tool):
    def __init__(self, config):
        super(MetagenemarkTool, self).__init__(config)
        self.metagenemark = '/bioinfo/metaGenomic/MetaGeneMark_linux_64/mgm/gmhmmp'
        self.mod_file = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod'
        self.set_environ(HOME=os.path.abspath(self.config.SOFTWARE_DIR + '/../'))
        self.genome = self.option('input_fasta').prop['path']
        self.sp_name = self.option('sample_name')
        if not self.sp_name:
            self.sp_name = os.path.basename(self.genome).lsplit('.')[0]

    def run(self):
        super(MetagenemarkTool, self).run()
        self.run_predict()
        self.cut_min()
        self.end()

    def run_predict(self):
        cmd = '{0} {1} -o {2}.gff -f 3 -D {2}.fna -A {2}.faa -m {3}'
        cmd = cmd.format(self.metagenemark, self.genome,
                         self.sp_name, self.mod_file)
        command = self.add_command('metagenemark_predict', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('metagenemark predcit done')
            self.reanme()
        else:
            self.set_error('wrong in metagenemark prediction')

    def reanme(self):
        cmd = '''/usr/bin/perl -i -lpe '$n = "{0}_";s/>(\\w+).*/>$n$1/' {0}.fna'''.format(self.sp_name)
        cmd2 = '''/usr/bin/perl -i -lpe '$n = "{0}_";s/>(\\w+).*/>$n$1\\_1/' {0}.faa'''.format(self.sp_name)
        command = self.add_command('rename_fna', cmd, shell=True).run()
        command2 = self.add_command('rename_faa', cmd2, shell=True).run()
        self.wait()
        if command.return_code == 0 and command2.return_code == 0:
            self.logger.info('rename_done')

    def cut_min(self):
        cmd = "perl -lne '$id = $1 and next if />(\\S+)/;$seq{{$id}} .= $_;" +\
            "END{{ map {{ print \">$_\\n$seq{{$_}}\" if length($seq{{$_}}) >= {0};}} keys %seq}}'" +\
            " {1}.{2} > output/{1}.metagene.more.{0}.{2}"
        cmd1 = cmd.format(self.option('min_gene'), self.sp_name, "fna")
        cmd2 = cmd.format(int(self.option('min_gene') / 3), self.sp_name, "faa")
        print(cmd1)
        print(cmd2)
        fa_path = '{}/{}.metagene.more.{}.fna'.format(self.output_dir, self.sp_name, self.option('min_gene'))
        faa_path = '{}/{}.metagene.more.{}.faa'.format(self.output_dir, self.sp_name,  int(self.option("min_gene") / 3))
        if os.system(cmd1) == 0 and os.system(cmd2) == 0:
            self.logger.info('cut min gene done')
            self.option('prot', faa_path)
            self.option('cds', fa_path)
        else:
            self.set_error('cut min failed')
