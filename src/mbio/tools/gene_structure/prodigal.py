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


class ProdigalTool(Tool):
    def __init__(self, config):
        super(ProdigalTool, self).__init__(config)
        self.prodigal = 'bioinfo/metaGenomic/Prodigal-2.6.3/prodigal'
        self.genome = self.option('input_fasta').prop['path']
        self.sp_name = self.option('sample_name')
        if not self.sp_name:
            self.sp_name = os.path.basename(self.genome).lsplit('.')[0]

    def run(self):
        super(ProdigalTool, self).run()
        self.run_predict()
        self.cut_min()
        self.set_output()
        self.end()

    def check_cmd(self, cmd_status):
        return os.path.exists(cmd_status)

    def run_predict(self):
        cmd = '{0} -i {1} -o {2}.gff -f gff -d {2}.fna -a {2}.faa -p {3}'
        cmd = cmd.format(self.prodigal, self.genome,
                         self.sp_name,
                         self.option('mode'))
        if self.check_cmd("prodigal_predict.status"):
            return
        command = self.add_command('prodigal_predict', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('prodigal predcit done')
            with open("prodigal_predict.status", 'w') as w:
                w.write('True')
            self.rename()
        else:
            self.set_error('wrong in prodigal prediction')

    def rename(self):
        cmd = '''/usr/bin/perl -i -lpe 's/>(\\S+)/>{0}_$1/' {0}.fna'''.format(self.sp_name)
        cmd2 = '''/usr/bin/perl -i -lpe 's/>(\\S+)/>{0}_$1\\_1/' {0}.faa'''.format(self.sp_name)
        command = self.add_command('rename_fa', cmd, shell=True).run()
        command2 = self.add_command('rename_faa', cmd2, shell=True).run()
        self.wait()
        if command.return_code == 0 and command2.return_code == 0:
            self.logger.info('renam_done')

    def cut_min(self):
        #cmd = "perl -lne '$id = $1 and next if />(\\S+)/;$seq{{$id}} .= $_;" +\
        #    "END{{ map {{ print \">$_\\n$seq{{$_}}\" if length($seq{{$_}}) >= {0};}} keys %seq}}'" +\
        #    " {1}.{2} > output/{1}.metagene.more.{0}.{2}"
        cmd = "perl -lne 'print \">$id\\n$seq\" if />/ and (length($seq) >= {0}); $seq = \"\" if />/;$id = $1 and next if />(\\S+)/;$seq .= $_;"+\
            "END{{ print \">$id\\n$seq\" if length($seq) >= {0} }}'" +\
            " {1}.{2} > output/{1}.metagene.more.{0}.{2}"
        cmd1 = cmd.format(self.option('min_gene'), self.sp_name, "fna")
        cmd2 = cmd.format(int(self.option('min_gene') / 3), self.sp_name, "faa")
        print(cmd1)
        print(cmd2)
        if os.system(cmd1) == 0 and os.system(cmd2) == 0:
            self.logger.info('cut min gene done')
        else:
            self.set_error('cut min failed')

    def set_output(self):
        prefix = self.sp_name + ".metagene.more"
        self.option("cds",
                    self.output_dir + "/{}.{}.fna".format(prefix, self.option("min_gene")))
        self.option("prot",
                    self.output_dir + "/{}.{}.faa".format(prefix, int(self.option("min_gene") / 3)))
