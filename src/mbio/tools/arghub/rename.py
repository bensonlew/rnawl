# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class RenameAgent(Agent):
    def __init__(self, parent):
        super(RenameAgent, self).__init__(parent)
        options = [
            {'name': 'file_list', 'type': 'infile', 'format': 'arghub.file_list'},
            {'name': 'cds', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'prot', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'gff', 'type': 'infile', 'format': 'gene_structure.gff3'},
            {'name': 'rename_cds', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'rename_prot', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'rename_gff', 'type': 'outfile', 'format': 'gene_structure.gff3'},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('file_list').is_set:
            raise OptionError('找不到输入文件')
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '4G'


class RenameTool(Tool):
    '''
    '''
    def __init__(self, config):
        super(RenameTool, self).__init__(config)

    def run(self):
        super(RenameTool, self).run()
        try:
            self.run_rename()
        except Exception as e:
            self.set_error('wrong in out format {}'.format(e))
        self.end()

    def run_rename(self):
        self.logger.info(self.option('file_list').path)
        self.logger.info(self.option('file_list').file_list)
        print "in run_name"
        with open('names.txt', 'w') as w:
            i = 1
            for f in self.option('file_list').file_list:
                self.logger.info(f)
                with open(f, 'r') as r:
                    for l in r:
                        # new_name \t ori_name
                        w.write('{}\tgene_{:06}\n'.format(l.strip(), i))
                        i += 1

        perl = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/bin/perl'
        cmd1 = perl + " -lane 'if ($ARGV[0]){{$i{{$F[0]}} = $F[1] }}" +\
                " else {{ s/>(\S+)/>$i{{$1}}/; print }}' names.txt {} > {}"
        cmd2 = perl + " -lane 'if ($ARGV[0]){{$i{{$F[0]}} = $F[1]; }}" +\
                "else {{ s/ID\\=([^\;]+)/ID\\=$i{{$1}}/; s/\\t(\S+)/\\t$i{{$1}}/; print }}' names.txt {} > {}"
        command1 = self.add_command('rename_cds', cmd1.format(self.option('cds').path, 'rename_cds'), shell=True).run()
        command2 = self.add_command('rename_prot', cmd1.format(self.option('prot').path, 'rename_prot'), shell=True).run()
        command3 = self.add_command('rename_gff', cmd2.format(self.option('gff').path, 'rename_gff'), shell=True).run()
        self.wait()
        err_c = ''
        for c in [command1, command2, command3]:
            if c.return_code != 0:
                err_c += c.name + ','
            else:
                self.option(c.name, self.work_dir + '/' + c.name)
                self.logger.info('{} 完成'.format(c.name))
        if err_c:
            self.set_error('{} 运行出错'.format(err_c))
