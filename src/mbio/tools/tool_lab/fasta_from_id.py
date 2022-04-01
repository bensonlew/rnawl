# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import re


class FastaFromIdAgent(Agent):
    def __init__(self, parent):
        super(FastaFromIdAgent, self).__init__(parent)
        options = [
            {'name': 'in_seq', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'id_list', 'type': 'string'},
            {'name': 'out_seq', 'type': 'string'}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('id_list'):
            raise OptionError('缺少输入id list 文件')
        if not self.option('in_seq').is_set:
            raise OptionError('缺少输入in seq 文件')

    def set_resource(self):
        self._cpu = 1
        self._memory = '1G'


class FastaFromIdTool(Tool):
    def __init__(self, config):
        super(FastaFromIdTool, self).__init__(config)
        self.seq = self.option('in_seq').prop['path']

    def run(self):
        super(FastaFromIdTool, self).run()
        self.do_it()
        self.set_output()
        self.end()

    def do_it(self):
        with open(self.option('id_list'), 'r') as r:
            ids = set(map(lambda x: x.strip(), r.readlines()))
        out_put = os.path.join(self.output_dir, self.option('out_seq'))
        with open(self.option('in_seq').prop['path'], 'r') as in_seq, open(out_put, 'w') as out:
            tmp_id = ''
            ids.add(tmp_id)
            for seq in in_seq:
                if seq.startswith('>'):
                    # ids.remove(tmp_id)
                    tmp_id = re.match(r'>(\S+)', seq).groups()[0]
                if tmp_id in ids:
                    out.write(seq)

    def set_output(self):
        pass
