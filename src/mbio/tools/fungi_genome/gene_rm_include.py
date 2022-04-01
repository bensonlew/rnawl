# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# version 1.0
# last_modify: 2018.06.013

import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from operator import itemgetter


class GeneRmIncludeAgent(Agent):
    """
    去除被含在另一基因内部的基因
    """

    def __init__(self, parent):
        super(GeneRmIncludeAgent, self).__init__(parent)
        options = [
            {"name": "gff", "type": "string", "default":""},
            {"name": "faa", "type": "string", "default":""},
            {"name": "sample", "type": "string","default":"sample"}
        ]

        self.add_option(options)

    def check_options(self):

        if not self.option("gff"):
            raise OptionError("必须设置参gff", code="32101701")
        if not self.option("faa"):
            raise OptionError("必须设置参数faa", code="32101702")



    def set_resource(self):
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        super(GeneRmIncludeAgent, self).end()


class GeneRmIncludeTool(Tool):
    def __init__(self, config):
        super(GeneRmIncludeTool, self).__init__(config)


    def run_rm_include_faa(self):
        gene_site = []
        len_scf = len('scaffold')
        with open(self.option('gff')) as f:
            for line in f:
                spl = line.split('\t')
                if len(spl)<8:
                    continue
                if spl[1] == 'maker':
                    if spl[2] == 'gene' :
                        gene_name = spl[8].split(';')[0][3:]
                        scaf_id = int(spl[0][len_scf:])
                        gene_site.append([gene_name, scaf_id, int(spl[3]), int(spl[4])])
        s2dim = sorted(gene_site, key=itemgetter(1,2))


        rm_gene_list = []
        ref = ['',-1,-1,-1]
        for i in s2dim:
            if i[1] == ref[1]:
                if i[2] >= ref[3]:
                    ref = i
                    continue

                if i[3] < ref[3]:
                    rm_gene_list.append(i[0])
                    continue
                rep = (i[3] - ref[2])*1.0
                len_i = i[3]-i[2]
                len_ref = ref[3]-ref[2]
                if rep/min(len_i,len_ref) < 0.3 :
                    ref = i
                    continue
                if len_i > len_ref:
                    rm_gene_list.append(ref[0])
                    ref = i
                else:
                    rm_gene_list.append(i[0])

            else:
                ref = i

        rmfile = open(self.work_dir + "/rmlist",'w')
        rmfile.write('\n'.join(rm_gene_list))

        with open(self.option('faa')) as f :
            fw = open(self.work_dir + '/rm_include.faa', 'w')
            len_end = len('-mRNA-1')
            mk_write = 1
            for line in f:
                if line.startswith('>'):
                    mk_write = 1
                    g_name = line.split(' ')[0][1:-len_end]
                    if g_name in rm_gene_list :
                        mk_write = 0
                if mk_write == 1 :
                    fw.write(line)


    def set_output(self):
        old_path = self.work_dir + '/rm_include.faa'
        new_path = self.output_dir + '/{}_rm_include.faa'.format(self.option("sample"))
        if os.path.exists(new_path):
            os.remove(new_path)
        os.link(old_path, new_path)


    def run(self):
        super(GeneRmIncludeTool, self).run()
        self.run_rm_include_faa()
        self.set_output()
        self.end()

