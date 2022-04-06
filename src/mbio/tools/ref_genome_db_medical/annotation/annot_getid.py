# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest


class AnnotGetidAgent(Agent):
    '''
    last_modify: 2020.08.10
    '''
    def __init__(self, parent):
        super(AnnotGetidAgent, self).__init__(parent)
        options = [
            {"name": "g2t2p", "type": "string", "default": None},
            {"name": "species", "type": "string", "default": None},
            {"name": "get_lines", "type": "string", "default": None}, #筛选列表
            {"name": "match_lines", "type": "string", "default": None},
            {"name": "genedb_vesion", "type": "string", "default": None},
            {"name": "ids_out", "type": "string", "default": "ids.tsv"},
            {"name": "ids_out_file", "type": "outfile", "format": "ref_rna_v2.common"},

        ]
        self.add_option(options)
        self.step.add_steps('annot_begin')
        self.on('start', self.step_start)
        self.on('end', self.step_end)


    def step_start(self):
        self.step.annot_begin.start()
        self.step.update()

    def step_end(self):
        self.step.annot_begin.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

        for opt in ['g2t2p', 'species', 'get_lines', 'match_lines']:
            if not self.option(opt):
                self.set_error("参数 {} 需指定".format(opt))

    def set_resource(self):
        self._cpu = 1
        self._memory = '16G'

    def end(self):
        super(AnnotGetidAgent, self).end()

class AnnotGetidTool(Tool):
    def __init__(self, config):
        super(AnnotGetidTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.gffread = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/cufflinks-2.2.1/gffread')
        self.getid = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/getid.py')
        self.id_mapping_dict = dict({
            "Homo_sapiens": self.config.SOFTWARE_DIR + "/database/gene_db/hsapiens.all_merged.xls",
            "Mus_musculus": self.config.SOFTWARE_DIR + "/database/gene_db/mmusculus.all_merged.xls",
            "Rattus_norvegicus": self.config.SOFTWARE_DIR + "/database/gene_db/rnorvegicus.all_merged.xls"
        })



    def run(self):
        super(AnnotGetidTool, self).run()
        self.run_getid()
        self.set_output()
        self.end()

    def run_getid(self):
        cmd = '{} {}'.format(self.python, self.getid)
        cmd += ' {}'.format(self.option("g2t2p"))
        cmd += ' {}'.format(self.id_mapping_dict[self.option('species')])
        cmd += ' {}'.format(self.option("match_lines"))
        cmd += ' {}'.format(self.option("get_lines"))
        cmd += ' {}'.format(self.option("ids_out"))
        cmd_name = 'run_get_id'
        self.run_code(cmd_name, cmd)


    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run %s, abord', variables=(n), code="33710102")

    def set_output(self):
        self.option("ids_out_file", self.work_dir + '/' + self.option("ids_out"))
        outfiles = ["ids.tsv"]
        for item in outfiles:
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + item, linkfile)
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        getline_list = [
            "ensembl_gene_id",
            "external_gene_name",
            "ensembl_transcript_id",
            "external_transcript_name",
            "gene_biotype",
            "transcript_biotype",
            "uniprot_gn_id",
            "uniprotswissprot",
            "reactome",
            "kegg",
            "eggnog",
            "entrezgene_id",
            "refseq_mrna",
            "refseq_peptide",
            "refseq_ncrna",
            "refseq_peptide_predicted",
            "refseq_ncrna_predicted",
            "refseq_mrna_predicted"
        ]
        matchline_list = [
            "ensembl_transcript_id",
            "refseq_mrna",
            "refseq_ncrna",
            "refseq_ncrna_predicted",
            "refseq_mrna_predicted"
        ]
        data = {
            'id': 'annot_getid_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_getid',
            'instant': False,
            'options': dict(
                g2t2p = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/Annotation_v2/g2t2p.txt",
                species = "Homo_sapiens",
                get_lines = ",".join(getline_list),
                match_lines = ",".join(matchline_list),
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
