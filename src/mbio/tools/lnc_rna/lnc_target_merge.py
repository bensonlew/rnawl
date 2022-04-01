# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import ConfigParser
import unittest
__author__ = 'liubinxu'

class LncTargetMergeAgent(Agent):
    """
    lnc target 预测
    """
    def __init__(self, parent):
        super(LncTargetMergeAgent, self).__init__(parent)
        options = [
            {'name': 'query', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'target', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'diff_summary', 'type': 'string', 'default': None},
            {'name': 'method', 'type': 'string'},
            {'name': 'p', 'type': 'int', "default": 30},
            {'name': 'query_list', 'type': 'string', "default": ''},
            {'name': 'query_list_type', 'type': 'string', "default": "G"},
            {'name': 'target_list', 'type': 'string', "default": ''},
            {'name': 'target_list_type', 'type': 'string', "default": "G"},
            {'type': 'infile', 'name': 'annotation', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'cis_known', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'cis_novol', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'exp_corr', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'g2t', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'novel_lnc_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'known_lnc_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'target_gtf', 'format': 'lnc_rna.gtf'},
        ]
        self.add_option(options)

    def check_options(self):
        '''
        if self.option('gtf').is_set or self.option('g2t').is_set:
            pass
        else:
            raise OptionError("gtf g2t 需要至少设置一个")
        '''
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('5')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(LncTargetMergeAgent, self).end()

class LncTargetMergeTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(LncTargetMergeTool, self).__init__(config)
        self.python_path = 'program/Python/bin/python'
        self.lncrna_target = self.config.PACKAGE_DIR + "/lnc_rna/lncrna_target.py"
        self.parafly = self.config.SOFTWARE_DIR + "/program/parafly-r2013-01-21/bin/bin/ParaFly"
        self.perl = '/program/perl-5.24.0/bin/'
        self.lnctar_lib = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/LncTar"
        self.boost = self.config.SOFTWARE_DIR + "/library/boost-1.69.0/lib"
        self.gcc_lib = self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64"
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.perl)
        self.set_environ(PERL5LIB=self.lnctar_lib)
        self.set_environ(LD_LIBRARY_PATH=self.boost)
        self.set_environ(LD_LIBRARY_PATH=self.gcc_lib)
        self.software_dict = dict(
            rnaplex = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/RNAplex",
            riblast = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/RIblast-1.1.3/RIblast",
            intarna = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/intaRNA-2.3.1/bin/IntaRNA",
            lnctar = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/LncTar/LncTar.pl",
        )
        self.query = ''
        self.target = ''

    def get_t2g(self):
        m_dict = dict()
        if self.option('g2t').is_set:
            g2t = self.option('g2t').prop['path']
            with open(g2t, 'r') as g2t_f:
                for line in g2t_f.readlines():
                    cols = line.strip().split("\t")
                    m_dict[cols[1]] = cols[0]
        else:
            m_dict = self.option('target_gtf').get_txpt_gene_dic()
            lnc_known_dict = self.option('known_lnc_gtf').get_txpt_gene_dic()
            lnc_novel_dict = self.option('novel_lnc_gtf').get_txpt_gene_dic()
            m_pos = self.option('target_gtf').get_txpt_pos_dic(type="G")
            lnc_known_pos = self.option('known_lnc_gtf').get_txpt_pos_dic()
            lnc_novel_pos= self.option('novel_lnc_gtf').get_txpt_pos_dic()

        return m_dict, lnc_known_dict, lnc_novel_dict, m_pos, lnc_known_pos, lnc_novel_pos

    def get_cis_target(self):
        target_dict = dict()
        if self.option("cis_known").is_set:
            with open(self.option("cis_known").prop['path']) as f:
                f.readline()
                for line in f:
                    cols = line.strip().split("\t")
                    target_dict[cols[0] + "|" + cols[3]] = [cols[-2], cols[-1]]
        if self.option("cis_novol").is_set:
            with open(self.option("cis_novol").prop['path']) as f:
                f.readline()
                for line in f:
                    cols = line.strip().split("\t")
                    target_dict[cols[0] + "|" + cols[3]] = [cols[-2], cols[-1]]
        return target_dict


    def add_annotation(self):
        t2g, lnc_k_t2g, lnc_n_t2g, m_pos, lnc_known_pos, lnc_novel_pos = self.get_t2g()
        t2n = self.get_trans2name()
        cis_dict = self.get_cis_target()

        with open(self.option("exp_corr").prop['path'], 'r') as corr_f, \
             open("cis_annot.xls", 'w') as cis_fo, \
             open("trans_annot.xls", 'w') as trans_fo:

            cis_fo.write("LncRNA id\tlnc_gene_id\tGene id\ttarget_type\tlnc_type\tTarget gene name\tcorr\tpvalue\tpadjust\tlocation\tdistance\tlnc_position\ttarget_position\tGene Description\n")
            trans_fo.write("LncRNA id\tlnc_gene_id\tGene id\ttarget_type\tlnc_type\tTarget gene name\tcorr\tpvalue\tpadjust\tGene Description\n")

            corr_f.readline()
            for line in corr_f:
                cols = line.strip().split("\t")
                if cols[0] in lnc_k_t2g:
                    lnc_type = "known"
                    loc_lnc = "{}:{}-{}".format(lnc_known_pos[cols[0]]["chr"],
                                                lnc_known_pos[cols[0]]["start"],
                                                lnc_known_pos[cols[0]]["end"])
                else:
                    lnc_type = "novel"
                    loc_lnc = "{}:{}-{}".format(lnc_novel_pos[cols[0]]["chr"],
                                                lnc_novel_pos[cols[0]]["start"],
                                                lnc_novel_pos[cols[0]]["end"])

                loc_m = "{}:{}-{}".format(m_pos[cols[1]]["chr"],
                                          m_pos[cols[1]]["start"],
                                          m_pos[cols[1]]["end"])

                if cols[0] in lnc_k_t2g:
                    lnc_gene = lnc_k_t2g[cols[0]]
                elif cols[0] in lnc_n_t2g:
                    lnc_gene = lnc_n_t2g[cols[0]]
                else:
                    lnc_gene = cols[0]

                if cols[0] + "|" + cols[1] in cis_dict:
                    cis_fo.write("\t".join([
                        cols[0],
                        lnc_gene,
                        cols[1],
                        "cis",
                        lnc_type,
                        t2n.get(cols[1], ["", ""])[0],
                        cols[2],
                        cols[3],
                        cols[4],
                    ] + cis_dict[cols[0] + "|" + cols[1]]  + [loc_lnc, loc_m]  + [t2n.get(cols[1], ["", ""])[1]]) + "\n")
                else:
                    trans_fo.write("\t".join([
                        cols[0],
                        lnc_gene,
                        cols[1],
                        "trans",
                        lnc_type,
                        t2n.get(cols[1], ["", ""])[0],
                        cols[2],
                        cols[3],
                        cols[4],
                        t2n.get(cols[1], ["", ""])[1]
                    ]) + "\n")

    def get_trans2name(self):
        trans2name = dict()
        if self.option('annotation').is_set:
            with open(self.option('annotation').prop['path'], 'r') as annot_f:
                annot_f.readline()
                for line in annot_f:
                    cols = line.strip('').split("\t")
                    trans2name[cols[0]] = [cols[3], cols[5]]
        return trans2name



    def set_output(self):
        all_files = os.listdir(self.work_dir)
        for each in all_files:
            if each.endswith('annot.xls'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

    def run(self):
        super(LncTargetMergeTool, self).run()
        self.add_annotation()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_lnc_rna'
        data = {
            "id": "LncTargetMerge_rna_plex" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "lnc_rna.lnc_target_trans",
            "instant": False,
            "options": dict(
                target = test_dir + '/test_data1/known_target.2.fa',
                query = test_dir + '/test_data1/known_linc.1.fa',
                method='rnaplex',
                diff_summary=test_dir + '/test_data1/diff_summary.xls'
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
