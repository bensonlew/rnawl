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

class LncTargetTransAgent(Agent):
    """
    lnc target 预测
    """
    def __init__(self, parent):
        super(LncTargetTransAgent, self).__init__(parent)
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
            {'type': 'infile', 'name': 'g2t', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'lnc_gtf', 'format': 'lnc_rna.gtf'},
        ]
        self.add_option(options)

    def check_options(self):
        if self.option('gtf').is_set or self.option('g2t').is_set:
            pass
        else:
            raise OptionError("gtf g2t 需要至少设置一个")
        pass

    def set_resource(self):
        self._cpu = self.option("p") + 1
        self._memory = "{}G".format('60')

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
        super(LncTargetTransAgent, self).end()

class LncTargetTransTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(LncTargetTransTool, self).__init__(config)
        self.python_path = 'miniconda2/bin/python'
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

    def get_gene_list_summary(self):
        gene_list = []
        if self.option("diff_summary"):
            with open(self.option("diff_summary"), 'r') as diff_f:
                diff_f.readline()
                diff_f.readline()
                for line in diff_f:
                    cols = line.strip().split("\t")
                    gene_list.append(cols[0])
        return gene_list

    def get_t2g(self):
        m_dict = dict()
        if self.option('g2t').is_set:
            g2t = self.option('g2t').prop['path']
            with open(g2t, 'r') as g2t_f:
                for line in g2t_f.readlines():
                    cols = line.strip().split("\t")
                    m_dict[cols[1]] = cols[0]
        else:
            m_dict = self.option('gtf').get_txpt_gene_dic()
            if self.option('lnc_gtf').is_set:
                m_dict.update(self.option('lnc_gtf').get_txpt_gene_dic())
        return m_dict

    def add_annotation(self):
        t2g = self.get_t2g()
        t2n = self.get_trans2name()
        all_files = os.listdir(self.work_dir)
        for each in all_files:
            if each.endswith('merge_out'):
                fname = os.path.basename(each)
                with open(fname + ".annot.xls", 'w') as trans_fo, open(fname, 'r') as trans_f:
                    trans_fo.write("LncRNA id\tGene id\tTarget mRNA id\tTarget gene id\tTarget gene name\tlncRNA start\tlncRNA end\tmRNA start\tmRNA end\tEnergy\n")
                    trans_f.readline()
                    for line in trans_f:
                        cols = line.strip().split("\t")
                        trans_fo.write("\t".join([
                            cols[0],
                            t2g[cols[0]],
                            cols[1],
                            t2g[cols[1]],
                            t2n.get(cols[1], ""),
                            cols[2],
                            cols[3],
                            cols[4],
                            cols[5],
                            cols[6],
                        ]) + "\n")

    def get_trans2name(self):
        trans2name = dict()
        if self.option('annotation').is_set:
            with open(self.option('annotation').prop['path'], 'r') as annot_f:
                annot_f.readline()
                for line in annot_f:
                    cols = line.strip('').split("\t")
                    trans2name[cols[1]] = cols[3]
        return trans2name

    def get_translist_from_genelist(self, gene_list, t2g):
        tran_list = []
        for t,g in t2g.items():
            if g in gene_list:
                tran_list.append(t)
        return tran_list

    def choose_seq(self):
        '''
        根据列表筛选mRNA 与 靶基因序列
        '''
        t2g = self.get_t2g()
        choose_list = list()
        print "**** {}".format(self.query)
        if os.path.isfile(self.option('query_list')):
            with open(self.option('query_list'), 'r') as l_in:
                choose_list = [x.strip() for x in l_in.readlines()]

            if self.option('query_list_type') == 'G':
                choose_list = self.get_translist_from_genelist(choose_list, t2g)
            self.option('query').choose_seq_by_list(choose_list, 'query_choose.fa')
            self.query = 'query_choose.fa'
        elif self.option("diff_summary"):
            choose_list = self.get_gene_list_summary()
            print "choose_list is {}".format(choose_list)
            if self.option('query_list_type') == 'G':
                choose_list = self.get_translist_from_genelist(choose_list, t2g)
            self.option('query').choose_seq_by_list(choose_list, 'query_choose.fa')
            self.query = 'query_choose.fa'
        else:
            self.query = self.option('query').prop['path']

        print "**** {}".format(self.query)

        if os.path.isfile(self.option('target_list')):
            with open(self.option('target_list'), 'r') as l_in:
                choose_list = [x.strip() for x in l_in.readlines()]
            if self.option('target_list_type') == 'G':
                choose_list = self.get_translist_from_genelist(choose_list, t2g)
            self.option('target').choose_seq_by_list(choose_list, 'target_choose.fa')
            self.target = 'target_choose.fa'
        elif self.option("diff_summary"):
            choose_list = self.get_gene_list_summary()
            if self.option('target_list_type') == 'G':
                choose_list = self.get_translist_from_genelist(choose_list, t2g)
            self.option('target').choose_seq_by_list(choose_list, 'target_choose.fa')
            self.target = 'target_choose.fa'
        else:
            self.target = self.option('target').prop['path']


    def run_LncTargetTrans(self):
        cmd = '{} {} '.format(self.python_path, self.lncrna_target)
        cmd += '-q {} '.format(self.query)
        cmd += '-t {} '.format(self.target)
        cmd += '-parafly {} '.format(self.parafly)
        cmd += '-m {} '.format(self.option('method').lower())
        cmd += '-p {} '.format(self.option('p'))
        cmd += '-{} {} '.format(self.option('method').lower(),
                                self.software_dict[self.option('method').lower()])

        cmd_name = 'lncrnatarget'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
        else:
            self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

    def set_output(self):
        all_files = os.listdir(self.work_dir)
        for each in all_files:
            if each.endswith('merge_out.annot.xls'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

    def run(self):
        super(LncTargetTransTool, self).run()
        self.choose_seq()
        self.run_LncTargetTrans()
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
            "id": "LncTargetTrans_rna_plex" + datetime.datetime.now().strftime('%H-%M-%S'),
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
