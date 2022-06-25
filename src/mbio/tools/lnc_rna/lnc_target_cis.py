# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'liubinxu'


class LncTargetCisAgent(Agent):
    """
    lnc rna target prediction
    """
    def __init__(self, parent):
        super(LncTargetCisAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'mrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'lncrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'mrna_bed', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'lncrna_bed', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'g2t', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'annotation', 'format': 'lnc_rna.common'},
            {'type': 'outfile', 'name': 'out', 'format': 'lnc_rna.common'},
            {'type': 'int', 'name': 'up_dis', 'default': 10},
            {'type': 'int', 'name': 'down_dis', 'default': 10},
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("mrna_gtf").is_set or self.option('mrna_bed').is_set:
            pass
        else:
            raise OptionError("mrna_gtf mrna_bed 需要至少设置一个")
        if self.option("lncrna_gtf").is_set or self.option('lncrna_bed').is_set:
            pass
        else:
            raise OptionError("lncrna_gtf lncrna_bed 需要至少设置一个")

    def set_resource(self):
        self._cpu = 1
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
        super(LncTargetCisAgent, self).end()


class LncTargetCisTool(Tool):
    """
    lnc rna target prediction
    """
    def __init__(self, config):
        super(LncTargetCisTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR

        self.python_path = 'miniconda2/bin/python'
        self.perl_path = 'miniconda2/bin/perl'
        self.bed2intersect_script = self.config.PACKAGE_DIR + '/lnc_rna/bed2intersect.py'
        self.gtf2bed_script = self.config.PACKAGE_DIR + '/lnc_rna/gtf2bed.pl'

    def gtf2bed(self, gtf_in, bed_out, cmd_name):
        cmd = '{} {} {} {}'.format(
            self.perl_path,
            self.gtf2bed_script,
            gtf_in,
            bed_out
        )
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
        return bed_out

    def get_t2g(self):
        m_dict = dict()
        if self.option('g2t').is_set:
            g2t = self.option('g2t').prop['path']
            with open(g2t, 'w') as g2t_f:
                for line in g2t_f.readlines():
                    cols = line.strip().split("\t")
                    m_dict[cols[1]] = cols[0]
        else:
            with open('g2t', 'w') as g2t_f:
                m_dict = self.option('mrna_gtf').get_txpt_gene_dic()
                lnc_dict = self.option('lncrna_gtf').get_txpt_gene_dic()
                m_dict.update(lnc_dict)
        return m_dict

    def add_annotation(self):
        t2g = self.get_t2g()
        t2n = self.get_trans2name()
        with open("lnc_rna_cistarget.xls", 'r') as cis_f, open("lnc_rna_cistarget.annot.xls", 'w') as cis_fo:
            cis_fo.write("LncRNA id\tGene id\tTarget mRNA id\tTarget gene id\tTarget gene name\tChromosome\tStrand\tLocation\tDistance\n")
            cis_f.readline()
            for line in cis_f:
                cols = line.strip().split("\t")
                cis_fo.write("\t".join([
                    cols[0],
                    t2g[cols[0]],
                    cols[1],
                    t2g[cols[1]],
                    t2n[cols[1]],
                    cols[2],
                    cols[6],
                    cols[7],
                    cols[8]
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


    def run_bed2intersect(self):
        if self.option('annotation').is_set:
            annot = self.option('annotation').prop['path']
        else:
            annot = 'None'
        if self.option('mrna_bed').is_set:
            mrna_bed = self.option('mrna_bed').prop['path']
        else:
            mrna_bed = self.gtf2bed(self.option('mrna_gtf').prop['path'], 'mrna_bed', 'mrna_gtf2bed')
        if self.option('lncrna_bed').is_set:
            lncrna_bed = self.option('lncrna_bed').prop['path']
        else:
            lncrna_bed = self.gtf2bed(self.option('lncrna_gtf').prop['path'], 'lncrna_bed', 'lncrna_gtf2bed')

        cmd = '{} {} '.format(self.python_path, self.bed2intersect_script)
        cmd += ' {} '.format(lncrna_bed)
        cmd += ' {} '.format(mrna_bed)
        cmd += ' {} '.format("lnc_rna_cistarget.xls")
        cmd += ' {} '.format(self.option("up_dis"))
        cmd += ' {} '.format(self.option("down_dis"))
        cmd_name = 'lnc_b2d_intersect'
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
        fname = "lnc_rna_cistarget.annot.xls"
        link = os.path.join(self.output_dir, fname)
        if os.path.exists(link):
            os.remove(link)
        os.link( self.work_dir + "/" +  fname, link)
        self.option("out", link)


    def run(self):
        super(LncTargetCisTool, self).run()
        self.run_bed2intersect()
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
        test_dir='/mnt/ilustre/users/isanger/sg-users/liubinxu/test_lnc_rna'
        data = {
            "id": "lnc_target_cis" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "lnc_rna.lnc_target_cis",
            "instant": False,
            "options": dict(
                lncrna_bed=test_dir + "/" + "lnc.tran.bed",
                mrna_bed=test_dir + "/" + "protein.tran.bed",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
